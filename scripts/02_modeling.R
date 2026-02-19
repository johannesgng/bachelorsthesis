# ==============================================================================
# File: 02_modeling.R
# Purpose: Perform the main statistical analysis:
#          1. Fourier Deseasonalization
#          2. Stationarity Tests
#          3. VAR Model Fitting (Pre-Whitening)
#          4. GARCH Model Fitting (Volatility)
#          5. Schaake Shuffle (Multivariate Reconstruction)
#          6. Forecast Reconstruction (Trend + Seasonality + Bias Correction)
#          7. Benchmark Generation (Climatology & Analog Ensemble)
# Output:  Saves forecast ensembles to 'output/models_and_ensembles.RData'
# ==============================================================================

# 1. Setup & Libraries ---------------------------------------------------------
library(dplyr)
library(vars)
library(rugarch)
library(lubridate)
library(zoo)
library(forecast)
library(tseries)
library(urca)

# Load helper functions (assumes running from project root 'Code_BA')
source("functions/fourier_utils.R")
source("functions/garch_utils.R")
source("functions/schaake_shuffle_utils.R")
source("functions/stationarity_utils.R")


# 2. Load Processed Data -------------------------------------------------------
# Loads: energy, energy_train, energy_test, threshold, split_date
load("output/data_processed_all.RData")

cat("Data loaded. Starting modeling pipeline...\n")

# 3. Fourier Decomposition (Seasonality & Trend) -------------------------------
series <- c("demand_log", "CF_solar_log", "CF_wind_log")
fourier_results <- list()

cat("Step 1/6: Optimizing Fourier Series orders (Grid Search)...\n")

for (s in series) {

  use_trend <- if(s == "demand_log") FALSE else TRUE 
  
  # Grid Search (AIC)
  aic_tab <- choose_best_P_combined(
    energy_train, col_name = s, 
    maxP_year = 15, maxP_month = 15, maxP_week = 15, 
    add_trend = use_trend
  )
  
  print(paste("Best Order for", s, ": Year=", aic_tab$P_year, "Month=", aic_tab$P_month, "Week=", aic_tab$P_week))
  
  # Fit Final Linear Model
  fourier_results[[s]] <- list(
    aic = aic_tab, 
    fitres = fit_and_deseason_combined(
      df = energy_train, 
      col_name = s, 
      P_year = aic_tab$P_year,   
      P_month = aic_tab$P_month,
      P_week = aic_tab$P_week,
      add_trend = use_trend
    )
  )
}

# Extract Deseasonalized Data for VAR
demand_dat <- fourier_results$demand_log$fitres$data
solar_dat  <- fourier_results$CF_solar_log$fitres$data
wind_dat   <- fourier_results$CF_wind_log$fitres$data

energy_train_deseasonalized <- data.frame(
  date = energy_train$date,
  demand_deseasonalized = demand_dat$deseasonalized,
  solar_deseasonalized  = solar_dat$deseasonalized,
  wind_deseasonalized   = wind_dat$deseasonalized
)

# 4. Stationarity Test ---------------------------------------------------------
cat("Stationarity Test... \n")

stationarity_check(energy_train_deseasonalized$demand_deseasonalized)
stationarity_check(energy_train_deseasonalized$solar_deseasonalized) 
stationarity_check(energy_train_deseasonalized$wind_deseasonalized)

johansen_test(energy_train_deseasonalized)

# 5. VAR Model (Pre-Whitening) -------------------------------------------------
cat("Step 2/6: Fitting VAR Model...\n")

# Lag Selection
lag_select <- VARselect(energy_train_deseasonalized[, -1], lag.max = 31, type = "const")
best_p_AIC <- as.integer(lag_select$selection["AIC(n)"])

# Fit VAR
model_var <- VAR(energy_train_deseasonalized[, -1], p = best_p_AIC, type = "const")

# Extract VAR Residuals for GARCH
var_residuals <- residuals(model_var)
res_demand <- var_residuals[, "demand_deseasonalized"]
res_solar  <- var_residuals[, "solar_deseasonalized"]
res_wind   <- var_residuals[, "wind_deseasonalized"]

# 6. GARCH Models (Volatility) -------------------------------------------------
cat("Step 3/6: Fitting GARCH Models (Grid Search). This takes time...\n")

# Grid Search & Fit (using the utility function)
# Note: scale_data=TRUE multiplies by 100 for numerical stability
uspec_demand <- univariate_garch_parameters(res_demand, scale_data = TRUE)
garch_demand <- ugarchfit(uspec_demand, data = res_demand * 100, solver = "hybrid")

uspec_solar <- univariate_garch_parameters(res_solar, scale_data = TRUE)
garch_solar <- ugarchfit(uspec_solar, data = res_solar * 100, solver = "hybrid")

uspec_wind <- univariate_garch_parameters(res_wind, scale_data = TRUE)
garch_wind <- ugarchfit(uspec_wind, data = res_wind * 100, solver = "hybrid")

# 7. Forecasting (VAR + GARCH + Schaake) ---------------------------------------
cat("Step 4/6: Generating Forecasts & Applying Schaake Shuffle...\n")

h <- nrow(energy_test)

# A) VAR Forecast (Conditional Mean)
fcst_var <- predict(model_var, n.ahead = h)
mu_demand <- fcst_var[["fcst"]][["demand_deseasonalized"]][,1]
mu_solar  <- fcst_var[["fcst"]][["solar_deseasonalized"]][,1]
mu_wind   <- fcst_var[["fcst"]][["wind_deseasonalized"]][,1]

# B) GARCH Forecast (Conditional Variance)
# Note: Need to divide sigma by 100 because we scaled inputs up
fcst_garch_d <- ugarchforecast(garch_demand, n.ahead = h)
sigma_d <- sigma(fcst_garch_d) / 100 

fcst_garch_s <- ugarchforecast(garch_solar, n.ahead = h)
sigma_s <- sigma(fcst_garch_s) / 100

fcst_garch_w <- ugarchforecast(garch_wind, n.ahead = h)
sigma_w <- sigma(fcst_garch_w) / 100

# C) Prepare Inputs for Schaake Shuffle
# Get standardized residuals from GARCH fit
e_std_d <- residuals(garch_demand, standardize = TRUE)
e_std_s <- residuals(garch_solar,  standardize = TRUE)
e_std_w <- residuals(garch_wind,   standardize = TRUE)

# Construct Rank Matrix (Empirical Copula)
error_std_mat <- cbind(as.numeric(e_std_d), as.numeric(e_std_s), as.numeric(e_std_w))
rank_mat <- get_rankmatrix(error_std_mat)

# Sample historical ranks (M = h, matching test set size for ensemble)
n_hist <- nrow(rank_mat)
M_ens  <- min(h, n_hist) 

set.seed(123) 
sample_idx <- sample(1:n_hist, size = M_ens, replace = TRUE)
rank_mat_use <- apply(rank_mat[sample_idx, ], 2, rank, ties.method = "first")

# D) Generate & Shuffle Ensembles
joint_ensembles <- vector("list", length = h)

for (k in seq_len(h)) {
  # 1. Generate Raw Ensemble (ECC-Q)
  raw_d <- r_empirical_dist(M_ens, mu_demand[k], sigma_d[k], as.numeric(e_std_d))
  raw_s <- r_empirical_dist(M_ens, mu_solar[k],  sigma_s[k], as.numeric(e_std_s))
  raw_w <- r_empirical_dist(M_ens, mu_wind[k],   sigma_w[k], as.numeric(e_std_w))
  
  raw_mat <- cbind(demand = raw_d, solar = raw_s, wind = raw_w)
  
  # 2. Apply Schaake Shuffle
  joint_ensembles[[k]] <- schaake_shuffle(raw_mat, rank_mat_use)
}
names(joint_ensembles) <- as.character(energy_test$date)

# 8. Reconstruction (Bias, Trend, Seasonality) ---------------------------------
cat("Step 5/6: Reconstructing Forecasts (Trend, Seasonality, Bias)...\n")

# Prepare Future Design Matrices
P_d <- fourier_results$demand_log$aic
P_s <- fourier_results$CF_solar_log$aic
P_w <- fourier_results$CF_wind_log$aic

X_fut_d <- future_X(energy_test$date, P_d$P_year, P_d$P_month, P_d$P_week)
X_fut_s <- future_X(energy_test$date, P_s$P_year, P_s$P_month, P_s$P_week)
X_fut_w <- future_X(energy_test$date, P_w$P_year, P_w$P_month, P_w$P_week)

lm_d <- fourier_results$demand_log$fitres$fit
lm_s <- fourier_results$CF_solar_log$fitres$fit
lm_w <- fourier_results$CF_wind_log$fitres$fit

# Calculate Deterministic Parts
det_d <- calc_annSeas_future(lm_d, X_fut_d) + add_trend_coef(lm_d, add_trend_future(lm_d, h))
det_s <- calc_annSeas_future(lm_s, X_fut_s) + add_trend_coef(lm_s, add_trend_future(lm_s, h))
det_w <- calc_annSeas_future(lm_w, X_fut_w) + add_trend_coef(lm_w, add_trend_future(lm_w, h))

# Adaptive Bias Correction (Last 60 Days)
calib_win <- 60
bias_d <- mean(tail(energy_train$demand_log, calib_win)) - mean(tail(predict(lm_d), calib_win))
bias_s <- mean(tail(energy_train$CF_solar_log, calib_win)) - mean(tail(predict(lm_s), calib_win))
bias_w <- mean(tail(energy_train$CF_wind_log, calib_win)) - mean(tail(predict(lm_w), calib_win))

det_d <- det_d + bias_d
det_s <- det_s + bias_s
det_w <- det_w + bias_w

# Combine Deterministic + Stochastic Parts & Transform back
ensembles_var_garch <- vector("list", length(joint_ensembles))
P_max_s_test <- energy_test$P_max_solar
P_max_w_test <- energy_test$P_max_wind

for (i in seq_along(joint_ensembles)) {
  # Add deterministic component
  log_d <- det_d[i] + joint_ensembles[[i]][, 1]
  log_s <- det_s[i] + joint_ensembles[[i]][, 2]
  log_w <- det_w[i] + joint_ensembles[[i]][, 3]
  
  # Exponentiate
  cf_s <- exp(log_s)
  cf_w <- exp(log_w)
  
  # Physical Cap (CF <= 1)
  cf_s[cf_s > 1] <- 1.0
  cf_w[cf_w > 1] <- 1.0
  
  # Convert to MWh
  ensembles_var_garch[[i]] <- cbind(
    demand = exp(log_d),
    solar  = cf_s * P_max_s_test[i],
    wind   = cf_w * P_max_w_test[i]
  )
}
names(ensembles_var_garch) <- as.character(energy_test$date)

# 9. Benchmarks ----------------------------------------------------------------
cat("Step 6/6: Generating Benchmarks...\n")

# A) Analog Ensemble (AnEn)
# Logic: +/- 14 days window (29 days total), matching wday, NO smoothing
historical_pool <- energy_train %>%
  dplyr::select(hist_date = date, hist_doy = doy, hist_wday = wday, 
                hist_demand = demand, hist_solar_cf = CF_solar, hist_wind_cf = CF_wind)

forecast_base <- data.frame(date = energy_test$date) %>%
  mutate(
    target_doy = ifelse(yday(date) == 366, 365, yday(date)),
    target_wday = wday(date, week_start = 1)
  )

p_max_future <- energy_test %>% dplyr::select(date, P_max_solar, P_max_wind)

ens_anen_df <- forecast_base %>%
  left_join(historical_pool, by = c("target_wday" = "hist_wday"), relationship = "many-to-many") %>%
  mutate(
    diff = abs(target_doy - hist_doy),
    dist = pmin(diff, 365 - diff)
  ) %>%
  filter(dist <= 14) %>%
  left_join(p_max_future, by = "date") %>%
  mutate(
    solar = pmin(hist_solar_cf * P_max_solar, P_max_solar),
    wind  = pmin(hist_wind_cf  * P_max_wind,  P_max_wind),
    demand = hist_demand
  ) %>%
  dplyr::select(date, demand, solar, wind)

ensembles_anen <- split(ens_anen_df, ens_anen_df$date)
ensembles_anen <- lapply(ensembles_anen, function(df) { df$date <- NULL; as.matrix(df) })

# B) Climatology Benchmark (Naive)
# Logic: Matching Calendar Week & Wday, NO smoothing
train_prep <- energy_train %>%
  dplyr::select(calendar_week, wday, doy, member_id = year, demand, solar_cf = CF_solar, wind_cf = CF_wind)

ens_naive_df <- forecast_base %>%
  mutate(calendar_week = ifelse(week(date) == 53, 52, week(date))) %>%
  # Join Demand (Week & Wday) -> HIER WAR DER FEHLER
  left_join(train_prep %>% dplyr::select(calendar_week, wday, member_id, demand), 
            by = c("calendar_week", "target_wday" = "wday"), relationship = "many-to-many") %>%
  # Join Solar/Wind (DOY & member_id to keep consistent weather scenarios)
  left_join(train_prep %>% dplyr::select(doy, member_id, solar_cf, wind_cf), 
            by = c("target_doy" = "doy", "member_id")) %>%
  left_join(p_max_future, by = "date") %>%
  mutate(
    solar = solar_cf * P_max_solar,
    wind  = wind_cf * P_max_wind
  ) %>%
  dplyr::select(date, demand, solar, wind)

ensembles_naive <- split(ens_naive_df, ens_naive_df$date)
ensembles_naive <- lapply(ensembles_naive, function(df) { df$date <- NULL; as.matrix(df) })

# 10. Save Results --------------------------------------------------------------
save(ensembles_var_garch, ensembles_anen, ensembles_naive, fourier_results,
     file = "output/models_and_ensembles.RData")

cat("Modeling complete. Results saved to 'output/models_and_ensembles.RData'.\n")