library(readr)
library(forecast)
library(vars)
library(tseries)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(RQuantLib)
library(lubridate)
library(urca)
library(rugarch)
library(FinTS)
library(RColorBrewer)
library(MultivCalibration)
library(scoringRules)

data_dir <- "data" # realtive path

energy_production <- read_delim(file.path(data_dir, "energy_production.csv"),
                                delim = ";", escape_double = F, col_types = cols(`Datum von` = col_date(format = "%Y-%m-%d"), 
                                                                                 `Datum bis` = col_date(format = "%Y-%m-%d")), 
                                locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                trim_ws = T)

energy_production$Wind_MWh <- energy_production$`Wind Onshore [MWh] Berechnete Auflösungen` + energy_production$`Wind Offshore [MWh] Berechnete Auflösungen`


energy_demand <- read_delim(file.path(data_dir, "energy_demand.csv"),
                            delim = ";", escape_double = F, col_types = cols(`Datum von` = col_date(format = "%d.%m.%Y"), 
                                                                             `Datum bis` = col_date(format = "%d.%m.%Y")), 
                            locale = locale(decimal_mark = ",", grouping_mark = "."), 
                            trim_ws = T)

installed_cap_daily <- read_delim(file.path(data_dir, "installed_capacity.csv"),
                                  delim = ";", escape_double = F, locale = locale(decimal_mark = ",", 
                                                                                  grouping_mark = "."), trim_ws = T)

threshold <- 0.06

energy <- data.frame(
  date = energy_production$`Datum von`,
  doy = yday(energy_production$`Datum von`), 
  month = month(energy_production$`Datum von`),
  year = year(energy_production$`Datum von`),
  wday = wday(energy_production$`Datum von`, week_start = 1),
  calendar_week = week(energy_production$`Datum von`),
  demand = energy_demand$`Netzlast [MWh] Berechnete Auflösungen`,
  solar = energy_production$`Photovoltaik [MWh] Berechnete Auflösungen`,
  wind = energy_production$Wind_MWh,
  residual_load = energy_demand$`Residuallast [MWh] Berechnete Auflösungen`,
  installed_solar = installed_cap_daily$`Photovoltaik [MW] Berechnete Auflösungen`,
  installed_wind = rowSums(cbind(installed_cap_daily$`Wind Offshore [MW] Berechnete Auflösungen`, 
                                 installed_cap_daily$`Wind Onshore [MW] Berechnete Auflösungen`), na.rm = T)
) %>% 
  as_tibble() %>%
  mutate(
    P_max_solar = installed_solar * 24,
    P_max_wind  = installed_wind * 24,
    P_max_total = P_max_solar + P_max_wind,
    CF_solar = solar / P_max_solar,
    CF_wind  = wind / P_max_wind,
    CF_total = (solar + wind) / P_max_total,
    demand_log   = log(demand),
    solar_log    = log(solar),
    wind_log     = log(wind),
    CF_solar_log = log(CF_solar),
    CF_wind_log  = log(CF_wind),
    h48_mean = rollmean(CF_total, k = 2, fill = NA, align = "right"),
    dark_doldrum = if_else(h48_mean < threshold, 1, 0)
  )

split_date <- as.Date("2022-12-31")
energy_train <- subset(energy, date <= split_date)
energy_test <- subset(energy, date > split_date)

# save data
save(energy, file = '/Users/johannesgoring/Desktop/Code_BA/00_Archiv/energy.RData')
save(energy_train, file = '/Users/johannesgoring/Desktop/Code_BA/00_Archiv/energy_train.RData')
save(energy_test, file = '/Users/johannesgoring/Desktop/Code_BA/00_Archiv/energy_test.RData')

## Deseasonalize with Fourier Series

# annual terms
build_fourier_annual <- function(doy, P = 3, period = 365) {
  n <- length(doy)
  if (P == 0) return(as.data.frame(matrix(nrow = n, ncol = 0)))
  X <- matrix(NA, nrow = n, ncol = 2 * P)
  colnames(X) <- unlist(lapply(1:P, function(p) c(paste0("C_a", p), paste0("S_a", p))))
  for (p in 1:P) {
    X[, 2 * p - 1] <- cos(2 * pi * p * doy / period)
    X[, 2 * p]     <- sin(2 * pi * p * doy / period)
  }
  as.data.frame(X)
}

# monthly terms
build_fourier_monthly <- function(month, P = 3, period = 12) {
  n <- length(month)
  if (P == 0) return(as.data.frame(matrix(nrow = n, ncol = 0)))
  X <- matrix(NA, nrow = n, ncol = 2 * P)
  colnames(X) <- unlist(lapply(1:P, function(p) c(paste0("C_m", p), paste0("S_m", p))))
  for (p in 1:P) {
    X[, 2 * p - 1] <- cos(2 * pi * p * month / period)
    X[, 2 * p]     <- sin(2 * pi * p * month / period)
  }
  as.data.frame(X)
}

# weekly terms
build_fourier_weekly <- function(wday, P = 3, period = 7) {
  n <- length(wday)
  if (P == 0) return(as.data.frame(matrix(nrow = n, ncol = 0)))
  X <- matrix(NA, nrow = n, ncol = 2 * P)
  colnames(X) <- unlist(lapply(1:P, function(p) c(paste0("C_w", p), paste0("S_w", p))))
  for (p in 1:P) {
    X[, 2 * p - 1] <- cos(2 * pi * p * wday / period) 
    X[, 2 * p] <- sin(2 * pi * p * wday / period) 
  }
  as.data.frame(X)
}

fit_fourier_lm_combined <- function(df, col, P_year = 6, P_month = 2, P_week = 4, add_trend = T) {
  
  # Simplification for leap years
  doy <- df$doy
  doy[doy == 366] <- 365
  
  # Design matrix
  Xy <- if(P_year > 0) build_fourier_annual(doy, P = P_year, period = 365) else as.data.frame(matrix(nrow=nrow(df), ncol=0))
  Xm <- if(P_month > 0) build_fourier_monthly(df$month, P = P_month, period = 12) else as.data.frame(matrix(nrow=nrow(df), ncol=0))
  Xw <- if(P_week > 0) build_fourier_weekly(df$wday, P = P_week, period = 7) else as.data.frame(matrix(nrow=nrow(df), ncol=0))
  X <- cbind(Xy, Xm, Xw)
  
  dat <- cbind(df, X)
  dat$.resp_col <- dat[[col]]
  
  rhs <- if(ncol(X) > 0) paste(colnames(X), collapse = " + ") else "1"
  if (add_trend) {
    dat$.tidx <- seq_len(nrow(dat))
    rhs <- paste(rhs, "+ .tidx")
  }
  form <- as.formula(paste(".resp_col ~", rhs))
  fit <- lm(formula = form, data = dat)
  list(fit = fit, data = dat, X = X)
}

choose_best_P_combined <- function(df, col_name, maxP_year = 8, maxP_month = 3, maxP_week = 3, add_trend = T, fallback = c(4,1)) {
  # empty container
  best <- list(AIC = Inf, P_year = NA, P_month = NA, P_week = NA)
  errors <- list()
  
  # grid search
  for(Py in 1:maxP_year){
    for(Pm in 0:maxP_month){
      for(Pw in 0:maxP_week){
        resP <- tryCatch({
          fit_fourier_lm_combined(df, col = col_name, P_year = Py, P_month = Pm, P_week = Pw, add_trend = add_trend)
        }, error = function(e) {
          errors[[paste0("Py",Py,"_Pm",Pm, "_Pw", Pw)]] <<- conditionMessage(e)
          NULL
        })
        if(is.null(resP)) next
        a <- AIC(resP$fit)
        if(a < best$AIC){ # using AIC to prevent overfitting of the seasonal-trend
          best$AIC <- a; best$P_year <- Py; best$P_month <- Pm; best$P_week <- Pw
        }
      }
    }
  }
  best$errors <- errors
  best
}

fit_and_deseason_combined <- function(df, col_name, P_year, P_month, P_week, add_trend = T) {
  
  res <- fit_fourier_lm_combined(df = df, col = col_name, P_year = P_year, P_month = P_month, P_week = P_week, add_trend = add_trend)
  fit <- res$fit
  dat <- res$data
  X <- res$X
  coef_all <- coef(fit)
  
  coef_X <- coef_all[colnames(X)]
  coef_X[is.na(coef_X)] <- 0
  annSeas_t <- as.matrix(X) %*% as.numeric(coef_X)
  
  dat$annSeas <- as.numeric(annSeas_t)
  dat$deseasonalized <- residuals(fit)
  dat$fitted_full <- predict(fit)
  list(fit = fit, data = dat, X = X)
}

series <- c("demand_log", "CF_solar_log", "CF_wind_log")
results <- list()

for (s in series) {
  
  
  trend <- if(s == "demand_log") FALSE else TRUE
  
  
  aic_tab <- choose_best_P_combined(energy_train, col_name = s, maxP_year = 15, maxP_month = 15, maxP_week = 15, add_trend = trend)
  
  # Extract best Fourier order for each time step
  bestP_year <- aic_tab$P_year
  bestP_month <- aic_tab$P_month
  bestP_week <- aic_tab$P_week
  
  # Results
  print(s)
  print(bestP_year)
  print(bestP_month)
  print(bestP_week)
  
  results[[s]] <- list(
    aic = aic_tab, 
    fitres = fit_and_deseason_combined(
      df = energy_train, 
      col_name = s, 
      P_year = bestP_year,   
      P_month = bestP_month,
      P_week = bestP_week,
      add_trend = trend
    )
  )
}

# demand
demand_dat <- results$demand$fitres$data
demand_fit <- results$demand$fitres$fit

# solar
solar_dat <- results$CF_solar_log$fitres$data
solar_fit <- results$CF_solar_log$fitres$fit

# wind
wind_dat <- results$CF_wind_log$fitres$data
wind_fit <- results$CF_wind_log$fitres$fit


# Data Frame with inputs for VAR-GARCH model
energy_train_deseasonalized <- data.frame(
  date = energy_train$date,
  demand_deseasonalized = demand_dat$deseasonalized,
  solar_deseasonalized = solar_dat$deseasonalized,
  wind_deseasonalized = wind_dat$deseasonalized
)


# Johansen test for cointegration
lag_select_levels_train <- VARselect(energy_train_deseasonalized[, -1], lag.max = 31, type = "const")
K_train <- lag_select_levels_train$selection["AIC(n)"]

johansen <- urca::ca.jo(energy_train_deseasonalized[, -1], type = "trace", spec = "longrun", K = K_train, ecdet = "const")

summary(johansen)

stationarity_check <- function(x){
  cat("\nADF (alt='stationary'):\n");  print(adf.test(x, alternative = "stationary"))
  cat("KPSS (null='Trend'):\n");      print(kpss.test(x, null = "Trend"))
}

stationarity_check(energy_train_deseasonalized$demand_deseasonalized)
stationarity_check(energy_train_deseasonalized$solar_deseasonalized) 
stationarity_check(energy_train_deseasonalized$wind_deseasonalized)

## Check, whether Fourier Deseasonalization worked
# Demand Log
(acf(energy_train$demand_log))
(pacf(energy_train$demand_log))
(acf(energy_train_deseasonalized$demand_deseasonalized))
(pacf(energy_train_deseasonalized$demand_deseasonalized))

# CF Solar Log
(acf(energy_train$CF_solar_log))
(pacf(energy_train$CF_solar_log))
(acf(energy_train_deseasonalized$solar_deseasonalized))
(pacf(energy_train_deseasonalized$solar_deseasonalized))

# CF Wind Log
(acf(energy_train$CF_wind_log))
(pacf(energy_train$CF_wind_log))
(acf(energy_train_deseasonalized$wind_deseasonalized))
(pacf(energy_train_deseasonalized$wind_deseasonalized))

## VAR
sel <- VARselect(energy_train_deseasonalized[,-1], lag.max = 40, type = "const")
sel$selection
(best_p_AIC <- as.integer(sel$selection["AIC(n)"]))

model_var_aic <- VAR(energy_train_deseasonalized[,-1], p = best_p_AIC, type = "const")

(model_roots_aic <- vars::roots(model_var_aic))
(max(model_roots_aic)) # < 1 -> VAR model is stationary

vars::serial.test(model_var_aic, lags.pt = best_p_AIC+3, type = "PT.asymptotic") # p-value > 0.05 -> VAR model captured correlation

vars::normality.test(model_var_aic) # suggests stylized facts

vars::arch.test(model_var_aic, lags.multi = 12) # implies heteroskedastic residuals -> GARCH models necessary

stability(model_var_aic)

## GARCH
var_res_aic <- residuals(model_var_aic)
res_demand_aic <- var_res_aic[, "demand_deseasonalized"]
res_solar_aic  <- var_res_aic[, "solar_deseasonalized"]
res_wind_aic   <- var_res_aic[, "wind_deseasonalized"]

univariate_garch_parameters <- function(res_data, scale_data = TRUE){
  
  if(scale_data) {
    res_data <- res_data * 100
  }
  
  # 2. Parameter Grid
  models_to_try <- c("sGARCH", "gjrGARCH", "eGARCH")
  orders_to_try <- list(c(1,1), c(1,2), c(2,1), c(2,2)) 
  distribution_to_try <- list("std", "sstd", "sged", "ged")
  
  results_df <- data.frame(
    model = character(),
    order = character(),
    dist = character(),
    aic = numeric(),
    arch_p = numeric(),
    lb_p = numeric(),
    nyblom_stable = logical(),
    stringsAsFactors = FALSE
  )
  
  # Storage list
  fitted_models <- list()
  
  total_combinations <- length(models_to_try) * length(orders_to_try) * length(distribution_to_try)
  cat("Starting grid search over", total_combinations, "combinations...\n")
  count <- 0
  
  for (m in models_to_try) {
    for (o in orders_to_try) {
      for(d in distribution_to_try){
        
        count <- count + 1
        model_id <- paste(m, paste(o, collapse=","), d, sep="_")
        
        # specification
        uspec_try <- try(ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
          variance.model = list(model = m, garchOrder = o), 
          distribution.model = d
        ), silent = TRUE)
        
        if (inherits(uspec_try, "try-error")) next
        
        fit_try <- try(ugarchfit(
          uspec_try, 
          data = res_data, 
          solver = "hybrid", 
          fit.control = list(scale = 1) 
        ), silent = TRUE)
        
        # Check Convergence
        if (inherits(fit_try, "try-error") || fit_try@fit[["convergence"]] != 0) {
          next
        }
        
        # AIC
        current_aic <- infocriteria(fit_try)["Akaike",]
        
        # Residual Diagnostics
        resid_std <- residuals(fit_try, standardize = TRUE)
        
        # ARCH Test
        arch_res <- try(ArchTest(resid_std, lags = 12), silent = TRUE) # Lags etwas reduzieren für Robustheit
        arch_p <- if (inherits(arch_res, "try-error")) 0 else arch_res$p.value
        
        # Ljung-Box
        lb_res <- try(Box.test(resid_std, lag = 12, type = "Ljung-Box"), silent = TRUE)
        lb_p <- if (inherits(lb_res, "try-error")) 0 else lb_res$p.value
        
        # Nyblom 
        nyblom_res <- try(nyblom(fit_try), silent = TRUE)
        nyblom_stable <- FALSE
        if (!inherits(nyblom_res, "try-error")) {
          stat <- nyblom_res$JointNyblomTest
          crit <- nyblom_res$JointCritical[2] # 5%
          if(length(stat) > 0 && length(crit) > 0 && stat < crit) {
            nyblom_stable <- TRUE
          }
        }
        
        # Results
        results_df[nrow(results_df) + 1, ] <- list(
          model = m,
          order = paste(o, collapse=","),
          dist = d,
          aic = current_aic,
          arch_p = arch_p,
          lb_p = lb_p,
          nyblom_stable = nyblom_stable
        )
        
        fitted_models[[model_id]] <- list(spec = uspec_try, fit = fit_try)
      }
    }
  }
  
  
  if (nrow(results_df) == 0) {
    cat("NO convergent models found.\n")
    return(NULL)
  }
  
  cat("\nSearch finished. Filtering results...\n")
  
  candidates <- results_df[results_df$arch_p > 0.05, ]
  
  if (nrow(candidates) == 0) {
    cat("Warning: No model fully eliminated ARCH effects (p > 0.05). Relaxing constraint to p > 0.01.\n")
    candidates <- results_df[results_df$arch_p > 0.01, ]
  }
  
  if (nrow(candidates) == 0) {
    cat("Warning: Severe ARCH effects remain in all models. Using all convergent models.\n")
    candidates <- results_df
  }
  
  candidates$score <- candidates$aic
  candidates$score[candidates$lb_p < 0.05] <- candidates$score[candidates$lb_p < 0.05] + 0.1
  
  # Best model
  best_row <- candidates[which.min(candidates$score), ]
  
  cat("\n--- Best Selected Model ---\n")
  cat("Model:", best_row$model, "| Order:", best_row$order, "| Dist:", best_row$dist, "\n")
  cat("AIC:", round(best_row$aic, 4), "| ARCH p-val:", round(best_row$arch_p, 4), "| LB p-val:", round(best_row$lb_p, 4), "\n")
  if(!best_row$nyblom_stable) cat("(Note: Nyblom stability test failed)\n")
  best_id <- paste(best_row$model, best_row$order, best_row$dist, sep="_")
  return(fitted_models[[best_id]]$spec)
}

uspec_demand_aic <- univariate_garch_parameters(res_demand_aic)
garch_demand_aic <- ugarchfit(uspec_demand_aic, data = res_demand_aic, solver = "hybrid")
uspec_solar_aic <- univariate_garch_parameters(res_solar_aic) 
garch_solar_aic <- ugarchfit(uspec_solar_aic, data = res_solar_aic, solver = "hybrid")
uspec_wind_aic <- univariate_garch_parameters(res_wind_aic)
garch_wind_aic <- ugarchfit(uspec_wind_aic, data = res_wind_aic, solver = "hybrid")


## Forecast
# Forecast of conditional variance & mean
h <- nrow(energy_test)
forecast_var_aic <- predict(model_var_aic, n.ahead = h)

# Demand
forecast_garch_demand_aic <- ugarchforecast(garch_demand_aic, n.ahead = h)
sigma_garch_demand_aic <- sigma(forecast_garch_demand_aic)
mu_demand_aic <- forecast_var_aic[["fcst"]][["demand_deseasonalized"]][,1]

# Solar
forecast_garch_solar_aic <- ugarchforecast(garch_solar_aic, n.ahead = h)
sigma_garch_solar_aic <- sigma(forecast_garch_solar_aic)
mu_solar_aic <- forecast_var_aic[["fcst"]][["solar_deseasonalized"]][,1]

# Wind
forecast_garch_wind_aic <- ugarchforecast(garch_wind_aic, n.ahead = h)
sigma_garch_wind_aic <- sigma(forecast_garch_wind_aic)
mu_wind_aic <- forecast_var_aic[["fcst"]][["wind_deseasonalized"]][,1]

# Demand
dist_garch_demand_aic <- garch_demand_aic@model$modeldesc$distribution
pars_demand <- coef(garch_demand_aic)
shape_demand <- pars_demand["shape"]
skew_demand  <- pars_demand["skew"]

# Solar 
dist_garch_solar_aic <- garch_solar_aic@model$modeldesc$distribution
pars_solar <- coef(garch_solar_aic)
shape_solar <- pars_solar["shape"]
skew_solar  <- pars_solar["skew"]

# Wind
dist_garch_wind_aic <- garch_wind_aic@model$modeldesc$distribution
pars_wind <- coef(garch_wind_aic)
shape_wind <- pars_wind["shape"]
skew_wind  <- pars_wind["skew"]

# Forecastind Distribution
forecast_dist_demand <- data.frame(
  horizon = 1:h,
  mu = mu_demand_aic,
  sigma = as.numeric(sigma_garch_demand_aic)
)

forecast_dist_solar <- data.frame(
  horizon = 1:h,
  mu = mu_solar_aic,
  sigma = as.numeric(sigma_garch_solar_aic)
)

forecast_dist_wind <- data.frame(
  horizon = 1:h,
  mu = mu_wind_aic,
  sigma = as.numeric(sigma_garch_wind_aic)
)

## Reconstruct multivariate Dependencies

# Rebuild Python function from GitHub Repo: https://github.com/FabianKaechele/Energy-Schaake/blob/main/EnergySchaake/functions_energy_schaake.py
# Standardized Residuals
e_demand_std <- residuals(garch_demand_aic, standardize = T)
e_solar_std  <- residuals(garch_solar_aic,  standardize = T)
e_wind_std   <- residuals(garch_wind_aic,   standardize = T)

error_std_mat <- cbind(
  demand = as.numeric(e_demand_std),
  solar  = as.numeric(e_solar_std),
  wind   = as.numeric(e_wind_std)
)

# Define Copula
get_rankmatrix <- function(data, param_dependence = F) { # default is the simpler empirical copola
  ranks <- apply(data, MARGIN = 2, rank, ties.method = "first") # Note: MARGIN = 2 => Funktion wird spaltenweise durchgeführt
  # Implied sorting of the values for each column and interchange value with index from the sorted order
  
  return(ranks) 
}

rank_mat <- get_rankmatrix(error_std_mat, param_dependence = F)

rank_mat <- apply(error_std_mat, MARGIN = 2, rank, ties.method = "first")

r_empirical_dist <- function(M, mu, sigma, z_hist) { # function 'density_forecast_nonparam' from functions_energy_schaake.py
  z_sorted <- sort(z_hist) # sorting the standardized residuals
  n_hist   <- length(z_sorted)
  
  # Quntiles
  probs <- (1:M) / (M + 1)
  
  z_quant <- approx(
    x = (1:n_hist) / (n_hist + 1),
    y = z_sorted,
    xout = probs,
    rule = 2
  )$y
  
  y_samples <- mu + sigma * z_quant
  
  return(as.numeric(y_samples))
}

# Shaake shuffle

schaake_shuffle <- function(forecast_mat, rank_mat) {
  stopifnot(all(dim(forecast_mat) == dim(rank_mat)))
  
  M <- nrow(forecast_mat)
  d <- ncol(forecast_mat)
  
  shuffled <- matrix(NA_real_, M, d) # initialize empty matrix
  colnames(shuffled) <- colnames(forecast_mat)
  
  # Rearrange data according to the rank-matrix
  for (j in seq_len(d)) {
    sorted_vals <- sort(forecast_mat[, j]) # sort ensemble values from small to large
    shuffled[, j] <- sorted_vals[rank_mat[, j]]
  }
  
  shuffled
}

n_hist <- nrow(rank_mat)
M_ens <- min(h, n_hist)

set.seed(123) 
sample_idx <- sample(1:n_hist, size = M_ens, replace = T)
rank_mat_use <- rank_mat[sample_idx, ]
rank_mat_use <- apply(rank_mat_use, MARGIN = 2, FUN = rank, ties.method = "first")

joint_ensembles <- vector("list", length = h)

for (k in seq_len(h)) {
  samp_demand_k <- r_empirical_dist(
    M = M_ens,
    mu = forecast_dist_demand$mu[k],
    sigma = forecast_dist_demand$sigma[k],
    z_hist = as.numeric(e_demand_std)
  )
  
  samp_solar_k <- r_empirical_dist(
    M = M_ens,
    mu = forecast_dist_solar$mu[k],
    sigma = forecast_dist_solar$sigma[k],
    z_hist = as.numeric(e_solar_std)
  )
  
  samp_wind_k <- r_empirical_dist(
    M = M_ens,
    mu = forecast_dist_wind$mu[k],
    sigma = forecast_dist_wind$sigma[k],
    z_hist = as.numeric(e_wind_std)
  )
  
  forecast_raw_k <- cbind(
    demand = samp_demand_k,
    solar = samp_solar_k,
    wind = samp_wind_k
  )
  
  joint_ensembles[[k]] <- schaake_shuffle(forecast_raw_k, rank_mat_use)
}

names(joint_ensembles) <- as.character(energy_test$date)

forecast_df <- data.frame(
  date = energy_test$date,
  step = 1:h
)

future_X <- function(dates, P_year, P_month, P_week){
  doy <- lubridate::yday(dates); doy[doy==366] <- 365
  Xy <- if(P_year>0) build_fourier_annual(doy, P=P_year, period=365) else as.data.frame(matrix(nrow=length(dates), ncol=0))
  Xm <- if(P_month>0) build_fourier_monthly(lubridate::month(dates), P=P_month, period=12) else as.data.frame(matrix(nrow=length(dates), ncol=0))
  Xw <- if(P_week>0) build_fourier_weekly(lubridate::wday(dates, week_start=1), P=P_week, period=7) else as.data.frame(matrix(nrow=length(dates), ncol=0))
  cbind(Xy, Xm, Xw)
}

calc_annSeas_future <- function(lmfit, Xfuture){
  coefs <- coef(lmfit)
  
  # Fourier Seasonality
  present <- intersect(names(coefs), colnames(Xfuture))
  seas <- if(length(present)==0) {
    rep(0, nrow(Xfuture))
  } else {
    as.numeric(as.matrix(Xfuture[, present, drop = FALSE]) %*% as.numeric(coefs[present]))
  }
  
  # Intercept
  intercept <- if("(Intercept)" %in% names(coefs)) coefs["(Intercept)"] else 0
  
  # Sum
  seas + intercept
}

add_trend_future <- function(lmfit, future_n){
  if(".tidx" %in% colnames(lmfit$model)){
    last_tidx <- max(lmfit$model$.tidx, na.rm = TRUE)
    return((last_tidx + seq_len(future_n)))
  } else return(rep(0, future_n))
}

add_trend_coef <- function(lmfit, trnd_vec){
  cfs <- coef(lmfit)
  if(".tidx" %in% names(cfs)){
    return(trnd_vec * cfs[".tidx"])
  } else return(rep(0, length(trnd_vec)))
}

# Fourier Orders are the same from training
P_demand <- results[["demand_log"]]$aic 
P_solar  <- results[["CF_solar_log"]]$aic
P_wind   <- results[["CF_wind_log"]]$aic

# Design Matrix
X_future_demand <- future_X(forecast_df$date, P_demand$P_year, P_demand$P_month, P_demand$P_week)
X_future_solar  <- future_X(forecast_df$date, P_solar$P_year,  P_solar$P_month,  P_solar$P_week)
X_future_wind   <- future_X(forecast_df$date, P_wind$P_year,   P_wind$P_month,   P_wind$P_week)

# OLS coefficients from training
demand_lm <- results[["demand_log"]]$fitres$fit
solar_lm  <- results[["CF_solar_log"]]$fitres$fit
wind_lm   <- results[["CF_wind_log"]]$fitres$fit

# Fourier series oos
forecast_df$demand_annSeas_future <- calc_annSeas_future(demand_lm, X_future_demand)
forecast_df$solar_annSeas_future  <- calc_annSeas_future(solar_lm,  X_future_solar)
forecast_df$wind_annSeas_future   <- calc_annSeas_future(wind_lm,   X_future_wind)

forecast_df$demand_trend_future <- add_trend_future(demand_lm, h)
forecast_df$solar_trend_future <- add_trend_future(solar_lm,  h)
forecast_df$wind_trend_future <- add_trend_future(wind_lm,   h)

demand_trend_future_coef <- add_trend_coef(demand_lm, forecast_df$demand_trend_future)
solar_trend_future_coef <- add_trend_coef(solar_lm,  forecast_df$solar_trend_future)
wind_trend_future_coef <- add_trend_coef(wind_lm,   forecast_df$wind_trend_future)
demand_det <- forecast_df$demand_annSeas_future + demand_trend_future_coef
solar_det <- forecast_df$solar_annSeas_future  + solar_trend_future_coef
wind_det <- forecast_df$wind_annSeas_future   + wind_trend_future_coef


## Adaptive Bias Correction
calib_window <- 60

# demand
obs_last_60_d  <- tail(energy_train$demand_log, calib_window)
pred_last_60_d <- tail(predict(demand_lm), calib_window)
(bias_d <- mean(obs_last_60_d) - mean(pred_last_60_d))

# solar
obs_last_60_s  <- tail(energy_train$CF_solar_log, calib_window)
pred_last_60_s <- tail(predict(solar_lm), calib_window)
(bias_s <- mean(obs_last_60_s) - mean(pred_last_60_s))

# wind
obs_last_60_w  <- tail(energy_train$CF_wind_log, calib_window)
pred_last_60_w <- tail(predict(wind_lm), calib_window)
(bias_w <- mean(obs_last_60_w) - mean(pred_last_60_w))

demand_det <- demand_det + bias_d
solar_det  <- solar_det  + bias_s
wind_det   <- wind_det   + bias_w

## Generate Ensembles
ensembles_var_garch <- vector("list", length(joint_ensembles))

P_max_solar_test <- energy_test$P_max_solar
P_max_wind_test  <- energy_test$P_max_wind

for (i in seq_along(joint_ensembles)) {
  demand_log <- demand_det[i] + joint_ensembles[[i]][, 1]
  solar_log <- solar_det[i] + joint_ensembles[[i]][, 2]
  wind_log  <- wind_det[i] + joint_ensembles[[i]][, 3]
  
  solar_cf <- exp(solar_log)
  wind_cf <- exp(wind_log)
  
  # Cut off values greater 1
  solar_cf[solar_cf > 1] <- 1.0
  wind_cf[wind_cf > 1]   <- 1.0
  
  ensembles_var_garch[[i]] <- cbind(
    demand = exp(demand_log),
    solar = solar_cf * P_max_solar_test[i],
    wind = wind_cf * P_max_wind_test[i]
  )
}

names(ensembles_var_garch) <- as.character(forecast_df$date)

## Benchmark Analog Ensemble

historical_pool <- energy_train %>%
  dplyr::select(
    hist_date = date,
    hist_year = year,
    hist_doy = doy,
    hist_wday = wday,
    hist_demand = demand,
    hist_solar_cf = CF_solar,
    hist_wind_cf = CF_wind
  )

forecast_dates <- energy_test$date 

forecast_base <- data.frame(date = forecast_dates) %>%
  mutate(
    target_doy = yday(date),
    target_doy = ifelse(target_doy == 366, 365, target_doy), # leap year
    target_wday = wday(date, week_start = 1)
  )

p_max_future <- energy %>% dplyr::select(date, P_max_solar, P_max_wind)

ensemble_forecast <- forecast_base %>%
  
  dplyr::left_join(historical_pool, by = c("target_wday" = "hist_wday"), relationship = "many-to-many") %>% 
  dplyr::mutate(
    diff = abs(target_doy - hist_doy),
    dist = pmin(diff, 365 - diff)
  ) %>%
  dplyr::filter(dist <= 14) %>% # window size 14
  
  dplyr::left_join(p_max_future, by = "date") %>%
  
  dplyr::mutate(
    demand = hist_demand,
    solar = pmin(hist_solar_cf * P_max_solar, P_max_solar),
    wind  = pmin(hist_wind_cf  * P_max_wind,  P_max_wind)
  ) %>%
  
  dplyr::select(date, demand, solar, wind)

ensembles_anen <- split(ensemble_forecast, ensemble_forecast$date)
ensembles_anen <- lapply(ensembles_anen, function(df) {
  df$date <- NULL
  return(as.matrix(df))
})

## Smoothed Climatology Ensemble Benchmark (before Rolling Mean)
p_max_future <- energy %>% dplyr::select(date, P_max_solar, P_max_wind)

train_prepared <- energy_train

# wday und calendar_week for demand
demand_scenarios <- train_prepared %>% dplyr::select(calendar_week, wday, member_id = year, demand) # week day information and calendar week (or doy) information is crucial for the energy demand as well

# doy for weather related energy data
weather_scenarios <- train_prepared %>% dplyr::select(doy, member_id = year, solar = CF_solar, wind = CF_wind) # since the produced energy is mainly dependet on meterological events, doy is a sufficient fiter criterion

# Ensemble forecast
target_years <- 2023:2025
date_seq <- seq(as.Date("2023-01-01"), as.Date("2025-09-30"), by = "day")

forecast_base <- data.frame(date = date_seq) %>%
  mutate(
    raw_doy = yday(date),
    doy = ifelse(raw_doy == 366, 365, raw_doy), # leap years
    
    raw_week = week(date),
    calendar_week = ifelse(raw_week == 53, 52, raw_week), # leap years
    
    wday = wday(date, week_start = 1)
  )

ensemble_forecast <- forecast_base %>%
  left_join(demand_scenarios, by = c("calendar_week", "wday"), relationship = "many-to-many") %>%
  left_join(weather_scenarios, by = c("doy", "member_id")) %>%
  left_join(p_max_future, by = "date") %>%
  mutate(
    solar = solar * P_max_solar,
    wind = wind * P_max_wind
  )

ensemble_forecast <- ensemble_forecast %>% dplyr::select(date, demand, solar, wind)

ensembles_naive <- split(ensemble_forecast, ensemble_forecast$date)
ensembles_naive <- lapply(ensembles_naive, function(df) {
  df$date <- NULL
  return(as.matrix(df))
})

## Evaluation

### Evaluation via ES & CRPS
scoring <- function(ensembels, data){
  
  results_df <- data.frame(
    date = data$date,
    es_trivar = NA,
    es_bivar = NA,
    crps_demand = NA,
    crps_solar = NA,
    crps_wind = NA
  )
  
  n_days <- nrow(data)
  
  for(i in 1:n_days){
    
    # Trivariate Energy Score
    ens_mat <- ensembels[[i]]
    
    obs_demand <- data$demand[i]
    obs_solar <- data$solar[i]
    obs_wind <- data$wind[i]
    
    obs_vec <- c(obs_demand, obs_solar, obs_wind)
    
    results_df$es_trivar[i] <- es_sample(y = obs_vec, dat = t(ens_mat))
    
    # Bivariate Energy Score
    ens_mat_bivar <- ensembels[[i]][,c("solar", "wind")]
    
    obs_vec_bivar <- c(obs_solar, obs_wind)
    
    results_df$es_bivar[i] <- es_sample(y = obs_vec_bivar, dat = t(ens_mat_bivar))
    
    # CRPS
    results_df$crps_demand[i] <- crps_sample(y = obs_demand, dat = as.numeric(ens_mat[,1, drop = F]))
    results_df$crps_solar[i] <- crps_sample(y = obs_solar, dat = as.numeric(ens_mat[,2, drop = F]))
    results_df$crps_wind[i] <- crps_sample(y = obs_wind, dat = as.numeric(ens_mat[,3, drop = F]))
  }
  
  results_df <- results_df %>% mutate(es_trivar_mean = mean(es_trivar),
                                      es_bivar_mean = mean(es_bivar),
                                      crps_demand_mean = mean(crps_demand),
                                      crps_solar_mean = mean(crps_solar),
                                      crps_wind_mean = mean(crps_wind))
  
  return(results_df)
}

scoring_var_garch <- scoring(ensembles_var_garch, energy_test)
scoring_anen <- scoring(ensembles_anen, energy_test)
scoring_naive <- scoring(ensembles_naive, energy_test)

## Diebold Mariano Test for comparison
# Trivariate
forecast::dm.test(scoring_var_garch$es_trivar, scoring_anen$es_trivar, alternative = "less", h = 1, power = 1)
forecast::dm.test(scoring_var_garch$es_trivar, scoring_naive$es_trivar, alternative = "less", h = 1, power = 1)
forecast::dm.test(scoring_anen$es_trivar, scoring_naive$es_trivar, alternative = "less", h = 1, power = 1)

# Bivariate
forecast::dm.test(scoring_var_garch$es_bivar, scoring_anen$es_bivar, alternative = "less", h = 1, power = 1)
forecast::dm.test(scoring_var_garch$es_bivar, scoring_naive$es_bivar, alternative = "less", h = 1, power = 1)
forecast::dm.test(scoring_anen$es_bivar, scoring_naive$es_bivar, alternative = "less", h = 1, power = 1)

## Dark Doldrum Evaluation
comparison_df <- scoring_var_garch %>%
  dplyr::select(date, es_bivar_schaake = es_bivar) %>%
  mutate(
    es_bivar_anen = scoring_anen$es_bivar,
    es_bivar_naive = scoring_naive$es_bivar
  ) %>%
  inner_join(energy_test, by = "date")

scoring_summary <- comparison_df %>%
  group_by(dark_doldrum) %>%
  summarise(
    n_days = n(),
    mean_ES_Var_Garch = mean(es_bivar_schaake, na.rm = T),
    mean_ES_AnEn = mean(es_bivar_anen, na.rm = T),
    mean_ES_Rolling = mean(es_bivar_naive, na.rm = T),
    ESS_Var_Garch_vs_AnEn_Pct = (mean_ES_AnEn - mean_ES_Var_Garch) / mean_ES_AnEn * 100,
    ESS_Var_Garch_vs_naive_Pct = (mean_ES_Rolling - mean_ES_Var_Garch) / mean_ES_Rolling * 100
  )

print(scoring_summary)

## Brier Score
h48_mean_scatter <- function(ensembles, obs_data, threshold_cf = threshold) {
  
  results_df <- data.frame(
    date = as.Date(names(ensembles)),
    prob = NA_real_,
    obs_event = NA_real_,
    M = NA_integer_
  )
  
  dates <- as.Date(names(ensembles))
  
  for (t in seq_along(dates)) {
    if (t == 1) next
    
    obs_prev <- obs_data %>% filter(date == dates[t-1])
    obs_curr <- obs_data %>% filter(date == dates[t])
    if(nrow(obs_prev) == 0 || nrow(obs_curr) == 0) next
    
    ens_mat = ensembles[[t]]
    
    M_current = nrow(ens_mat)
    
    df_curr <- as.data.frame(ens_mat)
    
    df_48h <- data.frame(
      solar_48h = (obs_prev$solar + df_curr$solar) / 2,
      wind_48h  = (obs_prev$wind + df_curr$wind) / 2
    )
    
    limit_48h <- threshold_cf * obs_curr$P_max_total
    df_48h$is_doldrum <- (df_48h$solar_48h + df_48h$wind_48h) < limit_48h
    
    obs_48h_val <- (obs_prev$solar + obs_curr$solar) / 2 + (obs_prev$wind + obs_curr$wind) / 2
    is_obs_doldrum <- as.numeric(obs_48h_val < limit_48h)
    
    results_df$prob[t] <- mean(df_48h$is_doldrum) # f_t
    results_df$obs_event[t] <- is_obs_doldrum # o_t   
    results_df$M[t] <- M_current
    
    obs_plot_df <- data.frame(
      solar_48h = (obs_prev$solar + obs_curr$solar) / 2,
      wind_48h  = (obs_prev$wind + obs_curr$wind) / 2
    )
  }
  
  results_df <- results_df %>% filter(!is.na(prob))
  return(results_df)
}

brier_df_an_en <- h48_mean_scatter(ensembles_anen, energy_test)
brier_df_naive <- h48_mean_scatter(ensembles_naive, energy_test)
brier_df_var_garch <- h48_mean_scatter(ensembles_var_garch, energy_test)

brier_score <- function(df){ # adjusted Brier-Score
  result_bs <- (df$prob - df$obs_event)^2
  print(mean(result_bs))
  return(mean(result_bs))
}

bs_anen <- brier_score(brier_df_an_en)
bs_naive <- brier_score(brier_df_naive)
bs_var_garch <- brier_score(brier_df_var_garch)

brier_df <- data.frame(
  AnEn = bs_anen,
  Climatology_Benchmark = bs_naive,
  Var_Garch = bs_var_garch,
  BSS_Var_Garch_vs_AnEn_Pct = ((-1) * ((bs_var_garch - bs_anen) / bs_anen))*100, 
  BSS_Var_Garch_vs_Climatology_Benchmark_Pct = ((-1) * ((bs_var_garch - bs_naive) / bs_naive))*100
)

rownames(brier_df) <- "Brier_Score"
brier_df

heatmap_df <- tibble(
  date = brier_df_an_en$date,
  prob_anen = brier_df_an_en$prob,
  prob_naive = brier_df_naive$prob,
  prob_var_garch = brier_df_var_garch$prob,
  doy = wday(brier_df_an_en$date),
  calendar_week = week(brier_df_an_en$date),
  month = month(brier_df_an_en$date),
  year = year(brier_df_an_en$date)
)

hm_all <- heatmap_df %>%
  mutate(
    month = factor(month, levels = 1:12, labels = month.abb),
    doy = factor(doy, levels = 1:7,
                 labels = c("Mo", "Di", "Mi", "Do", "Fr", "Sa", "So"))
  ) %>%
  dplyr::select(month, doy, prob_anen, prob_var_garch, prob_naive) %>%
  pivot_longer(
    starts_with("prob_"),
    names_to = "model",
    values_to = "prob"
  ) %>%
  group_by(model, month, doy) %>%
  summarise(prob = mean(prob, na.rm = TRUE), .groups = "drop")

ggplot(hm_all, aes(x = month, y = doy, fill = prob)) +
  geom_tile() +
  facet_wrap(~ model) +
  scale_fill_viridis_c(labels = scales::percent) +
  labs(
    x = "Month",
    y = "Week Day",
    fill = "Average Probability",
    title = "Dark Doldrum Probability per Week Day and Month"
  ) +
  theme_minimal()

energy_train <- energy_train %>%
  group_by(wday) %>%
  mutate(demand_quartile = ntile(demand, 4)) %>%
  ungroup()

quantiles_by_wday <- energy_train %>%
  group_by(wday) %>%
  summarise(
    q1 = quantile(demand, 0.25),
    q2 = quantile(demand, 0.50),
    q3 = quantile(demand, 0.75)
  )

energy_test <- energy_test %>%
  left_join(quantiles_by_wday, by = "wday") %>%
  mutate(
    demand_quartile = case_when(
      demand <= q1 ~ 1,
      demand <= q2 ~ 2,
      demand <= q3 ~ 3,
      T ~ 4
    )
  ) %>%
  dplyr::select(-q1, -q2, -q3)

scoring_var_garch <- as_tibble(scoring_var_garch)
scoring_anen <- as_tibble(scoring_anen)
scoring_naive <- as_tibble(scoring_naive)

comparison_var_garch <- scoring_var_garch %>%
  inner_join(energy_test %>% dplyr::select(date, dark_doldrum, demand_quartile), by = "date") %>%
  group_by(dark_doldrum, demand_quartile) %>%
  summarise(
    n = n(),
    es_trivar_mean = mean(es_trivar, na.rm = T),
    es_bivar_mean = mean(es_bivar, na.rm = T),
    crps_demand_mean = mean(crps_demand, na.rm = T),
    crps_solar_mean = mean(crps_solar, na.rm = T),
    crps_wind_mean = mean(crps_wind, na.rm = T),
    .groups = "drop"
  )

comparison_anen <- scoring_anen %>%
  inner_join(energy_test %>% dplyr::select(date, dark_doldrum, demand_quartile), by = "date") %>%
  group_by(dark_doldrum, demand_quartile) %>%
  summarise(
    n = n(),
    es_trivar_mean = mean(es_trivar, na.rm = T),
    es_bivar_mean = mean(es_bivar, na.rm = T),
    crps_demand_mean = mean(crps_demand, na.rm = T),
    crps_solar_mean = mean(crps_solar, na.rm = T),
    crps_wind_mean = mean(crps_wind, na.rm = T),
    .groups = "drop"
  )

comparison_naive <- scoring_naive %>%
  inner_join(energy_test %>% dplyr::select(date, dark_doldrum, demand_quartile), by = "date") %>%
  group_by(dark_doldrum, demand_quartile) %>%
  summarise(
    n = n(),
    es_trivar_mean = mean(es_trivar, na.rm = TRUE),
    es_bivar_mean = mean(es_bivar, na.rm = TRUE),
    crps_demand_mean = mean(crps_demand, na.rm = TRUE),
    crps_solar_mean = mean(crps_solar, na.rm = TRUE),
    crps_wind_mean = mean(crps_wind, na.rm = TRUE),
    .groups = "drop"
  )

improvement_df <- data.frame(
  dark_doldrum = comparison_var_garch$dark_doldrum,
  demand_quartile = comparison_var_garch$demand_quartile,
  N_Obs = comparison_var_garch$n,
  # Trivar(
  es_trivar_mean_var_garch_vs_anen_pct = ((comparison_anen$es_trivar_mean - comparison_var_garch$es_trivar_mean) / comparison_anen$es_trivar_mean) * 100,
  es_trivar_mean_var_garch_vs_naive_pct = ((comparison_naive$es_trivar_mean - comparison_var_garch$es_trivar_mean) / comparison_naive$es_trivar_mean) * 100,
  # Bivar
  es_bivar_mean_var_garch_vs_anen_pct = ((comparison_anen$es_bivar_mean - comparison_var_garch$es_bivar_mean) / comparison_anen$es_bivar_mean) * 100,
  es_bivar_mean_var_garch_vs_naive_pct = ((comparison_naive$es_bivar_mean - comparison_var_garch$es_bivar_mean) / comparison_naive$es_bivar_mean) * 100,
  # CRPS vs ANEN
  crps_demand_mean_var_garch_vs_anen_pct = ((comparison_anen$crps_demand_mean - comparison_var_garch$crps_demand_mean) / comparison_anen$crps_demand_mean) * 100,
  crps_solar_mean_var_garch_vs_anen_pct = ((comparison_anen$crps_solar_mean - comparison_var_garch$crps_solar_mean) / comparison_anen$crps_solar_mean) * 100,
  crps_wind_mean_var_garch_vs_anen_pct = ((comparison_anen$crps_wind_mean - comparison_var_garch$crps_wind_mean) / comparison_anen$crps_wind_mean) * 100,
  # CRPS vs Climatology_Benchmark
  crps_demand_mean_var_garch_vs_naive_pct = ((comparison_naive$crps_demand_mean - comparison_var_garch$crps_demand_mean) / comparison_naive$crps_demand_mean) * 100,
  crps_solar_mean_var_garch_vs_naive_pct = ((comparison_naive$crps_solar_mean - comparison_var_garch$crps_solar_mean) / comparison_naive$crps_solar_mean) * 100,
  crps_wind_mean_var_garch_vs_naive_pct = ((comparison_naive$crps_wind_mean - comparison_var_garch$crps_wind_mean) / comparison_naive$crps_wind_mean) * 100
)

print(improvement_df)

## Plot
impr_df_df <- improvement_df %>% filter(dark_doldrum == 1)

plot_data <- impr_df_df %>% 
  dplyr::select(demand_quartile,
                es_trivar_mean_var_garch_vs_anen_pct,
                es_trivar_mean_var_garch_vs_naive_pct) %>%
  pivot_longer(
    cols = -demand_quartile,
    names_to = "comparison",
    values_to = "improvement_pct"
  ) %>%
  mutate(
    comparison = case_when(
      comparison == "es_trivar_mean_var_garch_vs_anen_pct" ~ "VAR-GARCH vs AnEn",
      comparison == "es_trivar_mean_var_garch_vs_naive_pct" ~ "VAR-GARCH vs Climatology_Benchmark"
    )
  )

p <- ggplot(plot_data,
            aes(x = factor(demand_quartile),
                y = improvement_pct,
                fill = comparison)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "VAR-GARCH Performance vs. Benchmarks during Dark Doldrums (ES trivar)",
    subtitle = "",
    x = "Demand Quartile",
    y = "Improvement in %",
    fill = ""
  ) +
  theme_minimal()

print(p)

### PIT
plot_pit_univariate <- function(ensemble_list, obs_df, title_suffix = "") {
  n_days <- length(ensemble_list)
  ranks_df <- data.frame(
    date = obs_df$date[1:n_days],
    rank_demand = NA,
    rank_solar   = NA,
    rank_wind    = NA
  )
  
  get_pit <- function(obs, ens_vec) {
    M <- length(ens_vec)
    r <- sum(ens_vec <= obs, na.rm = TRUE) + 1
    return((r - 0.5) / (M + 1)) # Normalized PIT
  }
  
  for (i in 1:n_days) {
    ens_mat <- ensemble_list[[i]]
    if (is.null(ens_mat) || nrow(ens_mat) == 0) next
    
    ranks_df$rank_demand[i] <- get_pit(obs_df$demand[i], ens_mat[, "demand"])
    ranks_df$rank_solar[i]  <- get_pit(obs_df$solar[i],  ens_mat[, "solar"])
    ranks_df$rank_wind[i]   <- get_pit(obs_df$wind[i],   ens_mat[, "wind"])
  }
  
  n_bins <- 10 
  breaks_vec <- seq(0, 1, length.out = n_bins + 1)
  
  kit_green <- "#009682"
  palette_vec <- colorRampPalette(c("gray85", kit_green))(100)
  
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 4, 1))
  vars   <- c("rank_demand", "rank_solar", "rank_wind")
  titles <- c("Energy Demand", "Solar Power", "Wind Power")
  
  for (i in 1:3) {
    data_vec <- ranks_df[[vars[i]]]
    data_vec <- data_vec[!is.na(data_vec)]
    
    h <- hist(data_vec, breaks = breaks_vec, plot = FALSE)
    
    col_idx <- if (max(h$counts) > 0)
      floor((h$counts / max(h$counts)) * 99) + 1 else 1
    
    plot(h,
         col = palette_vec[col_idx],
         border = "white",
         main = titles[i],
         xlab = "Normalized Rank",
         ylab = "Density",
         freq = FALSE,
         xlim = c(0, 1),
         las = 1)
    
    abline(h = 1, col = "darkgrey", lwd = 2, lty = 2)
    
    box(bty = "l", col = "gray40")
  }
  
  mtext(
    paste("Univariate PIT Histogram", title_suffix),
    outer = TRUE, side = 3, line = 2, font = 2, cex = 1.3
  )
  
  par(mfrow = c(1, 1))
}
plot_pit_multivariate <- function(ensemble_list, obs_df, title_suffix = "") {
  n_days <- length(ensemble_list)
  multivar_pit <- numeric(n_days)
  
  for (i in 1:n_days) {
    ens_mat <- ensemble_list[[i]]
    if (is.null(ens_mat) || nrow(ens_mat) == 0) {
      multivar_pit[i] <- NA
      next
    }
    
    obs_vec <- c(obs_df$demand[i], obs_df$solar[i], obs_df$wind[i])
    M <- nrow(ens_mat)
    
    m_rank <- MultivCalibration::get_prerank(
      y = obs_vec,
      x = t(ens_mat),
      prerank = "average_rank",
      return_rank = TRUE
    )
    
    # Normalize
    multivar_pit[i] <- (m_rank - 0.5) / (M + 1)
  }
  
  multivar_pit <- multivar_pit[!is.na(multivar_pit)]
  
  n_bins <- 10
  breaks_vec <- seq(0, 1, length.out = n_bins + 1)
  
  kit_green <- "#009682"
  palette_vec <- colorRampPalette(c("gray85", kit_green))(100)
  
  par(mar = c(4, 4, 4, 2))
  
  h <- hist(multivar_pit, breaks = breaks_vec, plot = FALSE)
  col_idx <- if (max(h$counts) > 0)
    floor((h$counts / max(h$counts)) * 99) + 1 else 1
  
  plot(h,
       col = palette_vec[col_idx],
       border = "white",
       freq = FALSE, 
       main = paste("Multivariate PIT Histogram", title_suffix),
       xlab = "Normalized Rank",
       ylab = "Density",
       xlim = c(0, 1),
       las = 1)
  
  abline(h = 1, col = "darkgrey", lwd = 2, lty = 2)
  
  box(bty = "l", col = "gray40")
}

# VAR-GARCH (Schaake)
plot_pit_univariate(ensembles_var_garch, energy_test, "- VAR-GARCH")
plot_pit_multivariate(ensembles_var_garch, energy_test, "- VAR-GARCH")

# Analog Ensemble
plot_pit_univariate(ensembles_anen, energy_test, "- AnEn")
plot_pit_multivariate(ensembles_anen, energy_test, "- AnEn")

# Rolling Mean
plot_pit_univariate(ensembles_naive, energy_test, "- Climatology Benchmark")
plot_pit_multivariate(ensembles_naive, energy_test, "- Climatology Benchmark")

## Residuallast-Ensembles
residual_load_var_garch <- list()

for(i in 1:length(ensembles_var_garch)){
  residual_load_var_garch[[i]] <- ensembles_var_garch[[i]][,1] - ensembles_var_garch[[i]][,2] - ensembles_var_garch[[i]][,3]
}

names(residual_load_var_garch) <- names(ensembles_var_garch)

residual_load_anen <- list()

for(i in 1:length(ensembles_anen)){
  residual_load_anen[[i]] <- ensembles_anen[[i]][,1] - ensembles_anen[[i]][,2] - ensembles_anen[[i]][,3]
}

names(residual_load_anen) <- names(ensembles_anen)

residual_load_naive <- list()

for(i in 1:length(ensembles_naive)){
  residual_load_naive[[i]] <- ensembles_naive[[i]][,1] - ensembles_naive[[i]][,2] - ensembles_naive[[i]][,3]
}

names(residual_load_naive) <- names(ensembles_naive)

evaluate_ensembles <- function(ens_list, obs_df) {
  
  obs_vector <- as.numeric(obs_df$residual_load)
  n_days <- length(ens_list)
  
  results <- numeric(n_days)
  
  for(i in 1:n_days) {
    current_ens <- as.numeric(ens_list[[i]])
    current_obs <- obs_vector[i]
    results[i] <- crps_sample(y = current_obs, dat = current_ens)
  }
  
  return(data.frame(
    date = obs_df$date, 
    crps = results, 
    mean_crps = mean(results, na.rm = TRUE)
  ))
}

scoring_var_garch_res <- evaluate_ensembles(residual_load_var_garch, energy_test)
scoring_anen_res <- evaluate_ensembles(residual_load_anen, energy_test)
scoring_naive_res <- evaluate_ensembles(residual_load_naive, energy_test)

print(scoring_var_garch_res$mean_crps[1]); print(scoring_anen_res$mean_crps[1]); print(scoring_naive_res$mean_crps[1])

