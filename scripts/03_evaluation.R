# ==============================================================================
# File: 03_evaluation.R
# Purpose: Evaluate forecasts using Scoring Rules (CRPS, ES), Statistical Tests
#          (Diebold-Mariano), Brier Score for Dark Doldrums, and Visualizations
#          (Heatmaps, PIT Histograms, Improvement Plots).
# ==============================================================================

# 1. Setup & Libraries ---------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(scoringRules)
library(forecast)
library(MultivCalibration)
library(RColorBrewer)

# Load helper functions
source("functions/eval_utils.R")

# 2. Load Data & Models --------------------------------------------------------
cat("Loading data and model results...\n")
load("output/data_processed_all.RData")        # energy, energy_test, etc.
load("output/models_and_ensembles.RData")      # ensembles_var_garch, etc.

# 3. CRPS & Energy Score Evaluation --------------------------------------------
cat("Calculating Scores (CRPS, Energy Score)...\n")

scoring_var_garch <- scoring(ensembles_var_garch, energy_test)
scoring_anen      <- scoring(ensembles_anen, energy_test)
scoring_naive     <- scoring(ensembles_naive, energy_test) # formerly rollmean

# Convert to tibble for easier handling later
scoring_var_garch <- as_tibble(scoring_var_garch)
scoring_anen      <- as_tibble(scoring_anen)
scoring_naive     <- as_tibble(scoring_naive)

# 4. Diebold-Mariano Tests -----------------------------------------------------
cat("Performing Diebold-Mariano Tests...\n")

# Trivariate ES
cat("\n--- DM Test: Trivariate ES ---\n")
print(forecast::dm.test(scoring_var_garch$es_trivar, scoring_anen$es_trivar, alternative = "less", h = 1, power = 1))
print(forecast::dm.test(scoring_var_garch$es_trivar, scoring_naive$es_trivar, alternative = "less", h = 1, power = 1))
print(forecast::dm.test(scoring_anen$es_trivar,      scoring_naive$es_trivar, alternative = "less", h = 1, power = 1))

# Bivariate ES
cat("\n--- DM Test: Bivariate ES ---\n")
print(forecast::dm.test(scoring_var_garch$es_bivar, scoring_anen$es_bivar, alternative = "less", h = 1, power = 1))
print(forecast::dm.test(scoring_var_garch$es_bivar, scoring_naive$es_bivar, alternative = "less", h = 1, power = 1))
print(forecast::dm.test(scoring_anen$es_bivar,      scoring_naive$es_bivar, alternative = "less", h = 1, power = 1))

# 5. Dark Doldrum Analysis (Table) ---------------------------------------------
cat("Analyzing Dark Doldrum Performance...\n")

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
    mean_ES_Naive = mean(es_bivar_naive, na.rm = T),
    ESS_Var_Garch_vs_AnEn_Pct = (mean_ES_AnEn - mean_ES_Var_Garch) / mean_ES_AnEn * 100,
    ESS_Var_Garch_vs_Naive_Pct = (mean_ES_Naive - mean_ES_Var_Garch) / mean_ES_Naive * 100
  )

print(scoring_summary)

# 6. Forecaster's Dilemma ------------------------------------------------------
cat("Evaluation based on the Forecaster's Dilemma... \n")

# Define Color-Setting
my_colors <- c(
  "VAR-GARCH"   = "#1B9783",  
  "AnEn"        = "#377EB8",  
  "Climatology" = "#999999"  
)

# 6.1 Stratify on Forecasted Dark Doldrum Probability --------------------------
cat("Stratify on Forecasted Dark Doldrum Probability... \n")
probs_var <- h48_mean_scatter(ensembles_var_garch, energy_test, threshold) %>% dplyr::select(date, prob_var = prob)
probs_anen <- h48_mean_scatter(ensembles_anen, energy_test, threshold) %>%dplyr::select(date, prob_anen = prob)
probs_naive <- h48_mean_scatter(ensembles_naive, energy_test, threshold) %>% dplyr::select(date, prob_naive = prob)

mean_combined <- data.frame(
  date = probs_var$date,
  prob_combined = (probs_var$prob_var + probs_anen$prob_anen + probs_naive$prob_naive) / 3
)

strat_df <- scoring_var_garch %>%
  dplyr::select(date, es_var = es_bivar) %>%
  left_join(scoring_anen %>% dplyr::select(date, es_anen = es_bivar), by = "date") %>%
  left_join(scoring_naive %>% dplyr::select(date, es_naive = es_bivar), by = "date") %>%
  left_join(probs_var, by = "date") %>%
  left_join(probs_anen, by = "date") %>%
  left_join(probs_naive, by = "date") %>%
  left_join(mean_combined, by = "date") %>%
  inner_join(energy_test %>% dplyr::select(date, dark_doldrum), by = "date")

strat_df <- strat_df[-1,] # Remove first row, since there is NA due to 48h mean

calc_stratified_score <- function(data, prob_col, es_col, model_label) {
  data %>%
    dplyr::select(prob = all_of(prob_col), es = all_of(es_col)) %>%
    mutate(
      # Define risk bins for dark doldrums
      risk_bin = cut(prob, 
                     breaks = c(-Inf, 0.1, 0.5, Inf),  
                     labels = c("Low Risk (<10%)", "Medium Risk (10-50%)", "High Risk (>50%)"),
                     include.lowest = T) 
    ) %>%
    group_by(risk_bin) %>%
    summarise(
      Model = model_label,
      n_days = n(),
      mean_ES = mean(es, na.rm = T),
      .groups = "drop"
    )
}

# Apply to all models
(res_var <- calc_stratified_score(strat_df, "prob_combined", "es_var", "VAR-GARCH"))
(res_anen <- calc_stratified_score(strat_df, "prob_combined", "es_anen", "AnEn"))
(res_naive <- calc_stratified_score(strat_df, "prob_combined", "es_naive", "Climatology"))

# Combine results
stratified_results <- bind_rows(res_var, res_anen, res_naive)
print(stratified_results)

# Visualize
p_strat <- ggplot(stratified_results, aes(x = risk_bin, y = mean_ES, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Model Performance Conditioned on Forecasted Risk",
    x = "Forecasted Probability of Dark Doldrum",
    y = "Mean Bivariate Energy Score",
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_strat)
ggsave(filename = file.path('output/plots', "Forecasters_Dilemma_DD.pdf"), plot = p_strat, width = 10, height = 5, bg = "white")

# 6.2 Stratify on Forecasted Demand Quartiles ----------------------------------
cat("Stratify on Forecasted Demand Quartiles... \n")

get_ens_mean_demand <- function(ens_list) {
  sapply(ens_list, function(mat) mean(mat[, 1]))
}

get_ens_median_demand <- function(ens_list){
  sapply(ens_list, function(mat) median(mat[,1]))
}

# Median from each model
fc_demand_var_median <- as.data.frame(get_ens_median_demand(ensembles_var_garch))
fc_demand_anen_median <- as.data.frame(get_ens_median_demand(ensembles_anen))
fc_demand_naive_median <- as.data.frame(get_ens_median_demand(ensembles_naive))

fc_demand_median_mean <- (fc_demand_var_median + fc_demand_anen_median + fc_demand_naive_median) / 3
colnames(fc_demand_median_mean) <- "fc_demand_median_mean"


demand_eval_df <- tibble(
  date = energy_test$date,
  # Forecasted Demand
  fc_median_mean = fc_demand_median_mean$fc_demand_median_mean,
  # ES (Trivariate)
  es_var   = scoring_var_garch$es_trivar,
  es_anen  = scoring_anen$es_trivar,
  es_naive = scoring_naive$es_trivar
)

calc_demand_strat <- function(data, fc_col, es_col, model_label) {
  
  df_model <- data %>%
    dplyr::select(forecast = all_of(fc_col), es = all_of(es_col)) %>%
    filter(!is.na(forecast), !is.na(es)) 
  
  quantiles <- quantile(df_model$forecast, probs = c(0, 0.25, 0.5, 0.75, 1)) # Define Quartiles
  
  df_model %>%
    mutate(
      demand_bin = cut(forecast, 
                       breaks = quantiles,
                       labels = c("Q1 (Low Demand)", "Q2 (Medium-Low Demand)", "Q3 (Medium-High Demand)", "Q4 (High Demand)"),
                       include.lowest = T)
    ) %>%
    group_by(demand_bin) %>%
    summarise(
      Model = model_label,
      n_days = n(),
      mean_trivar_ES = mean(es, na.rm = T),
      .groups = "drop"
    )
}

# Apply to all models
(res_demand_var <- calc_demand_strat(demand_eval_df, "fc_median_mean", "es_var", "VAR-GARCH"))
(res_demand_anen <- calc_demand_strat(demand_eval_df, "fc_median_mean", "es_anen", "AnEn"))
(res_demand_naive <- calc_demand_strat(demand_eval_df, "fc_median_mean", "es_naive", "Climatology"))

demand_results <- bind_rows(res_demand_var, res_demand_anen, res_demand_naive)
print(demand_results)

# Visualize
p_demand <- ggplot(demand_results, aes(x = demand_bin, y = mean_trivar_ES, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Trivariate Performance Stratified by Forecasted Demand",
    x = "Forecasted Demand Quartile",
    y = "Mean Trivariate Energy Score",
    caption = "Quartiles are calculated individually per model forecast to avoid Forecaster's Dilemma."
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_demand)
ggsave(filename = file.path('output/plots', "Forecasters_Dilemma_Demand.pdf"), plot = p_demand, width = 10, height = 5, bg = "white")

# 6.3 Stratify on Forecasted Dark Doldrum Probability AND Forecasted Demand Quartiles
cat("Stratify on Forecasted Dark Doldrum Probability AND Forecasted Demand Quartiles... \n")

prepare_model_df <- function(dates, fc_demand, fc_prob, es_score, model_name) {
  tibble(
    date = dates[-1],
    forecast_demand = fc_demand[-1],
    forecast_prob = fc_prob,
    es = es_score[-1],
    Model = model_name
  )
}

# Apply to all models in Median Mean Demand and Mean Dark Doldrum Probability
df_var <- prepare_model_df(energy_test$date, fc_demand_median_mean$fc_demand_median_mean, strat_df$prob_combined, scoring_var_garch$es_trivar, "VAR-GARCH")
df_anen <- prepare_model_df(energy_test$date, fc_demand_median_mean$fc_demand_median_mean, strat_df$prob_combined, scoring_anen$es_trivar, "AnEn")
df_naive <- prepare_model_df(energy_test$date, fc_demand_median_mean$fc_demand_median_mean, strat_df$prob_combined, scoring_naive$es_trivar, "Climatology")

calc_double_strat <- function(df_model) {
  
  # Demand Quartiles
  quantiles_dem <- quantile(df_model$forecast_demand, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  
  df_model %>%
    filter(!is.na(forecast_demand), !is.na(forecast_prob), !is.na(es)) %>%
    mutate(
      demand_bin = cut(forecast_demand, 
                       breaks = quantiles_dem,
                       labels = c("Q1 (Low Demand)", "Q2 (Medium-Low Demand)", "Q3 (Medium-High Demand)", "Q4 (High Demand)"),
                       include.lowest = T),
      
      risk_bin = cut(forecast_prob,
                     breaks = c(-Inf, 0.1, 0.5, Inf),
                     labels = c("Low Risk (<10%)", "Medium Risk (10-50%)", "High Risk (>50%)"),
                     include.lowest = T)
    ) %>%
    # Gruoup by both factors
    group_by(demand_bin, risk_bin) %>%
    summarise(
      Model = unique(Model),
      n_days = n(),
      mean_ES = mean(es, na.rm = T),
      .groups = "drop"
    )
}

# Combine results
double_strat_results <- bind_rows(calc_double_strat(df_var), calc_double_strat(df_anen), calc_double_strat(df_naive))
print(double_strat_results %>% dplyr::select(Model, demand_bin, risk_bin, n_days) %>% arrange(risk_bin, demand_bin))

# Visualize
plot_data <- double_strat_results

p_double <- ggplot(plot_data, aes(x = demand_bin, y = mean_ES, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~risk_bin, ncol = 1, scales = "fixed") + 
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Trivariate ES Conditioned on Forecasted Demand and Forecasted Risk",
    x = "Forecasted Demand Quartile",
    y = "Mean Trivariate Energy Score",
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

print(p_double)
ggsave(filename = file.path('output/plots', "Forecasters_Dilemma_DD_Demand.pdf"), plot = p_double, width = 10, height = 5, bg = "white")

# 7. Days of Medium/ High Dark Doldrum Risk & High Energy Demand ---------------
extract_days <- function(input_df, risk_boolean){ # If risk_boolean = TRUE -> Medium and High Dark Doldrum Risk. Else: Only Medium Dark Doldrum Risk
  
  # Pre-Filter according to Forecasters Dilemma
  demand_quartiles <- quantile(input_df$forecast_demand, 0.75)
  
  lower_risk_bound <- 0.1
  
  if(!risk_boolean){
    upper_risk_bound <- 0.5
    input_df_subset <- subset(input_df, input_df$forecast_demand >= demand_quartiles & input_df$forecast_prob >= lower_risk_bound & input_df$forecast_prob <= upper_risk_bound)
  }
  else{
    input_df_subset <- subset(input_df, input_df$forecast_demand >= demand_quartiles & input_df$forecast_prob >= lower_risk_bound)
  }
  
  # Filter in Energy Test Data Set for high grid stress days (i.e. high demand & dark doldrum event)
  energy_test_subset <- filter(energy_test, energy_test$date %in% input_df_subset$date)
  
  # Calculate the portion of the residual load from the energy demand
  energy_test_subset$portion <- energy_test_subset$residual_load / energy_test_subset$demand
  
  print("Note: Portions in %")
  print(paste("Minimum Portion:", (min(energy_test_subset$portion))*100))
  print(paste("Maximum Portion:", (max(energy_test_subset$portion))*100))
  print(paste("Mean Portion:", (mean(energy_test_subset$portion))*100))
  
  return(energy_test_subset)
}

# Apply to respective models
residual_df_var <- extract_days(df_var, risk_boolean = F)
residual_df_anen <- extract_days(df_anen, risk_boolean = F)
residual_df_naive <- extract_days(df_naive, risk_boolean = T) # TRUE because of the little days in the High Demand & High Risk bin

# 8. Brier Score Evaluation ----------------------------------------------------
cat("Calculating Brier Scores...\n")

# Use threshold from data_processed_all.RData
brier_df_an_en     <- h48_mean_scatter(ensembles_anen, energy_test, threshold)
brier_df_naive     <- h48_mean_scatter(ensembles_naive, energy_test, threshold)
brier_df_var_garch <- h48_mean_scatter(ensembles_var_garch, energy_test, threshold)

bs_anen      <- brier_score(brier_df_an_en)
bs_naive     <- brier_score(brier_df_naive)
bs_var_garch <- brier_score(brier_df_var_garch)

brier_table <- data.frame(
  AnEn = bs_anen,
  Climatology_Benchmark = bs_naive,
  Var_Garch = bs_var_garch,
  BSS_Var_Garch_vs_AnEn_Pct = ((-1) * ((bs_var_garch - bs_anen) / bs_anen))*100, 
  BSS_Var_Garch_vs_Climatology_Pct = ((-1) * ((bs_var_garch - bs_naive) / bs_naive))*100
)
rownames(brier_table) <- "Brier_Score"
print(brier_table)

# Scatterplots
pdf_path <- file.path('output/plots', "daily_48h_analysis.pdf")
pdf(pdf_path, width = 13, height = 5.5)

for (t in 2:nrow(energy_test)) {
  
  current_date <- energy_test$date[t]
  obs_curr     <- energy_test[t, ]
  obs_prev     <- energy_test[t-1, ]
  
  limit_48h <- threshold * obs_curr$P_max_total
  if (length(limit_48h) == 0 || is.na(limit_48h)) next
  
  obs_plot_df <- data.frame(
    solar_48h = (obs_prev$solar + obs_curr$solar) / 2,
    wind_48h  = (obs_prev$wind + obs_curr$wind) / 2
  )
  
  calc_48h_data <- function(ens_mat, obs_p, label, limit) {
    df <- as.data.frame(ens_mat)
    colnames(df) <- c("demand", "solar", "wind")
    df_res <- data.frame(
      solar_48h = (obs_p$solar + df$solar) / 2,
      wind_48h  = (obs_p$wind + df$wind) / 2,
      Model = label
    )
    df_res$is_doldrum_point <- (df_res$solar_48h + df_res$wind_48h) < limit
    return(df_res)
  }
  
  df_var   <- calc_48h_data(ensembles_var_garch[[t]], obs_prev, "VAR-GARCH", limit_48h)
  df_anen  <- calc_48h_data(ensembles_anen[[t]],      obs_prev, "AnEn", limit_48h)
  df_naive <- calc_48h_data(ensembles_naive[[t]],     obs_prev, "Climatology", limit_48h)
  
  prob_df <- data.frame(
    Model = c("VAR-GARCH", "AnEn", "Climatology"),
    f_t   = c(mean(df_var$is_doldrum_point), 
              mean(df_anen$is_doldrum_point), 
              mean(df_naive$is_doldrum_point))
  )
  prob_df$label_text <- paste0("f_t = ", round(prob_df$f_t * 100, 1), "%")
  
  df_combined <- bind_rows(df_var, df_anen, df_naive)
  
  p <- ggplot(df_combined, aes(x = solar_48h, y = wind_48h)) +
    geom_point(aes(color = is_doldrum_point), alpha = 0.8, size = 0.8) +
    scale_color_manual(values = c("FALSE" = "darkgrey", "TRUE" = "red"), guide = "none") +
    
    geom_abline(intercept = limit_48h, slope = -1, linetype = "dashed", color = "black", linewidth = 0.7) +
    
    geom_point(data = obs_plot_df, color = "#009682", size = 4, shape = 16) +
    
    geom_text(data = prob_df, aes(label = label_text), 
              x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, 
              fontface = "bold", size = 4.5, color = "black") +
    
    facet_wrap(~ Model) +
    labs(
      title = paste("48h-Mean Analysis:", current_date),
      subtitle = paste("Actual Obs Dark Doldrum:", 
                       ifelse((obs_plot_df$solar_48h + obs_plot_df$wind_48h) < limit_48h, "True", "False")),
      x = "Solar (MWh)",
      y = "Wind (MWh)",
      caption = "f_t: Fraction of ensemble members below threshold (red points)"
    ) +
    theme_minimal() +
    scale_x_continuous(
      limits = c(0, max(energy_test$solar, na.rm = TRUE)),
      breaks = scales::breaks_pretty(n = 6),
      
      labels = function(x) {
        lbls <- scales::label_scientific(digits = 1)(x)
        lbls[seq_along(lbls) %% 2 != 0] <- "" 
        return(lbls)
      }
    ) +
    scale_y_continuous(
      limits = c(0, max(energy_test$wind, na.rm = TRUE)),
      labels = scales::label_scientific(digits = 1),
      breaks = scales::breaks_pretty(n = 5)
    )
  
  print(p)
}
dev.off()

# 9. Heatmap Plot --------------------------------------------------------------
cat("Generating Heatmap...\n")

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

p_heatmap <- ggplot(hm_all, aes(x = month, y = doy, fill = prob)) +
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

while (!is.null(dev.list())) dev.off()

print(p_heatmap)
ggsave(filename = file.path('output/plots', "heatmap_dark_doldrum.pdf"), plot = p_heatmap, width = 10, height = 5, bg = "white")

# 10. Improvement Bar Plot (by Demand Quartile) --------------------------------
cat("Generating Improvement Bar Plot...\n")

# Prepare Demand Quartiles
quantiles_by_wday <- energy_train %>%
  group_by(wday) %>%
  summarise(
    q1 = quantile(demand, 0.25),
    q2 = quantile(demand, 0.50),
    q3 = quantile(demand, 0.75)
  )

energy_test_quartiles <- energy_test %>%
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

# Aggregate Comparisons
comparison_var_garch <- scoring_var_garch %>%
  inner_join(energy_test_quartiles %>% dplyr::select(date, dark_doldrum, demand_quartile), by = "date") %>%
  group_by(dark_doldrum, demand_quartile) %>%
  summarise(
    n = n(),
    es_trivar_mean = mean(es_trivar, na.rm = T),
    .groups = "drop"
  )

comparison_anen <- scoring_anen %>%
  inner_join(energy_test_quartiles %>% dplyr::select(date, dark_doldrum, demand_quartile), by = "date") %>%
  group_by(dark_doldrum, demand_quartile) %>%
  summarise(es_trivar_mean = mean(es_trivar, na.rm = T), .groups = "drop")

comparison_naive <- scoring_naive %>%
  inner_join(energy_test_quartiles %>% dplyr::select(date, dark_doldrum, demand_quartile), by = "date") %>%
  group_by(dark_doldrum, demand_quartile) %>%
  summarise(es_trivar_mean = mean(es_trivar, na.rm = T), .groups = "drop")

# Merge and Calculate Improvement
improvement_df <- comparison_var_garch %>%
  left_join(comparison_anen, by = c("dark_doldrum", "demand_quartile"), suffix = c("", "_anen")) %>%
  left_join(comparison_naive, by = c("dark_doldrum", "demand_quartile"), suffix = c("", "_naive")) %>%
  mutate(
    es_trivar_mean_var_garch_vs_anen_pct = ((es_trivar_mean_anen - es_trivar_mean) / es_trivar_mean_anen) * 100,
    es_trivar_mean_var_garch_vs_naive_pct = ((es_trivar_mean_naive - es_trivar_mean) / es_trivar_mean_naive) * 100
  )

print(improvement_df)

# Plot for Dark Doldrums only
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
      comparison == "es_trivar_mean_var_garch_vs_naive_pct" ~ "VAR-GARCH vs Climatology"
    )
  )

p_imp <- ggplot(plot_data,
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

print(p_imp)
ggsave(filename = file.path('output/plots', "barplot_improvement.pdf"), plot = p_imp, width = 8, height = 5, bg = "white")

# 11. PIT Histograms -----------------------------------------------------------
cat("Generating PIT Histograms...\n")

# VAR-GARCH
plot_pit_univariate(ensembles_var_garch, energy_test, "- VAR-GARCH")
plot_pit_multivariate(ensembles_var_garch, energy_test, "- VAR-GARCH")

# Analog Ensemble
plot_pit_univariate(ensembles_anen, energy_test, "- AnEn")
plot_pit_multivariate(ensembles_anen, energy_test, "- AnEn")

# Climatology Benchmark
plot_pit_univariate(ensembles_naive, energy_test, "- Climatology Benchmark")
plot_pit_multivariate(ensembles_naive, energy_test, "- Climatology Benchmark")

# 12. Residual Load Evaluation -------------------------------------------------
cat("Evaluating Residual Load CRPS...\n")

# Calculate Residual Loads for all ensembles
calc_res_load <- function(ens_list) {
  lapply(ens_list, function(mat) mat[,1] - mat[,2] - mat[,3])
}

res_load_var_garch <- calc_res_load(ensembles_var_garch)
res_load_anen      <- calc_res_load(ensembles_anen)
res_load_naive     <- calc_res_load(ensembles_naive)

# Evaluate using helper function
score_res_var_garch <- evaluate_ensembles(res_load_var_garch, energy_test)
score_res_anen      <- evaluate_ensembles(res_load_anen, energy_test)
score_res_naive     <- evaluate_ensembles(res_load_naive, energy_test)

cat("\n--- Mean CRPS for Residual Load ---\n")
cat("VAR-GARCH:", score_res_var_garch$mean_crps[1], "\n")
cat("AnEn:     ", score_res_anen$mean_crps[1], "\n")
cat("Naive:    ", score_res_naive$mean_crps[1], "\n")

cat("\nEvaluation Script Complete.\n")
