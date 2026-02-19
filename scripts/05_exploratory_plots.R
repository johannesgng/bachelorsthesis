# ==============================================================================
# File: 05_exploratory_plots.R
# Purpose: Exploratory Data Analysis (EDA) including:
#          - Time Series Plots (Levels vs. Log)
#          - Descriptive Statistics (Moments, Jarque-Bera Test)
#          - Volatility Clustering Visualization
#          - Autocorrelation (ACF/PACF) and Cross-Correlation (CCF) Analysis
# Output:  Saves plots to 'output/plots/exploratory/' and tables as CSV.
# ==============================================================================

# 1. Setup & Libraries ---------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(tseries)
library(moments)
library(corrplot)
library(viridis)
library(zoo)
library(forecast)

# Define Output Directory
plot_dir <- "output/plots"

# 2. Load Data -----------------------------------------------------------------
cat("Loading processed data...\n")
load("output/data_processed_all.RData") 

# 3. Time Series Plots ---------------------------------------------------------
cat("Generating Time Series Plots...\n")

plot_time_series <- function(data, date_col, value_col, title, y_label = "", filename) {
  
  p <- ggplot(data, aes(x = .data[[date_col]], y = .data[[value_col]])) +
    geom_line(color = "#009682", linewidth = 0.6) + # KIT Green
    labs(
      title = title,
      x = "Date",
      y = y_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  ggsave(filename = file.path(plot_dir, filename), plot = p, width = 10, height = 5, bg = "white")
}

# Generate and Save Plots
plot_time_series(energy, "date", "demand", "Energy Demand", "MWh", "ts_demand.pdf")
plot_time_series(energy, "date", "solar", "Solar Power Production", "MWh", "ts_solar.pdf")
plot_time_series(energy, "date", "wind", "Wind Power Production", "MWh", "ts_wind.pdf")
plot_time_series(energy, "date", "demand_log", "Energy Demand (Log)", "Log(MWh)", "ts_demand_log.pdf")
plot_time_series(energy, "date", "CF_solar", "Capacity Factor Solar", "Ratio", "ts_cf_solar.pdf")
plot_time_series(energy, "date", "CF_wind", "Capacity Factor Wind", "Ratio", "ts_cf_wind.pdf")
plot_time_series(energy, "date", "CF_solar_log", "Capacity Factor Solar (Log)", "Log(Ratio)", "ts_cf_solar_log.pdf")
plot_time_series(energy, "date", "CF_wind_log", "Capacity Factor Wind (Log)", "Log(Ratio)", "ts_cf_wind_log.pdf")

# 4. Descriptive Statistics Table ----------------------------------------------
cat("Calculating Descriptive Statistics...\n")

vars_to_analyze <- c("demand", "solar", "wind", "demand_log", "CF_solar_log", "CF_wind_log")

descriptive_stats <- data.frame(
  min = numeric(length(vars_to_analyze)),
  max = numeric(length(vars_to_analyze)),
  mean = numeric(length(vars_to_analyze)),
  sd = numeric(length(vars_to_analyze)),
  skewness = numeric(length(vars_to_analyze)),
  kurtosis = numeric(length(vars_to_analyze)),
  JB_statistic = numeric(length(vars_to_analyze)),
  JB_p_val = numeric(length(vars_to_analyze)),
  row.names = vars_to_analyze
)

for (var in vars_to_analyze) {
  x <- energy[[var]]
  
  descriptive_stats[var, "min"] <- min(x)
  descriptive_stats[var, "max"] <- max(x)
  descriptive_stats[var, "mean"] <- mean(x)
  descriptive_stats[var, "sd"] <- sd(x)
  descriptive_stats[var, "skewness"] <- skewness(x)
  descriptive_stats[var, "kurtosis"] <- kurtosis(x)
  
  jb <- jarque.bera.test(x)
  descriptive_stats[var, "JB_statistic"] <- jb$statistic
  descriptive_stats[var, "JB_p_val"]     <- jb$p.value
}

# Save Table
write.csv(descriptive_stats, file.path("output", "table_descriptive_statistics.csv"))
print(descriptive_stats)

# 5. Volatility Clustering Visualization ---------------------------------------
cat("Visualizing Volatility Clustering...\n")

plot_volatility_clustering <- function(data, date_col, value_col, title_label, filename, n_intervals = 5, window = 30, yscale = "MWh") {
  
  plot_data <- data %>%
    dplyr::select(date = all_of(date_col), value = all_of(value_col)) %>%
    dplyr::arrange(date)
  
  # Calculate rolling variance
  plot_data$variance <- zoo::rollapply(
    plot_data$value, width = window, FUN = stats::var, fill = NA, align = "right"
  )
  
  plot_data <- plot_data %>% dplyr::filter(!is.na(variance))
  
  breaks <- quantile(plot_data$variance, probs = seq(0, 1, length.out = n_intervals + 1), na.rm = TRUE)
  
  plot_data$interval <- cut(
    plot_data$variance, breaks = breaks, labels = paste("Interval", 1:n_intervals), include.lowest = TRUE
  )
  
  # Prepare segments for plotting
  plot_data <- plot_data %>%
    dplyr::mutate(date_end = lead(date), value_end = lead(value)) %>%
    dplyr::filter(!is.na(date_end))
  
  p <- ggplot(plot_data) +
    geom_segment(aes(x = date, xend = date_end, y = value, yend = value_end, color = interval),
                 linewidth = 0.8, alpha = 0.9) +
    scale_color_viridis_d(name = "Variance Intervals") +
    labs(
      title = paste("Volatility Clustering:", title_label),
      x = "Date", y = yscale
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  print(p)
  ggsave(filename = file.path(plot_dir, filename), plot = p, width = 10, height = 6, bg = "white")
  
  return(breaks)
}

# Generate Volatility Plots and collect breaks
breaks_demand <- plot_volatility_clustering(energy, "date", "demand", "Energy Demand", "vol_clust_demand.pdf")
breaks_demand_log <- plot_volatility_clustering(energy, "date", "demand_log", "Energy Demand Log", "vol_clust_demand_log.pdf", yscale = "Log(MWh)")
breaks_solar <- plot_volatility_clustering(energy, "date", "solar", "Solar Power", "vol_clust_solar.pdf")
breaks_cf_solar_log <- plot_volatility_clustering(energy, "date", "CF_solar_log", "CF Solar (Log)", "vol_clust_cf_solar_log.pdf", yscale = "Log")
breaks_wind <- plot_volatility_clustering(energy, "date", "wind", "Wind Power", "vol_clust_wind.pdf")
breaks_cf_wind_log <- plot_volatility_clustering(energy, "date", "CF_wind_log", "CF Wind (Log)", "vol_clust_cf_wind_log.pdf", yscale = "Log")

# 6. Autocorrelation (ACF/PACF) ------------------------------------------------
cat("Generating ACF/PACF Plots...\n")

save_acf_pacf <- function(input_col, title, filename_base) {
  
  p1 <- forecast::ggAcf(input_col, lag.max = 30) +
    labs(title = paste(title, "- ACF")) +
    theme_minimal()
  
  p2 <- forecast::ggPacf(input_col, lag.max = 30) +
    labs(title = paste(title, "- PACF")) +
    theme_minimal()
  
  # Arrange plots side by side (using gridExtra logic or just saving separately)
  # Here saving separately for cleaner inclusion in LaTeX
  ggsave(file.path(plot_dir, paste0(filename_base, "_ACF.pdf")), p1, width = 6, height = 4, bg = "white")
  ggsave(file.path(plot_dir, paste0(filename_base, "_PACF.pdf")), p2, width = 6, height = 4, bg = "white")
}

save_acf_pacf(energy$demand, "Energy Demand", "acf_demand")
save_acf_pacf(energy$demand_log, "Energy Demand (Log)", "acf_demand_log")
save_acf_pacf(energy$CF_solar_log, "CF Solar (Log)", "acf_cf_solar_log")
save_acf_pacf(energy$CF_wind_log, "CF Wind (Log)", "acf_cf_wind_log")

# 7. Cross-Correlation Function (CCF) ------------------------------------------
cat("Generating Cross-Correlation Plots...\n")

save_ccf_plot <- function(col1, col2, title, filename) {
  p <- ggCcf(col1, col2, lag.max = 30) +
    labs(title = title) +
    theme_minimal()
  
  ggsave(file.path(plot_dir, filename), p, width = 8, height = 5, bg = "white")
}

save_ccf_plot(energy$demand_log, energy$CF_solar_log, "CCF: Demand (Log) vs. CF Solar (Log)", "ccf_demand_solar.pdf")
save_ccf_plot(energy$demand_log, energy$CF_wind_log, "CCF: Demand (Log) vs. CF Wind (Log)", "ccf_demand_wind.pdf")
save_ccf_plot(energy$CF_solar_log, energy$CF_wind_log, "CCF: CF Solar (Log) vs. CF Wind (Log)", "ccf_solar_wind.pdf")
