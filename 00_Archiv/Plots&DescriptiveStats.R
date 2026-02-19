library(readr)
library(tseries)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RQuantLib)
library(combinat)
library(lubridate)
library(moments)
library(knitr)
library(ExtremeRisks)
library(corrplot)
library(fBasics)
library(viridis)
library(zoo)
library(forecast)

# Check if you adjusted the project directory correctly in main.R
load("00_Archiv/energy.RData")

## Plots
plot_time_series <- function(data, date_col, value_col, title, y_label = "", line_color = "#009682") {
  
  ggplot(data, aes(x = .data[[date_col]], y = .data[[value_col]])) +
    geom_line(color = line_color, linewidth = 0.8) +
    labs(
      title = title,
      x = "Date",
      y = y_label
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_time_series(energy, "date", "demand", "Energy Demand", y_label = "MWh")
plot_time_series(energy, "date", "solar", "Solar Power Production", y_label = "MWh")
plot_time_series(energy, "date", "wind", "Wind Power Production", y_label = "MWh")
plot_time_series(energy, "date", "demand_log", "Energy Demand (Log)", y_label = "Log(MWh)")
plot_time_series(energy, "date", "CF_solar", "Capacity Factor Solar", y_label = "")
plot_time_series(energy, "date", "CF_wind", "Capacity Factor Wind", y_label = "")
plot_time_series(energy, "date", "CF_solar_log", "Capacity Factor Solar (Log)", y_label = "")
plot_time_series(energy, "date", "CF_wind_log", "Capacity Factor Wind (Log)", y_label = "")

rownames_df <- c("demand", "solar", "wind", "demand_log", "CF_solar_log", "CF_wind_log")

descriptive_stats <- data.frame(
  min = rep(NA_real_, length(rownames_df)),
  max = rep(NA_real_, length(rownames_df)),
  mean = rep(NA_real_, length(rownames_df)),
  sd = rep(NA_real_, length(rownames_df)),
  skewness = rep(NA_real_, length(rownames_df)),
  kurtosis = rep(NA_real_, length(rownames_df)),
  JB_statistic = rep(NA_real_, length(rownames_df)),
  JB_p_val = rep(NA_real_, length(rownames_df)),
  row.names = rownames_df
)

for (var in rownames(descriptive_stats)) {
  
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

descriptive_stats_t <- as.data.frame(t(descriptive_stats))

## Volatility Clustering
plot_volatility_clustering <- function(data, date_col, value_col, title_label, n_intervals = 5, window = 30, yscale = "MWh") {
  
  plot_data <- data %>%
    dplyr::select(
      date  = all_of(date_col),
      value = all_of(value_col)
    ) %>%
    dplyr::arrange(date)
  
  plot_data$variance <- zoo::rollapply(
    plot_data$value,
    width = window,
    FUN = stats::var,
    fill = NA,
    align = "right"
  )
  
  plot_data <- plot_data %>% dplyr::filter(!is.na(variance))
  
  breaks <- quantile(
    plot_data$variance,
    probs = seq(0, 1, length.out = n_intervals + 1),
    na.rm = TRUE
  )
  
  print(breaks)
  
  plot_data$interval <- cut(
    plot_data$variance,
    breaks = breaks,
    labels = paste("Interval", 1:n_intervals),
    include.lowest = TRUE
  )
  
  plot_data <- plot_data %>%
    dplyr::mutate(
      date_end  = lead(date),
      value_end = lead(value)
    ) %>%
    dplyr::filter(!is.na(date_end))
  
  ggplot(plot_data) +
    geom_segment(
      aes(
        x = date, xend = date_end,
        y = value, yend = value_end,
        color = interval
      ),
      linewidth = 0.8,
      alpha = 0.9
    ) +
    scale_color_viridis_d(name = "Variance\nIntervals") +
    labs(
      title = paste("Volatility Clustering:", title_label),
      x = "Date",
      y = yscale
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  return(breaks)
}

p1 <- plot_volatility_clustering(energy, "date", "demand", "Energy Demand")
p2 <- plot_volatility_clustering(energy, "date", "demand_log", "Energy Demand Log", yscale = "Log(MWh)")
p3 <- plot_volatility_clustering(energy, "date", "solar", "Solar Power Production")
p4 <- plot_volatility_clustering(energy, "date", "CF_solar_log", "Capacity Factor Solar (Log)", yscale = "")
p5 <- plot_volatility_clustering(energy, "date", "wind", "Wind Power Production")
p6 <- plot_volatility_clustering(energy, "date", "CF_wind_log", "Capacity Factor Wind (Log)", yscale = "")

quantiles_table <- data.frame(
  demand = as.data.frame(p1),
  demand_log = as.data.frame(p2),
  solar = as.data.frame(p3),
  cf_solar = as.data.frame(plot_volatility_clustering(energy, "date", "CF_solar", "", yscale = "")),
  cf_solar_log = as.data.frame(p4),
  wind = as.data.frame(p5),
  cf_wind = as.data.frame(plot_volatility_clustering(energy, "date", "CF_wind", "", yscale = "")),
  cf_wind_log = as.data.frame(p6)
)

colnames(quantiles_table) <- c("demand", "demand_log", "solar", "cf_solar","cf_solar_log", "wind", "cf_wind","cf_wind_log")
quantiles_table["relative_length",] <- ((quantiles_table[nrow(quantiles_table),] / quantiles_table[1,])*100)
quantiles_table

## ACF / PACF
lag_plot_uni <- function(input_col, title){

  p1 <- forecast::ggAcf(input_col, lag.max = 15) +
    labs(title = title) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank()
    )
  
  p2 <- forecast::ggPacf(input_col, lag.max = 15) +
    labs(title = title) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank()
    )
  
  print(p1); print(p2)
}

lag_plot_uni(energy$demand, title = "Energy Demand")
lag_plot_uni(energy$solar, title = "Solar Power Production")
lag_plot_uni(energy$wind, title = "Wind Power Production")
lag_plot_uni(energy$demand_log, title = "Energy Demand (Log)")
lag_plot_uni(energy$CF_solar_log, title = "Capacity Factor Solar (Log)")
lag_plot_uni(energy$CF_wind_log, title = "Capacity Factor Wind (Log)")

lag_plot_bi <- function(input_col1, input_col2, title) {
  ggCcf(input_col1, input_col2, lag.max = 15, type = "correlation") +
    labs(title = title) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
    )
}

lag_plot_bi(energy$demand_log, energy$CF_solar_log, "CCF: Demand (Log) vs. CF Solar (Log)")
lag_plot_bi(energy$demand_log, energy$CF_wind_log, "CCF: Demand (Log) vs. CF Wind (Log)")
lag_plot_bi(energy$CF_solar_log, energy$CF_wind_log, "CCF: CF Solar (Log) vs. CF Wind (Log)")


