# ==============================================================================
# File: eval_utils.R
# Purpose: Functions for evaluating forecast performance including 
#          Scoring Rules (CRPS, ES), Brier Score for Dark Doldrums, 
#          and PIT Histograms (Calibration).
# Dependencies: scoringRules, dplyr, MultivCalibration
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. General Scoring (CRPS & Energy Score)
# ------------------------------------------------------------------------------

scoring <- function(ensembels, data){
  # Requires 'scoringRules' and 'dplyr'
  
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
    
    # Note: es_sample expects observation vector y and ensemble matrix dat (d x m)
    results_df$es_trivar[i] <- scoringRules::es_sample(y = obs_vec, dat = t(ens_mat))
    
    # Bivariate Energy Score
    ens_mat_bivar <- ensembels[[i]][,c("solar", "wind")]
    
    obs_vec_bivar <- c(obs_solar, obs_wind)
    
    results_df$es_bivar[i] <- scoringRules::es_sample(y = obs_vec_bivar, dat = t(ens_mat_bivar))
    
    # CRPS
    results_df$crps_demand[i] <- scoringRules::crps_sample(y = obs_demand, dat = as.numeric(ens_mat[,1, drop = F]))
    results_df$crps_solar[i] <- scoringRules::crps_sample(y = obs_solar, dat = as.numeric(ens_mat[,2, drop = F]))
    results_df$crps_wind[i] <- scoringRules::crps_sample(y = obs_wind, dat = as.numeric(ens_mat[,3, drop = F]))
  }
  
  results_df <- results_df %>% dplyr::mutate(
    es_trivar_mean = mean(es_trivar),
    es_bivar_mean = mean(es_bivar),
    crps_demand_mean = mean(crps_demand),
    crps_solar_mean = mean(crps_solar),
    crps_wind_mean = mean(crps_wind)
  )
  
  return(results_df)
}

# ------------------------------------------------------------------------------
# 2. Brier Score & Dark Doldrum Detection
# ------------------------------------------------------------------------------

h48_mean_scatter <- function(ensembles, obs_data, threshold_cf = 0.06) {
  # Note: Default threshold set to 0.06 to avoid dependency on global variable
  
  results_df <- data.frame(
    date = as.Date(names(ensembles)),
    prob = NA_real_,
    obs_event = NA_real_,
    M = NA_integer_
  )
  
  dates <- as.Date(names(ensembles))
  
  for (t in seq_along(dates)) {
    if (t == 1) next
    
    obs_prev <- obs_data %>% dplyr::filter(date == dates[t-1])
    obs_curr <- obs_data %>% dplyr::filter(date == dates[t])
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
  }
  
  results_df <- results_df %>% dplyr::filter(!is.na(prob))
  return(results_df)
}

brier_score <- function(df){ 
  # Calculates Brier Score
  result_bs <- (df$prob - df$obs_event)^2
  print(mean(result_bs))
  return(mean(result_bs))
}

# ------------------------------------------------------------------------------
# 3. PIT Histograms (Calibration)
# ------------------------------------------------------------------------------

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
  # Requires 'MultivCalibration'
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

# ------------------------------------------------------------------------------
# 4. Residual Load Evaluation
# ------------------------------------------------------------------------------

evaluate_ensembles <- function(ens_list, obs_df) {
  
  obs_vector <- as.numeric(obs_df$residual_load)
  n_days <- length(ens_list)
  
  results <- numeric(n_days)
  
  for(i in 1:n_days) {
    current_ens <- as.numeric(ens_list[[i]])
    current_obs <- obs_vector[i]
    results[i] <- scoringRules::crps_sample(y = current_obs, dat = current_ens)
  }
  
  return(data.frame(
    date = obs_df$date, 
    crps = results, 
    mean_crps = mean(results, na.rm = TRUE)
  ))
}