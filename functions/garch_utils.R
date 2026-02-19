# ==============================================================================
# File: garch_utils.R
# Purpose: Functions for stationarity testing and GARCH model selection 
#          (Grid Search).
# Dependencies: tseries, rugarch, FinTS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Stationarity Checks (ADF & KPSS)
# ------------------------------------------------------------------------------

stationarity_check <- function(x){
  # Requires 'tseries' library
  cat("\nADF (alt='stationary'):\n");  print(tseries::adf.test(x, alternative = "stationary"))
  cat("KPSS (null='Trend'):\n");      print(tseries::kpss.test(x, null = "Trend"))
}

# ------------------------------------------------------------------------------
# 2. GARCH Grid Search
# ------------------------------------------------------------------------------

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
        uspec_try <- try(rugarch::ugarchspec(
          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
          variance.model = list(model = m, garchOrder = o), 
          distribution.model = d
        ), silent = TRUE)
        
        if (inherits(uspec_try, "try-error")) next
        
        fit_try <- try(rugarch::ugarchfit(
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
        current_aic <- rugarch::infocriteria(fit_try)["Akaike",]
        
        # Residual Diagnostics
        resid_std <- residuals(fit_try, standardize = TRUE)
        
        # ARCH Test
        # Requires 'FinTS' library
        arch_res <- try(FinTS::ArchTest(resid_std, lags = 12), silent = TRUE) # Lags etwas reduzieren für Robustheit
        arch_p <- if (inherits(arch_res, "try-error")) 0 else arch_res$p.value
        
        # Ljung-Box
        lb_res <- try(Box.test(resid_std, lag = 12, type = "Ljung-Box"), silent = TRUE)
        lb_p <- if (inherits(lb_res, "try-error")) 0 else lb_res$p.value
        
        # Nyblom 
        nyblom_res <- try(rugarch::nyblom(fit_try), silent = TRUE)
        nyblom_stable <- FALSE
        if (!inherits(nyblom_res, "try-error")) {
          stat <- nyblom_res[["JointStat"]]
          crit <- nyblom_res[["JointCritical"]][["5%"]] # 5%
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