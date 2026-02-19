# ==============================================================================
# File: stationarity_utils.R
# Purpose: Functions for evaluating the stationarity of the multivariate
#          time series
# Dependencies: tseries, VAR, urca
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Univariate Analysis of Stationarity
# ------------------------------------------------------------------------------

stationarity_check <- function(x){
  cat("\nADF (alt='stationary'):\n");  print(adf.test(x, alternative = "stationary"))
  cat("KPSS (null='Trend'):\n");      print(kpss.test(x, null = "Trend"))
}

# ------------------------------------------------------------------------------
# 2. Johansen test for cointegration
# ------------------------------------------------------------------------------

johansen_test <- function(x){
  lag_selecion <- VARselect(x[,-1], lag.max = 31, type = "const")
  K_train <- lag_selecion$selection["AIC(n)"]
  johansen <- urca::ca.jo(x[,-1], type = "trace", spec = "longrun", K = K_train, ecdet = "const")
  return(summary(johansen))
}

