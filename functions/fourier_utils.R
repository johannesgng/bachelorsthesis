# ==============================================================================
# File: fourier_utils.R
# Purpose: Functions for creating Fourier series, fitting linear models, 
#          finding optimal seasonality orders, and reconstructing forecasts.
# ==============================================================================

# Note: Ensure 'lubridate' is loaded in the main script or load it here
# library(lubridate)

# ------------------------------------------------------------------------------
# 1. Fourier Basis Functions
# ------------------------------------------------------------------------------

# Annual terms (Period = 365)
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

# Monthly terms (Period = 12)
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

# Weekly terms (Period = 7)
build_fourier_weekly <- function(wday, P = 3, period = 7) {
  n <- length(wday)
  if (P == 0) return(as.data.frame(matrix(nrow = n, ncol = 0)))
  X <- matrix(NA, nrow = n, ncol = 2 * P)
  colnames(X) <- unlist(lapply(1:P, function(p) c(paste0("C_w", p), paste0("S_w", p))))
  for (p in 1:P) {
    X[, 2 * p - 1] <- cos(2 * pi * p * wday / period) 
    X[, 2 * p]     <- sin(2 * pi * p * wday / period) 
  }
  as.data.frame(X)
}

# ------------------------------------------------------------------------------
# 2. Model Fitting & Grid Search
# ------------------------------------------------------------------------------

# Fit combined LM with Trend and Fourier Series
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

# Grid Search for best Fourier Order (AIC)
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

# Wrapper to fit and return deseasonalized data
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

# ------------------------------------------------------------------------------
# 3. Future Forecasting / Reconstruction Helpers
# ------------------------------------------------------------------------------

# Build Design Matrix for Future Dates
future_X <- function(dates, P_year, P_month, P_week){
  # Ensure lubridate functions are available or use ::
  doy <- lubridate::yday(dates); doy[doy==366] <- 365
  Xy <- if(P_year>0) build_fourier_annual(doy, P=P_year, period=365) else as.data.frame(matrix(nrow=length(dates), ncol=0))
  Xm <- if(P_month>0) build_fourier_monthly(lubridate::month(dates), P=P_month, period=12) else as.data.frame(matrix(nrow=length(dates), ncol=0))
  Xw <- if(P_week>0) build_fourier_weekly(lubridate::wday(dates, week_start=1), P=P_week, period=7) else as.data.frame(matrix(nrow=length(dates), ncol=0))
  cbind(Xy, Xm, Xw)
}

# Calculate Seasonality Component for Future
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

# Calculate Trend Index for Future
add_trend_future <- function(lmfit, future_n){
  if(".tidx" %in% colnames(lmfit$model)){
    last_tidx <- max(lmfit$model$.tidx, na.rm = TRUE)
    return((last_tidx + seq_len(future_n)))
  } else return(rep(0, future_n))
}

# Calculate Trend Component (Value) for Future
add_trend_coef <- function(lmfit, trnd_vec){
  cfs <- coef(lmfit)
  if(".tidx" %in% names(cfs)){
    return(trnd_vec * cfs[".tidx"])
  } else return(rep(0, length(trnd_vec)))
}