# ==============================================================================
# File: schaake_utils.R
# Purpose: Functions for constructing the empirical copula, generating 
#          ensembles via ECC-Q (Quantization), and applying the Schaake Shuffle.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Rank Matrix Construction
# ------------------------------------------------------------------------------

get_rankmatrix <- function(data, param_dependence = F) { 
  # Constructs the rank structure from historical data (Empirical Copula)
  # param_dependence is kept for compatibility, default uses empirical ranks
  
  # Note: MARGIN = 2 => Function is applied column-wise
  ranks <- apply(data, MARGIN = 2, rank, ties.method = "first") 
  
  # Implied sorting of the values for each column and interchange value 
  # with index from the sorted order
  
  return(ranks) 
}

# ------------------------------------------------------------------------------
# 2. Empirical Distribution Sampling (ECC-Q)
# ------------------------------------------------------------------------------

r_empirical_dist <- function(M, mu, sigma, z_hist) { 
  # Re-implementation of 'density_forecast_nonparam' from Python code
  # Generates M ensemble members based on equidistant quantiles of the 
  # historical standardized residuals (z_hist)
  
  z_sorted <- sort(z_hist) # sorting the standardized residuals
  n_hist   <- length(z_sorted)
  
  # Equidistant Quantiles
  probs <- (1:M) / (M + 1)
  
  z_quant <- approx(
    x = (1:n_hist) / (n_hist + 1),
    y = z_sorted,
    xout = probs,
    rule = 2
  )$y
  
  # Scale by predicted mean and volatility
  y_samples <- mu + sigma * z_quant
  
  return(as.numeric(y_samples))
}

# ------------------------------------------------------------------------------
# 3. The Schaake Shuffle
# ------------------------------------------------------------------------------

schaake_shuffle <- function(forecast_mat, rank_mat) {
  # Reorders the ensemble members (forecast_mat) to match the rank structure 
  # provided by rank_mat.
  
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