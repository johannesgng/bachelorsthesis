# ==============================================================================
# File: requirements.R
# Purpose: Install and load all required R packages for the project.
#          Ensures reproducibility across different environments.
# ==============================================================================

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste("Package already installed:", pkg))
  }
}

# 1. CRAN Packages -------------------------------------------------------------

cran_packages <- c(
  # Data Manipulation & Utils
  "tidyverse",      # Includes dplyr, tidyr, readr, ggplot2, etc.
  "lubridate",      # Date handling
  "zoo",            # Rolling windows
  "knitr",          # Reporting
  "devtools",       # GitHub Download
  
  # Modeling & Time Series
  "forecast",       # ARIMA, DM-Test
  "vars",           # Vector Autoregression
  "urca",           # Cointegration / Unit Root Tests (used in VAR)
  "rugarch",        # GARCH Modeling
  "tseries",        # ADF / KPSS Tests
  "FinTS",          # ARCH Test
  
  # Evaluation & Scoring
  "scoringRules",   # CRPS, Energy Score
  
  # Visualization & Exploration
  "ggplot2",        # Plotting
  "RColorBrewer",   # Color Palettes
  "viridis",        # Color Palettes (Exploratory plots)
  "corrplot",       # Correlation Plots
  "moments"         # Skewness/Kurtosis
)

invisible(lapply(cran_packages, install_if_missing))

# 2. GitHub Packages -----------------------------------------------------------

# 'devtools' is required to install packages from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 'MultivCalibration' is needed for Rank Histogram
if (!requireNamespace("MultivCalibration", quietly = TRUE)) {
  message("Installing 'MultivCalibration' from GitHub")
  devtools::install_github("sallen12/MultivCalibration")
}

# 4. Load Packages (Sanity Check) ----------------------------------------------
all_packages <- c(cran_packages, "MultivCalibration")

loaded <- lapply(all_packages, require, character.only = TRUE)
success <- all(unlist(loaded))

if (success) {
  cat("\n==================================================================\n")
  cat(" SUCCESS: All required packages are installed and loaded.\n")
  cat("==================================================================\n")
} else {
  cat("\n==================================================================\n")
  cat(" ERROR: Some packages could not be loaded. Please check logs.\n")
  cat("==================================================================\n")
}
