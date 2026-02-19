# ==============================================================================
# File: main.R
# Purpose: Master script to orchestrate the entire analysis pipeline.
#          Controls execution of Data Prep, Modeling, and Evaluation.
# ==============================================================================

# 1. Project Setup -------------------------------------------------------------
project_path <- "[SET PROJECT PATH]"
setwd(project_path)

# 2. Execution Flags -----------------------------------------------------------
# Set these to TRUE or FALSE to control which parts of the pipeline to run.

RUN_DATA_PREP <- T  # Step 01: Load CSVs, clean data, split train/test
RUN_MODELING <- T  # Step 02: Fit VAR/GARCH, run Schaake Shuffle 
RUN_EVAL <- T  # Step 03: Calc Scores, create Plots & Tables
RUN_DESCRIPTIVE <- T  # Step 04: Descriptive Analysis
RUN_EXPLORATORY <- T  # Step 05: Exploratory Analysis

# 3. Execution Pipeline --------------------------------------------------------

# --- Step 1: Data Preparation ---
if (RUN_DATA_PREP) {
  cat("\n==================================================================\n")
  cat("Running Step 01: Data Preparation...\n")
  cat("==================================================================\n")
  source("scripts/01_data_prep.R", echo = TRUE)
} else {
  cat("\nSkipping Step 01 (Data Prep).\n")
}

# --- Step 2: Modeling (The Heavy Lifting) ---
if (RUN_MODELING) {
  cat("\n==================================================================\n")
  cat("Running Step 02: Modeling (VAR, GARCH, Benchmarks)...\n")
  cat("This may take a while due to GARCH grid search.\n")
  cat("==================================================================\n")
  source("scripts/02_modeling.R", echo = TRUE)
} else {
  cat("\nSkipping Step 02 (Modeling).\n")
}

# --- Step 3: Evaluation & Plots ---
if (RUN_EVAL) {
  cat("\n==================================================================\n")
  cat("Running Step 03: Evaluation & Visualization...\n")
  cat("==================================================================\n")
  source("scripts/03_evaluation.R", echo = TRUE)
} else {
  cat("\nSkipping Step 03 (Evaluation).\n")
}

# --- Step 4: Descriptive Analysis ---
if (RUN_DESCRIPTIVE) {
  cat("\n==================================================================\n")
  cat("Running Step 04: Descriptive Analysis...\n")
  cat("==================================================================\n")
  source("scripts/04_descriptive_analysis.R", echo = TRUE)
} else {
  cat("\nSkipping Step 04 (Descriptive Analysis).\n")
}

if (RUN_EXPLORATORY) {
  cat("\n==================================================================\n")
  cat("Running Step 05: Exploratory Analysis (Plots & Stats)...\n")
  cat("==================================================================\n")
  source("scripts/05_exploratory_plots.R", echo = TRUE)
} else {
  cat("\nSkipping Step 05 (Exploratory Analysis).\n")
}

cat("\n------------------------------------------------------------------\n")
cat("Pipeline finished successfully.\n")
cat("------------------------------------------------------------------\n")