# Codebase: Multivariate probabilistic modeling of electricity demand and generation

This repository contains the R code for the statistical analysis, modeling, and evaluation of electricity demand, solar power, and wind power generation. The focus of this Bachelor Thesis is the generation of multivariate scenarios (ensembles) considering dependency structures (Schaake Shuffle) and volatility clustering (GARCH).

## Project Structure

The project follows a modular structure to separate data preparation, modeling, and evaluation.

<pre> ```text
Code_BA/
│
├── main.R                     # Master Control Script (Orchestrator)
├── requirements.R             # Package installation script
├── README.md                  # Documentation (this file)
│
├── data/                      # Raw Data (CSV)
│   ├── energy_production.csv
│   ├── energy_demand.csv
│   └── installed_capacity.csv
│
├── functions/                 # Helper Functions
│   ├── eval_utils.R           # Scoring Rules (CRPS, ES), Brier Score
│   ├── fourier_utils.R        # Seasonality (Fourier Series)
│   ├── garch_utils.R          # GARCH Specification & Grid Search
│   └── schaake_utils.R        # Schaake Shuffle & Rank Correlation
│
├── scripts/                   # Analysis Modules
│   ├── 01_data_prep.R         # Data Import, Cleaning, Feature Engineering
│   ├── 02_modeling.R          # VAR-GARCH Modeling & Benchmarks
│   ├── 03_evaluation.R        # Scoring, DM-Tests, PIT-Histograms, Heatmaps
│   ├── 04_descriptive_analysis.R # Dark Doldrum Statistics & Time Series
│   └── 05_exploratory_plots.R    # ACF/PACF, Volatility Clustering, Stationarity
│
├── output/                    # Automatically generated results
│   ├── data_processed_all.RData   # Processed R objects
│   ├── models_and_ensembles.RData # Saved forecasts/ensembles
│   └── plots/                     # Generated figures (PDF/PNG)
│
└── 00_Archiv/
    ├── EntireAnalyis.R        # Development script that contains everything in one script (VAR-GARCH, Benchmarks, Scoring)
    ├── Plots&DescriptiveStats.R # Contains Calcluations and Plots for chapter 3 "Data"
    └── Dunkelflaute_Mockert.R # Dark Doldrum Analysis of complete data
``` </pre>   
~~~

## Installation & Prerequisites

### 1. R Version
The code was developed and tested using **R Version 4.5.1**.

### 2. Install Packages
To install all required dependencies automatically (CRAN & GitHub), run the `requirements.R` script once:

~~~r
source("requirements.R")
~~~

This installs necessary packages such as `tidyverse`, `rugarch`, `vars`, `forecast`, `scoringRules`, and `MultivCalibration` (via GitHub).

## Pipeline

The entire analysis is controlled centrally via the **`main.R`** file. You do not need to run the scripts in the `scripts/` folder individually.

### Steps:

1.  Open `main.R` in RStudio.
2.  Adjust the working directory path at the top of the file:
    ~~~r
    setwd('/Path/to/your/Code_BA')
    ~~~
3.  Use the **Flags (TRUE/FALSE)** to decide which parts of the pipeline to execute:

    ~~~r
    RUN_DATA_PREP   <- FALSE  # Step 1: Load/Clean data
    RUN_MODELING    <- FALSE  # Step 2: Fit models? (Time consuming!)
    RUN_EVAL        <- TRUE   # Step 3: Calculate scores & plots
    RUN_DESCRIPTIVE <- TRUE   # Step 4: Dark Doldrum statistics
    RUN_EXPLORATORY <- FALSE  # Step 5: Exploratory analysis
    ~~~

4.  Click **Run** to execute the pipeline.

## Module Descriptions

### 1. Data Preparation (`01_data_prep.R`)
* Imports raw CSV data.
* Aggregates wind power (Onshore + Offshore).
* Calculates Capacity Factors (CF) based on installed capacity.
* Identifies "Dark Doldrums" based on a specific threshold.
* Splits data into Training and Test sets (Split Date: 2022-12-31).

### 2. Modeling (`02_modeling.R`)
* **Seasonality:** Fourier Series decomposition (Trend + Seasonality).
* **VAR:** Vector Autoregression for modeling linear dependencies (Pre-Whitening).
* **GARCH:** Univariate GARCH models for residual volatility.
* **Schaake Shuffle:** Reconstruction of the multivariate dependency structure (Copula) based on historical ranks.
* **Benchmarks:** Generates an Analog Ensemble (AnEn) and a naive Climatology benchmark.

### 3. Evaluation (`03_evaluation.R`)
* Calculates multivariate scores (Energy Score), univariate scores (CRPS) and Brier Scores.
* Performs Diebold-Mariano tests for significance.
* Generates PIT Histograms (Probability Integral Transform) to check calibration.
* Visualizes results specifically for Dark Doldrums (Heatmaps, Improvement Plots).

### 4. Descriptive Analysis (`04_descriptive_analysis.R`)
* Visualizes the frequency of Dark Doldrums per year and month.
* Creates time series plots highlighting events below the threshold.

### 5. Exploratory Analysis (`05_exploratory_plots.R`)
* Analyzes statistical properties (Moments, Stationarity).
* Visualizes Volatility Clustering and Autocorrelation (ACF/PACF).
