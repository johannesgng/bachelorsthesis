# ==============================================================================
# File: 01_data_prep.R
# Purpose: Load raw CSV data, perform feature engineering, handle missing values,
#          calculate capacity factors, identify dark doldrums, and split data.
# Output:  Saves processed data objects to .RData files.
# ==============================================================================

# 1. Libraries -----------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)

# 2. Settings & Paths ----------------------------------------------------------
# Data Directory
data_dir <- 'data'
# Output Directory
output_dir <- 'output'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Parameter
threshold <- 0.06
split_date <- as.Date("2022-12-31")

# 3. Load Data -----------------------------------------------------------------

# Energy Production
energy_production <- read_delim(
  file.path(data_dir, "energy_production.csv"),
  delim = ";", 
  escape_double = FALSE, 
  col_types = cols(
    `Datum von` = col_date(format = "%Y-%m-%d"), 
    `Datum bis` = col_date(format = "%Y-%m-%d")
  ), 
  locale = locale(decimal_mark = ",", grouping_mark = "."), 
  trim_ws = TRUE
)

# Calculate Total Wind (Onshore + Offshore)
energy_production$Wind_MWh <- energy_production$`Wind Onshore [MWh] Berechnete Auflösungen` + 
                              energy_production$`Wind Offshore [MWh] Berechnete Auflösungen`

# Energy Demand
energy_demand <- read_delim(
  file.path(data_dir, "energy_demand.csv"),
  delim = ";", 
  escape_double = FALSE, 
  col_types = cols(
    `Datum von` = col_date(format = "%d.%m.%Y"), 
    `Datum bis` = col_date(format = "%d.%m.%Y")
  ), 
  locale = locale(decimal_mark = ",", grouping_mark = "."), 
  trim_ws = TRUE
)

# Installed Capacities
installed_cap_daily <- read_delim(
  file.path(data_dir, "installed_capacity.csv"),
  delim = ";", 
  escape_double = FALSE, 
  locale = locale(decimal_mark = ",", grouping_mark = "."), 
  trim_ws = TRUE
)

# 4. Feature Engineering & DataFrame Construction ------------------------------

energy <- data.frame(
  date = energy_production$`Datum von`,
  doy = yday(energy_production$`Datum von`), 
  month = month(energy_production$`Datum von`),
  year = year(energy_production$`Datum von`),
  wday = wday(energy_production$`Datum von`, week_start = 1),
  calendar_week = week(energy_production$`Datum von`),
  demand = energy_demand$`Netzlast [MWh] Berechnete Auflösungen`,
  solar = energy_production$`Photovoltaik [MWh] Berechnete Auflösungen`,
  wind = energy_production$Wind_MWh,
  residual_load = energy_demand$`Residuallast [MWh] Berechnete Auflösungen`,
  
  # Capacity Data
  installed_solar = installed_cap_daily$`Photovoltaik [MW] Berechnete Auflösungen`,
  installed_wind = rowSums(cbind(
    installed_cap_daily$`Wind Offshore [MW] Berechnete Auflösungen`, 
    installed_cap_daily$`Wind Onshore [MW] Berechnete Auflösungen`
  ), na.rm = TRUE)
) %>% 
  as_tibble() %>%
  mutate(
    # Max Generation (MW -> MWh per day)
    P_max_solar = installed_solar * 24,
    P_max_wind  = installed_wind * 24,
    P_max_total = P_max_solar + P_max_wind,
    
    # Capacity Factors
    CF_solar = solar / P_max_solar,
    CF_wind  = wind / P_max_wind,
    CF_total = (solar + wind) / P_max_total,
    
    # Log Transformations
    demand_log   = log(demand),
    solar_log    = log(solar),
    wind_log     = log(wind),
    CF_solar_log = log(CF_solar),
    CF_wind_log  = log(CF_wind),
    
    # Dark Doldrum Identification
    h48_mean = rollmean(CF_total, k = 2, fill = NA, align = "right"),
    dark_doldrum = if_else(h48_mean < threshold, 1, 0)
  )

# 5. Split Train / Test --------------------------------------------------------

energy_train <- subset(energy, date <= split_date)
energy_test <- subset(energy, date > split_date)

# 6. Save Data -----------------------------------------------------------------

# Save individual files
save(energy, file = file.path(output_dir, 'energy.RData'))
save(energy_train, file = file.path(output_dir, 'energy_train.RData'))
save(energy_test, file = file.path(output_dir, 'energy_test.RData'))

# Save all in one file
save(energy, energy_train, energy_test, threshold, split_date, 
     file = file.path(output_dir, 'data_processed_all.RData'))

cat("Data preparation complete. Files saved to:", output_dir, "\n")