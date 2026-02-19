# ==============================================================================
# File: 04_descriptive_analysis.R
# Purpose: Descriptive statistics and visualization of Dark Doldrum events.
#          - Time Series Plot with Threshold Clipping
#          - Event Frequency per Month and Year
# Output:  Saves plots to 'output/plots/'
# ==============================================================================

# 1. Setup & Libraries ---------------------------------------------------------
library(dplyr)
library(lubridate)
library(ggplot2)

# Define Output Directory for Plots
plot_dir <- "output/plots"


# 2. Load Data -----------------------------------------------------------------
cat("Loading processed data...\n")
load("output/data_processed_all.RData") 


# 3. Time Series Plot (Base R with Clipping) -----------------------------------
cat("Generating Time Series Plot (Threshold Clipping)...\n")

# Plot Setup
plot(energy$date, energy$h48_mean, type = "l", col = "black", lwd = 1.5, 
     xlab = "Date", ylab = "48h Mean Capacity Factor",
     main = "Dark Doldrum Events (Below Threshold)", ylim = c(0, 1))

# Threshold Line
abline(h = threshold, col = "red", lwd = 2, lty = 2)

# Clipping Logic
usr <- par("usr") 
clip(usr[1], usr[2], usr[3], threshold)
lines(energy$date, energy$h48_mean, col = "#009682", lwd = 2)
do.call("clip", as.list(usr)) # Reset clipping

# Legend
legend("topright", legend = c("CF > Threshold", "Dark Doldrum (CF < Threshold)", "Threshold"),
       col = c("black", "#009682", "red"), lty = c(1, 1, 2), lwd = c(1.5, 2, 2))

# 4. Event Statistics ----------------------------------------------------------
cat("Calculating Event Statistics...\n")

# Count total moments (ignoring NA)
cat("Number of time steps classified as Dark Doldrums:", sum(energy$dark_doldrum, na.rm = TRUE), "\n")

# Identify distinct EVENTS (consecutive days count as 1 event)
# Logic: A day is a start of an event if it is 1 AND the day before was 0
energy_events <- energy %>%
  arrange(date) %>%
  mutate(
    is_event_start = ifelse(dark_doldrum == 1 & lag(dark_doldrum, default = 0) == 0, 1, 0)
  )

# --- A) Events per Month ---
events_per_month <- energy_events %>% 
  group_by(month) %>% 
  summarise(number_of_events = sum(is_event_start, na.rm = TRUE)) %>%
  mutate(
    month_name = factor(month.name[month], levels = month.name)
  )

p_month <- ggplot(events_per_month, aes(x = month_name, y = number_of_events)) +
  geom_col(fill = "#009682", color = "black", width = 0.7) +
  labs(
    title = "Frequency of Dark Doldrums per Month",
    subtitle = paste("Threshold CF <", threshold),
    x = "Month",
    y = "Number of Distinct Events"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_month)
ggsave(file.path(plot_dir, "barplot_events_per_month.pdf"), p_month, width = 8, height = 5, dpi = 300)


# --- B) Events per Year ---
events_per_year <- energy_events %>% 
  group_by(year) %>% 
  summarise(number_of_events = sum(is_event_start, na.rm = TRUE)) %>%
  mutate(
    year_label = as.character(year),
    year_label = ifelse(year == 2025, "2025 (ongoing)", year_label)
  )

p_year <- ggplot(events_per_year, aes(x = year_label, y = number_of_events)) +
  geom_col(fill = "#009682", color = "black", width = 0.6) +
  labs(
    title = "Frequency of Dark Doldrums per Year",
    x = "Year",
    y = "Number of Distinct Events"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_year)
ggsave(file.path(plot_dir, "barplot_events_per_year.pdf"), p_year, width = 8, height = 5, dpi = 300)
