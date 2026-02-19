library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

# Check if you adjusted the project directory correctly in main.R
load("00_Archiv/energy.RData")
load("00_Archiv/energy_train.RData")
load("00_Archiv/energy_test.RData")

summary(energy$CF_total)
threshold <- 0.06

plot(energy$date, energy$h48_mean, type = "l", col = "black", lwd = 2, xlab = "Date", ylab = "48h Mean CF")
abline(h = threshold, col = "red", lwd = 2)
usr <- par("usr") 
clip(usr[1], usr[2], usr[3], threshold)
lines(energy$date, energy$h48_mean, col = "#009682", lwd = 2)
do.call("clip", as.list(usr))

print(paste("Number of Dark Doldrums in Test Set:", sum(energy$dark_doldrum[-1]))) # ignores NA

events_per_month <- energy %>% group_by(month) %>% summarise(number_of_events = sum(dark_doldrum == 1 & lag(dark_doldrum, default = 0) == 0, na.rm = T))
events_per_month$month <- month.name[events_per_month$month]
events_per_month$month <- factor(events_per_month$month, levels = month.name)
barplot(
  height = events_per_month$number_of_events,
  names.arg = events_per_month$month,
  col = "#009682",
  xlab = "Month",
  ylab = "Number of Events",
  main = "Number of Events per Month"
)

ggplot(events_per_month, aes(x = month, y = number_of_events)) +
  geom_col(fill = "#009682", color = "black") +
  labs(
    title = "Number of Events per Month",
    x = "Month",
    y = "Number of Events"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


events_per_year <- energy %>% group_by(year) %>% summarise(number_of_events = sum(dark_doldrum == 1 & lag(dark_doldrum, default = 0) == 0, na.rm = T)) 
events_per_year$year <- as.character(events_per_year$year)
events_per_year$year[events_per_year$year == "2025"] <- "2025 (ongoing)"

barplot(
  height = events_per_year$number_of_events,
  names.arg = events_per_year$year,
  col = "#009682",
  xlab = "Year",
  ylab = "Number of Events",
  main = "Number of Events per Year",
  cex.names = 0.672
)

ggplot(events_per_year, aes(x = year, y = number_of_events)) +
  geom_col(fill = "#009682", color = "black") +
  labs(
    title = "Number of Events per Year",
    x = "Year",
    y = "Number of Events"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

