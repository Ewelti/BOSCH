# wd for Ellen
# setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

#Load packages
library(tidyverse)
library(plotrix)
library(dplyr)

# Read the data into a data frame
Yr_temps <- read.csv("RawData/Temp_data/BioData_linkedto_TempData_updatedR1.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert the Date column to the date format
Yr_temps$Date_new <- as.Date(Yr_temps$Date_new, format = "%Y-%m-%d")
Yr_temps$Year <- as.numeric(format(Yr_temps$Date_new, "%Y"))

# Plot the graph using base plotting
# svg("plots/temperature_vs_time_updatedR1.svg", width = 7, height = 3.1, pointsize = 12)
tiff("plots/temperature_vs_time_updatedR1.tiff", width = 10, height = 5, units = "in", res = 600, pointsize = 12, compression = 'lzw')

# Reduce margins to minimize white space (bottom, left, top, right)
par(mar = c(5, 4, 1, 1))

# Set up the plot area
plot(Yr_temps$Year, Yr_temps$Yryly_Temp, ylim = c(8.5, 12.8), type = "n", xaxt = "n", xlab = "", ylab = "", las=1)

# Plot the temperatures for each site with lines (FIXED: consistent colors)
lines(Yr_temps$Year[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], lwd = 2, col = alpha(1, 0.6), lty = 1)
lines(Yr_temps$Year[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], lwd = 2, col = alpha(2, 0.6), lty = 1)
lines(Yr_temps$Year[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], lwd = 2, col = alpha(3, 0.6), lty = 1) # significant
lines(Yr_temps$Year[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], lwd = 2, col = alpha(4, 0.6), lty = 1) # significant

# Add solid circles as point markers (FIXED: consistent colors)
points(Yr_temps$Year[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], pch = 21, bg = alpha(1, 0.6), col = alpha(1, 0.6), lwd = 2, cex = 1)
points(Yr_temps$Year[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], pch = 22, bg = alpha(2, 0.6), col = alpha(2, 0.6), lwd = 2, cex = 1)
points(Yr_temps$Year[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], pch = 23, bg = alpha(3, 0.6),col = alpha(3, 0.6), lwd = 2, cex = 1)
points(Yr_temps$Year[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], pch = 24, bg = alpha(4, 0.6),col = alpha(4, 0.6), lwd = 2, cex = 1)

# Find the year range for each site and convert to dates for plotting
min_year_Auba <- Yr_temps |> filter(Site == "Auba") |> pull(Year) |> min(na.rm = TRUE)
max_year_Auba <- Yr_temps |> filter(Site == "Auba") |> pull(Year) |> max(na.rm = TRUE)
min_year_Bieb <- Yr_temps |> filter(Site == "Bieb") |> pull(Year) |> min(na.rm = TRUE)
max_year_Bieb <- Yr_temps |> filter(Site == "Bieb") |> pull(Year) |> max(na.rm = TRUE)
min_year_O3 <- Yr_temps |> filter(Site == "O3") |> pull(Year) |> min(na.rm = TRUE)
max_year_O3 <- Yr_temps |> filter(Site == "O3") |> pull(Year) |> max(na.rm = TRUE)
min_year_W1 <- Yr_temps |> filter(Site == "W1") |> pull(Year) |> min(na.rm = TRUE)
max_year_W1 <- Yr_temps |> filter(Site == "W1") |> pull(Year) |> max(na.rm = TRUE)

# Add regression lines for each site clipped to their respective date ranges
# Note: ablineclip needs Date_new for x-axis alignment, but we'll use Year-based models for analysis
ablineclip(lm(Yryly_Temp ~ as.numeric(Year), data = subset(Yr_temps, Site == "Auba" & Year >= min_year_Auba)), x1 = min_year_Auba, x2 = max_year_Auba, lwd = 1.5, col = alpha(1, 0.6), lty = 2)
ablineclip(lm(Yryly_Temp ~ as.numeric(Year), data = subset(Yr_temps, Site == "Bieb" & Year >= min_year_Bieb)), x1 = min_year_Bieb, x2 = max_year_Bieb, lwd = 1.5, col = alpha(2, 0.6), lty = 2)
ablineclip(lm(Yryly_Temp ~ as.numeric(Year), data = subset(Yr_temps, Site == "O3" & Year >= min_year_O3)), x1 = min_year_O3, x2 = max_year_O3, lwd = 2, col = alpha(3, 0.6), lty = 1) # significant
ablineclip(lm(Yryly_Temp ~ as.numeric(Year), data = subset(Yr_temps, Site == "W1" & Year >= min_year_W1)), x1 = min_year_W1, x2 = max_year_W1, lwd = 2, col = alpha(4, 0.6), lty = 1) # significant

# Format x-axis as years
# axis.Date(1, at = seq(min(Yr_temps$Date_new), max(Yr_temps$Date_new), by = "year"), format = "%Y", las = 2)
axis(1, at = seq(min(Yr_temps$Year), max(Yr_temps$Year), by = 1), las = 2)

# Add title and axis labels
title(xlab = "Year", ylab = expression("Temperature (" * degree * "C)"))

# Add legend
legend("topleft", legend = c("Aubach", "Bieber", "Kinzig 1", "Kinzig 2"),
       bg = c(alpha(1, 0.6), alpha(2, 0.6), alpha(3, 0.6), alpha(4, 0.6)),
       col = c(alpha(1, 0.6), alpha(2, 0.6), alpha(3, 0.6), alpha(4, 0.6)),
       pch = c(21:24), bty = "n", cex = 1.5, lwd = 2)

# finalise saved plot
dev.off()

# Extract linear regression results (using Year as predictor)
Auba_lm <- lm(Yryly_Temp ~ Year, data = subset(Yr_temps, Site == "Auba" & Year >= min_year_Auba))
summary(Auba_lm)
Bieb_lm <- lm(Yryly_Temp ~ Year, data = subset(Yr_temps, Site == "Bieb" & Year >= min_year_Bieb))
summary(Bieb_lm)
O3_lm <- lm(Yryly_Temp ~ Year, data = subset(Yr_temps, Site == "O3" & Year >= min_year_O3))
summary(O3_lm)
W1_lm <- lm(Yryly_Temp ~ Year, data = subset(Yr_temps, Site == "W1" & Year >= min_year_W1))
summary(W1_lm)

# Calculate temperature changes through time (using Year as predictor)
# Auba
Auba <- subset(Yr_temps, Site == "Auba" & Year >= min_year_Auba)
Auba_est <- as.numeric(coef(Auba_lm)["Year"])  # Annual slope coefficient (°C/year)
Auba_ave <- mean(Auba$Yryly_Temp, na.rm = TRUE)
Auba_dCChange_perYr <- Auba_est  # Annual change in °C/year (direct from slope)
Auba_percChange_perYr <- (Auba_dCChange_perYr/Auba_ave) * 100  # Annual % change
# Total change from first to last year
Auba_years_span <- max_year_Auba - min_year_Auba
Auba_total_dCChange <- Auba_est * Auba_years_span  # Total °C change over period
Auba_total_percChange <- (Auba_total_dCChange/Auba_ave) * 100  # Total % change over period
print(paste("Auba annual change:", round(Auba_dCChange_perYr, 4), "°C/year"))
print(paste("Auba annual % change:", round(Auba_percChange_perYr, 4), "%/year"))
print(paste("Auba total change (", min_year_Auba, "-", max_year_Auba, "):", round(Auba_total_dCChange, 4), "°C"))
print(paste("Auba total % change (", min_year_Auba, "-", max_year_Auba, "):", round(Auba_total_percChange, 4), "%"))

# Bieb
Bieb <- subset(Yr_temps, Site == "Bieb" & Year >= min_year_Bieb)
Bieb_est <- as.numeric(coef(Bieb_lm)["Year"])
Bieb_ave <- mean(Bieb$Yryly_Temp, na.rm = TRUE)
Bieb_dCChange_perYr <- Bieb_est
Bieb_percChange_perYr <- (Bieb_dCChange_perYr/Bieb_ave) * 100
# Total change from first to last year
Bieb_years_span <- max_year_Bieb - min_year_Bieb
Bieb_total_dCChange <- Bieb_est * Bieb_years_span
Bieb_total_percChange <- (Bieb_total_dCChange/Bieb_ave) * 100
print(paste("Bieb annual change:", round(Bieb_dCChange_perYr, 4), "°C/year"))
print(paste("Bieb annual % change:", round(Bieb_percChange_perYr, 4), "%/year"))
print(paste("Bieb total change (", min_year_Bieb, "-", max_year_Bieb, "):", round(Bieb_total_dCChange, 4), "°C"))
print(paste("Bieb total % change (", min_year_Bieb, "-", max_year_Bieb, "):", round(Bieb_total_percChange, 4), "%"))

# O3
O3 <- subset(Yr_temps, Site == "O3" & Year >= min_year_O3)
O3_est <- as.numeric(coef(O3_lm)["Year"])
O3_ave <- mean(O3$Yryly_Temp, na.rm = TRUE)
O3_dCChange_perYr <- O3_est
O3_percChange_perYr <- (O3_dCChange_perYr/O3_ave) * 100
# Total change from first to last year
O3_years_span <- max_year_O3 - min_year_O3
O3_total_dCChange <- O3_est * O3_years_span
O3_total_percChange <- (O3_total_dCChange/O3_ave) * 100
print(paste("O3 annual change:", round(O3_dCChange_perYr, 4), "°C/year"))
print(paste("O3 annual % change:", round(O3_percChange_perYr, 4), "%/year"))
print(paste("O3 total change (", min_year_O3, "-", max_year_O3, "):", round(O3_total_dCChange, 4), "°C"))
print(paste("O3 total % change (", min_year_O3, "-", max_year_O3, "):", round(O3_total_percChange, 4), "%"))

# W1
W1 <- subset(Yr_temps, Site == "W1" & Year >= min_year_W1)
W1_est <- as.numeric(coef(W1_lm)["Year"])
W1_ave <- mean(W1$Yryly_Temp, na.rm = TRUE)
W1_dCChange_perYr <- W1_est
W1_percChange_perYr <- (W1_dCChange_perYr/W1_ave) * 100
# Total change from first to last year
W1_years_span <- max_year_W1 - min_year_W1
W1_total_dCChange <- W1_est * W1_years_span
W1_total_percChange <- (W1_total_dCChange/W1_ave) * 100
print(paste("W1 annual change:", round(W1_dCChange_perYr, 4), "°C/year"))
print(paste("W1 annual % change:", round(W1_percChange_perYr, 4), "%/year"))
print(paste("W1 total change (", min_year_W1, "-", max_year_W1, "):", round(W1_total_dCChange, 4), "°C"))
print(paste("W1 total % change (", min_year_W1, "-", max_year_W1, "):", round(W1_total_percChange, 4), "%"))

