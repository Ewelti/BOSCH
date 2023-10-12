# wd for Ellen
# setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

#Load packages
library(tidyverse)
library(plotrix)

# Read the data into a data frame
Yr_temps <- read.csv("RawData/Temp_data/BioData_linkedto_TempData.csv", header = TRUE, stringsAsFactors = FALSE)
Yr_temps <- Yr_temps %>% select(-1)

# Convert the Date column to the date format
Yr_temps$Date_new <- as.Date(Yr_temps$Date_new, format = "%Y-%m-%d")
Yr_temps$Year <- as.numeric(format(Yr_temps$Date_new, "%Y"))

# Plot the graph using base plotting
svg("plots/temperature_vs_time_updated.svg", width = 16, height = 4, pointsize = 12)
tiff("plots/temperature_vs_time_updated.tiff", width = 16, height = 4, units = "in", res = 600, pointsize = 12)
# Set up the plot area
plot(Yr_temps$Date_new, Yr_temps$Yryly_Temp, ylim = c(8.5, 12.8), type = "n", xaxt = "n", xlab = "", ylab = "", las=1)

# Plot the Yr_temps for each site with lines
# lines(Yr_temps$Date_new[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], col = "blue", lwd = 2)
# lines(Yr_temps$Date_new[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], col = "red", lwd = 2)
# lines(Yr_temps$Date_new[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], col = "green", lwd = 2)
# lines(Yr_temps$Date_new[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], col = "orange", lwd = 2)
lines(Yr_temps$Date_new[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], lwd = 2, col = alpha(1, 0.6), lty = 1)
lines(Yr_temps$Date_new[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], lwd = 2, col = alpha(2, 0.6), lty = 1)
lines(Yr_temps$Date_new[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], lwd = 2, col = alpha(3, 0.6), lty = 1)
lines(Yr_temps$Date_new[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], lwd = 2, col = alpha(5, 0.6), lty = 1)

# Add solid circles as point markers
# points(Yr_temps$Date_new[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], col = "blue", pch = 15)
# points(Yr_temps$Date_new[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], col = "red", pch = 16)
# points(Yr_temps$Date_new[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], col = "green", pch = 17)
# points(Yr_temps$Date_new[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], col = "orange", pch = 18)
points(Yr_temps$Date_new[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], pch = 21, bg = alpha(1, 0.6), col = alpha(1, 0.6), lwd = 2, cex = 2.5)
points(Yr_temps$Date_new[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], pch = 22, bg = alpha(2, 0.6), col = alpha(2, 0.6), lwd = 2, cex = 2.5)
points(Yr_temps$Date_new[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], pch = 23, bg = alpha(3, 0.6),col = alpha(3, 0.6), lwd = 2, cex = 2.5)
points(Yr_temps$Date_new[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], pch = 24, bg = alpha(4, 0.6),col = alpha(4, 0.6), lwd = 2, cex = 2.5)

# Find the first year of data for each site
min_year_Auba <- Yr_temps[24, 6]
min_year_Bieb <- Yr_temps[55, 6]
min_year_O3 <- Yr_temps[85, 6]
min_year_W1 <- Yr_temps[112, 6]

# Clip ablines to temp observation start years
as.numeric(Yr_temps$Date_new[Yr_temps$Site == "Auba"])
clip(x1=15180, x2=18082, y1 = -0, y2 = 100)

# Add regression lines for each site starting from the first year with temperature Yr_temps
# ablineclip(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Auba" & Date_new >= min_year_Auba)), col = "blue", lwd = 1, lty = 2)
# abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Bieb" & Date_new >= min_year_Bieb)), col = "red", lwd = 1, lty = 2)
# abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "O3" & Date_new >= min_year_O3)), col = "green", lwd = 1, lty = 2)
# abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "W1" & Date_new >= min_year_W1)), col = "orange", lwd = 1, lty = 2)
ablineclip(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Auba" & Date_new >= min_year_Auba)), lwd = 1, col = alpha(1, 0.6), lty = 2)
abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Bieb" & Date_new >= min_year_Bieb)), lwd = 1, col = alpha(2, 0.6), lty = 2)
abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "O3" & Date_new >= min_year_O3)), lwd = 1, col = alpha(3, 0.6), lty = 2)
abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "W1" & Date_new >= min_year_W1)), lwd = 1, col = alpha(4, 0.6), lty = 2)

# Format x-axis as years
axis.Date(1, at = seq(min(Yr_temps$Date_new), max(Yr_temps$Date_new), by = "year"), format = "%Y", las = 2)

# Add title and axis labels
title(xlab = "Observation period", ylab = expression("Temperature (" * degree * "C)"))

# Add legend
legend("bottomleft", legend = c("Aubach", "Bieber", "Kinzig 1", "Kinzig 2"),
       bg = c(alpha(1, 0.6), alpha(2, 0.6), alpha(3, 0.6), alpha(4, 0.6)), 
       col = c(alpha(1, 0.6), alpha(2, 0.6), alpha(3, 0.6), alpha(4, 0.6)), 
       pch = c(21:24), bty = "n", cex = 2, lwd = 2)

# finalise saved plot
dev.off()

# extract linear regression results
Auba_lm <- lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Auba" & Date_new >= min_year_Auba))
summary(Auba_lm)
Bieb_lm <- lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Bieb" & Date_new >= min_year_Bieb))
summary(Bieb_lm)
O3_lm <- lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "O3" & Date_new >= min_year_O3))
summary(O3_lm)
W1_lm <- lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "W1" & Date_new >= min_year_W1))
summary(W1_lm)

# Calculate percentage change through time
# Auba
## annual change in degrees C through time
Auba <- subset(Yr_temps, Site == "Auba" & Date_new >= min_year_Auba) # potentially remove NAs: & !is.na(Yryly_Temp)
Auba_est <- as.numeric(coef(Auba_lm)["as.numeric(Date_new)"])
Auba_ave <- mean(Auba$Yryly_Temp, na.rm = T)
Auba_dCChange_perYr <- (Auba_est)*Auba_ave
print(Auba_dCChange_perYr)
## % change through time
Auba <- subset(Yr_temps, Site == "Auba" & Date_new >= min_year_Auba) # potentially remove NAs: & !is.na(Yryly_Temp)
Auba_est <- as.numeric(coef(Auba_lm)["as.numeric(Date_new)"])
Auba_ave <- mean(Auba$Yryly_Temp, na.rm = T)
Auba_percChange_perYr <- (Auba_est/Auba_ave)*100
print(Auba_percChange_perYr)
# Bieb
## annual change in degrees C through time
Bieb <- subset(Yr_temps, Site == "Bieb" & Date_new >= min_year_Bieb) # potentially remove NAs: & !is.na(Yryly_Temp)
Bieb_est <- as.numeric(coef(Bieb_lm)["as.numeric(Date_new)"])
Bieb_ave <- mean(Bieb$Yryly_Temp, na.rm = T)
Bieb_dCChange_perYr <- (Bieb_est)*Bieb_ave
print(Bieb_dCChange_perYr)
## % change through time
Bieb <- subset(Yr_temps, Site == "Bieb" & Date_new >= min_year_Bieb) # potentially remove NAs: & !is.na(Yryly_Temp)
Bieb_est <- as.numeric(coef(Bieb_lm)["as.numeric(Date_new)"])
Bieb_ave <- mean(Bieb$Yryly_Temp, na.rm = T)
Bieb_percChange_perYr <- (Bieb_est/Bieb_ave)*100
print(Bieb_percChange_perYr)
# O3
## annual change in degrees C through time
O3 <- subset(Yr_temps, Site == "O3" & Date_new >= min_year_O3) # potentially remove NAs: & !is.na(Yryly_Temp)
O3_est <- as.numeric(coef(O3_lm)["as.numeric(Date_new)"])
O3_ave <- mean(O3$Yryly_Temp, na.rm = T)
O3_dCChange_perYr <- (O3_est)*O3_ave
print(O3_dCChange_perYr)
# % change through time
O3 <- subset(Yr_temps, Site == "O3" & Date_new >= min_year_O3) # potentially remove NAs: & !is.na(Yryly_Temp)
O3_est <- as.numeric(coef(O3_lm)["as.numeric(Date_new)"])
O3_ave <- mean(O3$Yryly_Temp, na.rm = T)
O3_percChange_perYr <- (O3_est/O3_ave)*100
print(O3_percChange_perYr)
# W1
## annual change in degrees C through time
W1 <- subset(Yr_temps, Site == "W1" & Date_new >= min_year_W1) # potentially remove NAs: & !is.na(Yryly_Temp)
W1_est <- as.numeric(coef(W1_lm)["as.numeric(Date_new)"])
W1_ave <- mean(W1$Yryly_Temp, na.rm = T)
W1_dCChange_perYr <- (W1_est)*W1_ave
print(W1_dCChange_perYr)
## % change through time
W1 <- subset(Yr_temps, Site == "W1" & Date_new >= min_year_W1) # potentially remove NAs: & !is.na(Yryly_Temp)
W1_est <- as.numeric(coef(W1_lm)["as.numeric(Date_new)"])
W1_ave <- mean(W1$Yryly_Temp, na.rm = T)
W1_percChange_perYr <- (W1_est/W1_ave)*100
print(W1_percChange_perYr)
