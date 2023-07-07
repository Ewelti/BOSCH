#wd for Ellen
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

#Load packages
library(tidyverse)
library(plotrix)

# Read the data into a data frame
Yr_temps <- read.csv("output_data/BioData_linkedto_TempData.csv", header = TRUE, stringsAsFactors = FALSE)
Yr_temps <- Yr_temps %>% select(-1)

# Convert the Date column to the date format
Yr_temps$Date_new <- as.Date(Yr_temps$Date_new, format = "%m/%d/%Y")
Yr_temps$Year <- as.numeric(format(Yr_temps$Date_new, "%Y"))

# Plot the graph using ggplot2
ggplot(Yr_temps, aes(x = Date_new, y = Yryly_Temp, group = Site, color = Site)) +
  geom_line() +
  labs(x = "Date", y = "Temperature") +
  theme_minimal()

# Plot the graph using base plotting
svg("plots/temperature_vs_time.svg", width = 16, height = 3, pointsize = 12)
# tiff("plots/temperature_vs_time.tiff", width = 16, height = 4, units = "in", res = 600, pointsize = 12)
# Set up the plot area
plot(Yr_temps$Date_new, Yr_temps$Yryly_Temp, type = "n", xaxt = "n", xlab = "", ylab = "", las=1)

# Plot the Yr_temps for each site with lines
lines(Yr_temps$Date_new[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], col = "blue", lwd = 2)
lines(Yr_temps$Date_new[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], col = "red", lwd = 2)
lines(Yr_temps$Date_new[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], col = "green", lwd = 2)
lines(Yr_temps$Date_new[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], col = "orange", lwd = 2)

# Add solid circles as point markers
points(Yr_temps$Date_new[Yr_temps$Site == "Auba"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Auba"], col = "blue", pch = 16)
points(Yr_temps$Date_new[Yr_temps$Site == "Bieb"], Yr_temps$Yryly_Temp[Yr_temps$Site == "Bieb"], col = "red", pch = 16)
points(Yr_temps$Date_new[Yr_temps$Site == "O3"], Yr_temps$Yryly_Temp[Yr_temps$Site == "O3"], col = "green", pch = 16)
points(Yr_temps$Date_new[Yr_temps$Site == "W1"], Yr_temps$Yryly_Temp[Yr_temps$Site == "W1"], col = "orange", pch = 16)

# Find the first year of data for each site
min_year_Auba <- Yr_temps[24, 6]
min_year_Bieb <- Yr_temps[55, 6]
min_year_O3 <- Yr_temps[85, 6]
min_year_W1 <- Yr_temps[112, 6]

as.numeric(Yr_temps$Date_new[Yr_temps$Site == "Auba"])
# Add regression lines for each site starting from the first year with temperature Yr_temps
clip(x1=15180, x2=18082, y1 = -0, y2 = 100)
ablineclip(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Auba" & Date_new >= min_year_Auba)), col = "blue", lwd = 1, lty = 2)
abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "Bieb" & Date_new >= min_year_Bieb)), col = "red", lwd = 1, lty = 2)
abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "O3" & Date_new >= min_year_O3)), col = "green", lwd = 1, lty = 2)
abline(lm(Yryly_Temp ~ as.numeric(Date_new), data = subset(Yr_temps, Site == "W1" & Date_new >= min_year_W1)), col = "orange", lwd = 1, lty = 2)

# Add legend
legend("bottomleft", legend = c("Aubach", "Bieber", "Kinzig 1", "Kinzig 2"),
       col = c("blue", "red", "green", "orange"), pch = 16, bty="n")

# Format x-axis as years
axis.Date(1, at = seq(min(Yr_temps$Date_new), max(Yr_temps$Date_new), by = "year"), format = "%Y", las = 2)

# Add title and axis labels
title(main = "", xlab = "", ylab = expression("Temperature (" * degree * "C)"))

# finalise saved plot
dev.off()