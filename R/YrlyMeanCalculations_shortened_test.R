library(tidyverse)

## Load data
Bio_dates <- read.csv("RawData/Temp_data/RMO_BioDataDates.csv", header = TRUE, sep = ",")
Temp <- read.csv("RawData/Temp_data/Interpolated_dailyTemp.csv", header = TRUE, sep = ",")

## Ensure date columns in both datasets are correct
Bio_dates$Date_new <- as.Date(Bio_dates$Date, format = "%m/%d/%Y")
Temp$Date_new <- as.Date(Temp$Date, format = "%m/%d/%Y")

## Function to calculate average temperature and missing days
calculate_temp_missing_days <- function(Bio_dates_subset, site_name, gap_col) {
  for (i in 1:nrow(Bio_dates_subset)) {
    bio_sampling_date <- as.Date(Bio_dates_subset$Date_new[i])
    start_date <- bio_sampling_date - days(1) - months(12)
    
    filtered_data <- Temp %>%
      filter(Date_new > start_date & Date_new <= bio_sampling_date) %>%
      slice(-n())
    
    num_entries <- sum(!is.na(filtered_data[[gap_col]]))
    missing_days <- 365 - num_entries
    
    if (missing_days <= 5) {
      average_temperature <- mean(filtered_data[[gap_col]], na.rm = TRUE)
    } else {
      average_temperature <- NA
    }
    
    Bio_dates_subset$Yryly_Temp[i] <- average_temperature
    Bio_dates_subset$Missing_Days[i] <- missing_days
  }
  return(Bio_dates_subset)
}

## Site names and corresponding gap columns
sites <- c("Auba", "Bieb", "O3", "W1")
gap_columns <- c("au_gap", "bieb_gap", "kio_gap", "kiw_gap")

## Calculate temperature and missing days for each site
BioData_linkedto_TempData <- lapply(seq_along(sites), function(i) {
  site_name <- sites[i]
  gap_col <- gap_columns[i]
  Bio_dates_subset <- Bio_dates[grepl(site_name, Bio_dates$Site), ]
  Bio_dates_subset <- Bio_dates_subset[order(Bio_dates_subset$Date_new), ]
  calculate_temp_missing_days(Bio_dates_subset, site_name, gap_col)
})

## Combine everything into a big table for exporting
BioData_linkedto_TempData2 <- do.call(rbind, BioData_linkedto_TempData)
write.csv(BioData_linkedto_TempData, "output_data/BioData_linkedto_TempData.csv", row.names = FALSE)

##### CLEAN UP --------------------
library(pacman)
# Clear data
rm(list = ls())  # Removes all objects from environment
# Clear packages
p_unload(all)  # Remove all contributed packages
# Clear plots
graphics.off()  # Clears plots, closes all graphics devices
# Clear console
cat("\014")  # Mimics ctrl+L
# Clear mind :)
