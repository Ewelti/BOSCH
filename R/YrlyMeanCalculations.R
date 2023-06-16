library(tidyverse)
library(dplyr)

## Load data
# Biological data sampling dates
Bio_dates <- read.csv("RawData/Temp_data/RMO_BioDataDates.csv", header = T, sep = ",")

# interpolated temperature data - daily means
Temp <- read.csv("RawData/Temp_data/Interpolated_dailyTemp.csv", header = T, sep = ",")

## Ensure date columns in both datasets are correct
Bio_dates$Date_new <- as.Date(Bio_dates$Date, format = "%m/%d/%Y")
Temp$Date_new <- as.Date(Temp$Date, format = "%m/%d/%Y")

## make subsets for each site
Bio_dates_auba <- Bio_dates[grepl("Auba", Bio_dates$Site), ]
Bio_dates_bieb <- Bio_dates[grepl("Bieb", Bio_dates$Site), ]
Bio_dates_o3 <- Bio_dates[grepl("O3", Bio_dates$Site), ]
Bio_dates_w1 <- Bio_dates[grepl("W1", Bio_dates$Site), ]

## sort the datasets by the date
Bio_dates_auba <- Bio_dates_auba[order(Bio_dates_auba$Date_new), ]
Bio_dates_bieb <- Bio_dates_bieb[order(Bio_dates_bieb$Date_new), ]
Bio_dates_o3 <- Bio_dates_o3[order(Bio_dates_o3$Date_new), ]
Bio_dates_w1 <- Bio_dates_w1[order(Bio_dates_w1$Date_new), ]

## Create function to average temp data from 12 months prior to biological sampling
## AUBA
# Iterate over each biological sampling
for (i in 1:nrow(Bio_dates_auba)) {
  bio_sampling_date <- as.Date(Bio_dates_auba$Date_new[i])
  start_date <- bio_sampling_date - days(1) - months(12) # Subtract 1 day and 1 year to start from the day before the sampling date
  
  # Filter temperature data for the 12-month period
  filtered_data <- Temp %>%
    filter(Date_new > start_date & Date_new <= bio_sampling_date) %>%
    slice(-n())  # Remove the last row from filtered_data
    
  # Calculate the number of missing days
  num_entries <- sum(!is.na(filtered_data$au_gap))
  missing_days <- 365 - num_entries
  
  # Calculate average temperature
  if (missing_days < 1) {
    average_temperature <- mean(filtered_data$au_gap, na.rm = T)
  } else {
    average_temperature <- NA
  }
  
  # Assign the average temperature to the sampling data
  Bio_dates_auba$Yryly_Temp[i] <- average_temperature
  Bio_dates_auba$Missing_Days[i] <- missing_days
}

## BIEB
# Iterate over each biological sampling
for (i in 1:nrow(Bio_dates_bieb)) {
  bio_sampling_date <- as.Date(Bio_dates_bieb$Date_new[i])
  start_date <- bio_sampling_date - days(1) - months(12) # Subtract 1 day and 1 year to start from the day before the sampling date
  
  # Filter temperature data for the 12-month period
  filtered_data <- Temp %>%
    filter(Date_new > start_date & Date_new <= bio_sampling_date) %>%
    slice(-n())  # Remove the last row from filtered_data
  
  # Calculate the number of missing days
  num_entries <- sum(!is.na(filtered_data$bieb_gap))
  missing_days <- 365 - num_entries
  
  # Calculate average temperature
  if (missing_days < 1) {
    average_temperature <- mean(filtered_data$bieb_gap, na.rm = T)
  } else {
    average_temperature <- NA
  }
  
  # Assign the average temperature to the sampling data
  Bio_dates_bieb$Yryly_Temp[i] <- average_temperature
  Bio_dates_bieb$Missing_Days[i] <- missing_days
}

## O3
# Iterate over each biological sampling
for (i in 1:nrow(Bio_dates_o3)) {
  bio_sampling_date <- as.Date(Bio_dates_o3$Date_new[i])
  start_date <- bio_sampling_date - days(1) - months(12) # Subtract 1 day and 1 year to start from the day before the sampling date
  
  # Filter temperature data for the 12-month period
  filtered_data <- Temp %>%
    filter(Date_new > start_date & Date_new <= bio_sampling_date) %>%
    slice(-n())  # Remove the last row from filtered_data
  
  # Calculate the number of missing days
  num_entries <- sum(!is.na(filtered_data$kio_gap))
  missing_days <- 365 - num_entries
  
  # Calculate average temperature
  if (missing_days < 1) {
    average_temperature <- mean(filtered_data$kio_gap, na.rm = T)
  } else {
    average_temperature <- NA
  }
  
  # Assign the average temperature to the sampling data
  Bio_dates_o3$Yryly_Temp[i] <- average_temperature
  Bio_dates_o3$Missing_Days[i] <- missing_days
}

## O3
# Iterate over each biological sampling
for (i in 1:nrow(Bio_dates_w1)) {
  bio_sampling_date <- as.Date(Bio_dates_w1$Date_new[i])
  start_date <- bio_sampling_date - days(1) - months(12) # Subtract 1 day and 1 year to start from the day before the sampling date
  
  # Filter temperature data for the 12-month period
  filtered_data <- Temp %>%
    filter(Date_new > start_date & Date_new <= bio_sampling_date) %>%
    slice(-n())  # Remove the last row from filtered_data
  
  # Calculate the number of missing days
  num_entries <- sum(!is.na(filtered_data$kiw_gap))
  missing_days <- 365 - num_entries
  
  # Calculate average temperature
  if (missing_days < 1) {
    average_temperature <- mean(filtered_data$kiw_gap, na.rm = T)
  } else {
    average_temperature <- NA
  }
  
  # Assign the average temperature to the sampling data
  Bio_dates_w1$Yryly_Temp[i] <- average_temperature
  Bio_dates_w1$Missing_Days[i] <- missing_days
}

## Combine everything into some big table for exporting
BioData_linkedto_TempData <- rbind(Bio_dates_auba, Bio_dates_bieb, Bio_dates_o3, Bio_dates_w1)
write.csv(BioData_linkedto_TempData, "output_data/BioData_linkedto_TempData.csv")

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
