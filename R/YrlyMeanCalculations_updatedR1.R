library(tidyverse)

## Load data
Bio_dates <- read.csv("RawData/Temp_data/RMO_BioDataDates_updatedR1.csv", header = TRUE, sep = ",") |> 
    arrange(Site, Year, Month, Day)
CHELSA_Temp <- read.csv("RawData/Temp_data/CHELSA_temperature_data_20250804.csv", header = TRUE, sep = ",") |> 
  rename(Year = year,
         Month = month) |> 
  pivot_wider(
    names_from = site_code,
    values_from = temperature_C,
    id_cols = c(Year, Month)) |> 
  mutate(
    Day = 1,  # Set day to first day of month
    Date = paste(Month, Day, Year, sep = "/")  # Create Date column with M/D/Y format
  )  |> 
  arrange(Year, Month)

## Ensure date columns in both datasets are correct
Bio_dates$Date_new <- as.Date(Bio_dates$Date, format = "%m/%d/%Y")
CHELSA_Temp$Date_new <- as.Date(CHELSA_Temp$Date, format = "%m/%d/%Y")

## Function to calculate average temperature and missing days

## Function to calculate average temperature for monthly data (strict version)
calculate_temp_monthly <- function(Bio_dates_subset, site_name, temp_col) {
  for (i in 1:nrow(Bio_dates_subset)) {
    bio_sampling_date <- as.Date(Bio_dates_subset$Date_new[i])
    
    # Get exactly 12 months of data leading up to (and including) the sampling date
    start_year <- year(bio_sampling_date - months(11))
    start_month <- month(bio_sampling_date - months(11))
    end_year <- year(bio_sampling_date)
    end_month <- month(bio_sampling_date)
    
    filtered_data <- CHELSA_Temp |>
      filter(
        (Year == start_year & Month >= start_month) |
        (Year > start_year & Year < end_year) |
        (Year == end_year & Month <= end_month)
      )
    
    num_entries <- sum(!is.na(filtered_data[[temp_col]]))
    average_temperature <- mean(filtered_data[[temp_col]], na.rm = TRUE)
    
    Bio_dates_subset$Yryly_Temp[i] <- average_temperature
    Bio_dates_subset$Num_Months[i] <- num_entries
  }
  return(Bio_dates_subset)
}

## Site names and corresponding temperature columns
sites <- c("Auba", "Bieb", "O3", "W1")
temp_columns <- c("Auba", "Bieb", "KiO3", "KiW1")

## Calculate temperature and number of months for each site
BioData_linkedto_TempData <- lapply(seq_along(sites), function(i) {
  site_name <- sites[i]
  temp_col <- temp_columns[i]
  Bio_dates_subset <- Bio_dates[grepl(site_name, Bio_dates$Site), ]
  Bio_dates_subset <- Bio_dates_subset[order(Bio_dates_subset$Date_new), ]
  calculate_temp_monthly(Bio_dates_subset, site_name, temp_col)
})

## Combine everything into a big table for exporting
BioData_linkedto_TempData2 <- do.call(rbind, BioData_linkedto_TempData)

## Add Sample column as requested (Site_Year_Month_Day)
BioData_linkedto_TempData2 <- BioData_linkedto_TempData2 |>
  mutate(Sample = paste(Site, Year, Month, Day, sep = "_")) |>
  select(Sample, everything()) 

write.csv(BioData_linkedto_TempData2, "RawData/Temp_data/BioData_linkedto_TempData_updatedR1.csv", row.names = FALSE)

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
