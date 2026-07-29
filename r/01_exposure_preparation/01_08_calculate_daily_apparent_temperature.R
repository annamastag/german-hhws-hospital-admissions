# ==============================================================================
# Script: 01_08_calculate_daily_apparent_temperature.R
# Original file: 00_calculate_AT.R
#
# Purpose:
#   Calculates district-level daily apparent temperature by combining prepared
#   air-temperature, dew-point temperature, and wind-speed data.
#
# Inputs:
#   - Annual district-level daily air-temperature CSV files
#   - Annual district-level daily dew-point temperature CSV files
#   - Annual district-level daily wind-speed CSV files
#
# Outputs:
#   - One annual CSV file containing AGS, date, daily mean apparent temperature,
#     and daily maximum apparent temperature
#
# Method:
#   Apparent temperature is calculated using the Steadman equation implemented
#   in the original analysis.
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   dplyr, readr
#
# Repository note:
#   Paths are project-relative and assume that the working directory is set to
#   the repository root. Only path strings and explanatory comments were
#   standardized for repository publication.
# ==============================================================================

# Load required packages

library(dplyr)
library(readr)

####################
## SETUP & CONFIG ##
####################

path_temp <- "data/intermediate/exposure/daily_air_temperature"
path_dew  <- "data/intermediate/exposure/daily_dew_point_temperature"
path_wind <- "data/intermediate/exposure/daily_wind_speed"
path_out  <- "data/intermediate/exposure/daily_apparent_temperature"

if(!dir.exists(path_out)) dir.create(path_out)

years <- 2000:2022

# FUNKTION: Australian Apparent Temp
calc_at <- function(temp, dew, wind) {
  # Formel nach Steadman 1994
  e <- 6.112 * exp((17.67 * dew) / (dew + 243.5))
  at <- temp + (0.33 * e) - (0.70 * wind) - 4.00
  return(at)
}

##########
## LOOP ##
##########

for (yr in years) {
  print(paste("Erstelle Master-Tabelle für Jahr:", yr))
  
  f_temp <- file.path(path_temp, paste0("Temperature_Kreise_Daily_", yr, ".csv"))
  f_dew  <- file.path(path_dew,  paste0("Dewpoint_Kreise_Daily_", yr, ".csv"))
  f_wind <- file.path(path_wind, paste0("Wind_Kreise_Daily_", yr, ".csv"))
  
  if (!all(file.exists(f_temp, f_dew, f_wind))) {
    warning(paste("Dateien fehlen für Jahr", yr))
    next
  }
  
  # Einlesen
  d_temp <- read_csv(f_temp, show_col_types = FALSE)
  d_dew  <- read_csv(f_dew, show_col_types = FALSE)
  d_wind <- read_csv(f_wind, show_col_types = FALSE)
  
  # Joinen
  df_full <- d_temp %>%
    inner_join(d_dew, by = c("AGS", "Datum")) %>%
    inner_join(d_wind, by = c("AGS", "Datum"))
  
  # Apparent Temperature berechnen
  df_final <- df_full %>%
    mutate(
      
      # Durchschnittliche AT
      AT_Mean = calc_at(temp = Temp_Mean, 
                        dew  = Dew_Mean, 
                        wind = Wind_Mean),
      # Maximale AT 
      AT_Max = calc_at(temp = Temp_Max, 
                       dew  = Dew_Mean, 
                       wind = Wind_Mean)
    )%>%
    
  select(AGS, Datum, AT_Mean, AT_Max)
  
  # Speichern
  write_csv(df_final, file.path(path_out, paste0("AT_Kreise_Daily_", yr, ".csv")))
  
  # Speicher freigeben
  rm(d_temp, d_dew, d_wind, df_full, df_final)
  gc()
  
}
