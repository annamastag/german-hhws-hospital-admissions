# ==============================================================================
# Script: 01_07_calculate_population_weighted_exposures.R
# Original file: 00_pop_weighted_mean.R
#
# Purpose:
#   Calculates daily population-weighted district-level exposure estimates for
#   air temperature, dew-point temperature, and wind speed.
#
# Inputs:
#   - Hourly HOSTRADA NetCDF files for air temperature, dew-point temperature,
#     and wind speed
#   - Harmonized district boundary shapefile containing AGS identifiers
#   - Population raster derived from the 2011 German Census
#
# Outputs:
#   - One CSV file per meteorological variable containing AGS, date, and the
#     population-weighted daily mean exposure estimate
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   terra, sf, exactextractr, dplyr, lubridate, tidyr
#
# Repository note:
#   Paths are project-relative and assume that the working directory is set to
#   the repository root. Only path strings and explanatory comments were
#   standardized for repository publication.
# ==============================================================================

# Load required packages

library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(lubridate)
library(tidyr)

####################
## SETUP & CONFIG ##
####################

path_shapes <- "data/raw/spatial/vg2500/vg2500"
path_pop    <- "data/raw/population/census_2011/population_1km.tif"

folder_list <- c(
  "temp"     = "data/raw/hostrada/air_temperature",
  "dewpoint" = "data/raw/hostrada/dew_point_temperature",
  "wind"     = "data/raw/hostrada/wind_speed"
)

target_years <- 2000:2022

# SETUP shapefiles
kreise <- st_read(path_shapes, quiet = TRUE)
kreise$AGS <- sprintf("%05d", as.numeric(kreise$AGS))
pop    <- rast(path_pop)

# Check
if (st_crs(kreise)$wkt != crs(pop)) {
  kreise <- st_transform(kreise, crs(pop))
}

# Funktion
process_daily_pop_weighted <- function(folder_path, var_name, pop_raster, shapes) {
  
  nc_files <- list.files(folder_path, pattern = "\\.nc$", full.names = TRUE)
  
  output_csv <- paste0(
    "data/intermediate/exposure/population_weighted/",
    var_name,
    "_DAILY_pop_weighted.csv"
  )

  header_df <- data.frame(AGS = character(), Datum = character(), Value = numeric())
  colnames(header_df) <- c("AGS", "Datum", paste0(var_name, "_pop_w"))
  write.csv(header_df, output_csv, row.names = FALSE)
  
  print(paste(">>> Starte TÄGLICHE Berechnung für:", var_name))
  
  for (f_path in nc_files) {
    fname <- basename(f_path)
    
    tryCatch({
      nc <- rast(f_path)
      times <- time(nc)
      if (any(is.na(times))) { 
      }
      
      file_years <- unique(year(times))
      if (length(intersect(file_years, target_years)) == 0) { next }
      
      print(paste("   Verarbeite Datei:", fname))
      
      r_daily <- tapp(nc, index = "days", fun = mean, na.rm = TRUE)
      daily_dates <- unique(as.Date(times))
      
      crs(r_daily) <- "epsg:3034"
    
      if (crs(r_daily) != crs(pop_raster)) {
        r_daily <- project(r_daily, pop_raster)
      }
      
      if (!compareGeom(r_daily, pop_raster, stopOnError = FALSE)) {
        r_daily <- resample(r_daily, pop_raster, method = "bilinear")
      }
      
      # EXTRAKTION
      extract_df <- exact_extract(r_daily, shapes, weights = pop_raster, fun = "weighted_mean", progress = FALSE)
      
      # FORMATIERUNG
      extract_df$AGS <- shapes$AGS
      colnames(extract_df)[1:length(daily_dates)] <- as.character(daily_dates)
 
      df_long <- extract_df %>%
        pivot_longer(
          cols = all_of(as.character(daily_dates)), 
          names_to = "Datum", 
          values_to = paste0(var_name, "_pop_w")
        )
      
      # SPEICHERN
      write.table(df_long, output_csv, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      
      rm(nc, r_daily, extract_df, df_long)
      gc()
      
    }, error = function(e) {
      print(paste("FEHLER bei Datei", fname, ":", e$message))
    })
  }
}

# START
for (v_name in names(folder_list)) {
  try({
    process_daily_pop_weighted(folder_list[[v_name]], v_name, pop, kreise)
  })
  gc()
}
