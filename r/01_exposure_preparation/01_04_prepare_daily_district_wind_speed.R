# ==============================================================================
# Script: 01_04_prepare_daily_district_wind_speed.R
# Original file: 00_windspeed_vorbereitung.R
#
# Purpose:
#   Aggregates hourly HOSTRADA wind-speed raster data to daily mean values and
#   extracts district-level spatial summaries.
#
# Inputs:
#   - Hourly HOSTRADA NetCDF files containing wind speed
#   - Harmonized district boundary shapefile containing AGS identifiers
#
# Outputs:
#   - One annual CSV file containing AGS, date, and district-level spatial mean,
#     minimum, 10th percentile, and standard deviation of daily mean wind speed
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   terra, sf
#
# Repository note:
#   Paths are project-relative and assume that the working directory is set to
#   the repository root. Only path strings and explanatory comments were
#   standardized for repository publication. 
# ==============================================================================

# Load required packages

library(terra)
library(sf)

####################
## SETUP & CONFIG ##
####################

ordner_pfad_nc <- "data/raw/hostrada/wind_speed"
pfad_shp       <- "data/raw/spatial/vg2500/vg2500"
output_ordner  <- "data/intermediate/exposure/daily_wind_speed"

if (!dir.exists(output_ordner)) dir.create(output_ordner)

# Load input data

print("Lade Shapefile...")
vg_shp <- st_read(pfad_shp, quiet = TRUE)
kreise_terra <- vect(vg_shp)

print("Lade NetCDF Indizes...")
dateien_liste <- list.files(path = ordner_pfad_nc, pattern = "\\.nc$", full.names = TRUE)
if(length(dateien_liste) == 0) stop("Keine NetCDF Dateien gefunden!") 

r_stunden <- rast(dateien_liste)

# Kreise anpassen falls nötig
if (crs(kreise_terra) != crs(r_stunden)) {
  kreise_terra <- project(kreise_terra, crs(r_stunden))
}

# ZEITRAUM
print("Verarbeite Zeitstempel...")
zeit_stempel_raw <- time(r_stunden)

# Datum EINMAL berechnen
alle_tage_vektor <- as.Date(zeit_stempel_raw) 
alle_jahre_vektor <- as.numeric(format(alle_tage_vektor, "%Y"))

alle_jahre_unique <- unique(alle_jahre_vektor)

# Training Zeitraum festlegen
jahre_zu_berechnen <- alle_jahre_unique[alle_jahre_unique >= 2000 & alle_jahre_unique <= 2022] 

print(paste("Bearbeite Jahre:", paste(jahre_zu_berechnen, collapse=", ")))

##########
## LOOP ##
##########

for (jahr in jahre_zu_berechnen) {
  
  file_name <- file.path(output_ordner, paste0("Wind_Kreise_Daily_", jahr, ".csv"))
  
  if (file.exists(file_name)) {
    print(paste("Jahr", jahr, "schon fertig. Next..."))
    next
  }
  
  print(paste(">>> Starte Analyse für Jahr:", jahr, "<<<"))
  
  # Tage filtern
  tage_im_jahr <- unique(alle_tage_vektor[alle_jahre_vektor == jahr])
  
  jahres_ergebnisse <- list()
  pb <- txtProgressBar(min = 0, max = length(tage_im_jahr), style = 3)
  
  for (i in seq_along(tage_im_jahr)) {
    
    tag_datum <- tage_im_jahr[i]
    
    # Index finden
    idx_tag <- which(alle_tage_vektor == tag_datum)
    
    r_heute <- r_stunden[[idx_tag]]
    
    # Zeitliche Aggregation (Stunde zu Tag) 
    r_daily_mean <- mean(r_heute, na.rm = TRUE)
    
    # Räumliche Extraktion (Raster zu Kreise) 
    
    # 1. Mean
    ex_mean <- terra::extract(r_daily_mean, kreise_terra, fun = mean, na.rm = TRUE, ID = FALSE)
    
    # 2. Min (Wichtig für Stagnation)
    ex_min  <- terra::extract(r_daily_mean, kreise_terra, fun = min, na.rm = TRUE, ID = FALSE)
    
    # 3. Q10 (Wichtig für Stagnation in Tälern)
    ex_q10  <- terra::extract(r_daily_mean, kreise_terra, 
                              fun = function(x) quantile(x, probs=0.1, na.rm=TRUE), 
                              ID = FALSE)
    
    # 4. SD
    ex_sd   <- terra::extract(r_daily_mean, kreise_terra, fun = sd, na.rm = TRUE, ID = FALSE)
    
    
    # Tabelle bauen
    if("AGS" %in% names(kreise_terra)) { ags_col <- "AGS" } else { ags_col <- names(kreise_terra)[1] }
    ags_vec <- values(kreise_terra)[[ags_col]]
    
    df_tag <- data.frame(
      AGS = ags_vec,
      Datum = rep(tag_datum, length(ags_vec)),
      Wind_Mean = as.vector(ex_mean[,1]), 
      Wind_Min  = as.vector(ex_min[,1]),  
      Wind_Q10  = as.vector(ex_q10[,1]), 
      Wind_SD   = as.vector(ex_sd[,1])
    )
    
    jahres_ergebnisse[[i]] <- df_tag
    
    # Speicher freiräumen
    rm(r_heute, r_daily_mean, ex_mean, ex_min, ex_q10, ex_sd, df_tag)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Speichern
  df_jahr_komplett <- do.call(rbind, jahres_ergebnisse)
  write.csv(df_jahr_komplett, file_name, row.names = FALSE)
  
  rm(jahres_ergebnisse, df_jahr_komplett)
  gc()
}
