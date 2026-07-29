# ==============================================================================
# Script: 01_03_prepare_daily_district_relative_humidity.R
# Original file: 00_humidity_relative_vorbereitung.R
#
# Purpose:
#   Aggregates hourly HOSTRADA relative-humidity raster data to daily metrics
#   and extracts district-level spatial summaries.
#
# Inputs:
#   - Hourly HOSTRADA NetCDF files containing relative humidity
#   - Harmonized district boundary shapefile containing AGS identifiers
#
# Outputs:
#   - One annual CSV file containing AGS, date, and district-level daily mean,
#     maximum, minimum, spatial 90th percentile, and spatial standard deviation
#     of relative humidity
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

ordner_pfad_nc <- "data/raw/hostrada/relative_humidity"
pfad_shp       <- "data/raw/spatial/vg2500/vg2500"
output_ordner  <- "data/intermediate/exposure/daily_relative_humidity"

if (!dir.exists(output_ordner)) dir.create(output_ordner)

# Load input data

vg_shp <- st_read(pfad_shp, quiet = TRUE)
kreise_terra <- vect(vg_shp)

dateien_liste <- list.files(path = ordner_pfad_nc, pattern = "\\.nc$", full.names = TRUE)
if(length(dateien_liste) == 0) stop("Keine .nc Dateien gefunden! Bitte Pfad prüfen.")

r_stunden <- rast(dateien_liste)

# Kreise anpassen falls nötig
if (crs(kreise_terra) != crs(r_stunden)) {
  print("Projeziere Kreise um...")
  kreise_terra <- project(kreise_terra, crs(r_stunden))
}

# Define the processing period

zeit_stempel_raw <- time(r_stunden)
alle_jahre <- unique(as.numeric(format(zeit_stempel_raw, "%Y")))

jahre_zu_berechnen <- alle_jahre[alle_jahre >= 2000 & alle_jahre <= 2022]

print(paste("Folgende Jahre werden bearbeitet:", paste(jahre_zu_berechnen, collapse=", ")))

##########
## LOOP ##
##########

for (jahr in jahre_zu_berechnen) {
  
  file_name <- file.path(output_ordner, paste0("RelHum_Kreise_Daily_", jahr, ".csv"))
  
  if (file.exists(file_name)) {
    print(paste("Datei für Jahr", jahr, "existiert schon. Überspringe..."))
    next
  }
  
  print(paste(">>> Starte Berechnung für Jahr:", jahr, "<<<"))
  
  idx_jahr_start <- which(as.numeric(format(zeit_stempel_raw, "%Y")) == jahr)[1]
  tage_im_jahr <- unique(as.Date(zeit_stempel_raw[as.numeric(format(zeit_stempel_raw, "%Y")) == jahr]))
  
  jahres_ergebnisse <- list()
  
  pb <- txtProgressBar(min = 0, max = length(tage_im_jahr), style = 3)
  
  # INNERER LOOP (TAG FÜR TAG)
  
  for (i in seq_along(tage_im_jahr)) {
    
    tag_datum <- tage_im_jahr[i]
    
    idx_tag <- which(as.Date(zeit_stempel_raw) == tag_datum)
    
    r_heute_24h <- r_stunden[[idx_tag]]
    
    r_mean <- mean(r_heute_24h, na.rm = TRUE) 
    r_max  <- max(r_heute_24h, na.rm = TRUE) 
    r_min  <- min(r_heute_24h, na.rm = TRUE)  
    
    # Auf Kreisebene bringen
    
    # Mean
    ex_mean <- terra::extract(r_mean, kreise_terra, fun = mean, na.rm = TRUE, ID = FALSE)
    
    # Max 
    ex_max  <- terra::extract(r_max, kreise_terra, fun = max, na.rm = TRUE, ID = FALSE)
    
    # Min 
    ex_min  <- terra::extract(r_min, kreise_terra, fun = min, na.rm = TRUE, ID = FALSE)
    
    # Q90 
    ex_q90  <- terra::extract(r_max, kreise_terra, 
                              fun = function(x) quantile(x, probs=0.9, na.rm=TRUE), 
                              ID = FALSE)
    
    # SD
    ex_sd   <- terra::extract(r_max, kreise_terra, fun = sd, na.rm = TRUE, ID = FALSE)
    
    # Tabelle bauen
    
    if("AGS" %in% names(kreise_terra)) { ags_col <- "AGS" } else { ags_col <- names(kreise_terra)[1] }
    ags_vec <- values(kreise_terra)[[ags_col]]
    
    df_tag <- data.frame(
      AGS = ags_vec,
      Datum = rep(tag_datum, length(ags_vec)),
      RelHum_Mean = as.vector(ex_mean[,1]),
      RelHum_Max  = as.vector(ex_max[,1]),
      RelHum_Min  = as.vector(ex_min[,1]),
      RelHum_Q90  = as.vector(ex_q90[,1]),
      RelHum_SD   = as.vector(ex_sd[,1])
    )
    
    jahres_ergebnisse[[i]] <- df_tag
    
    # RAM räumen 
    rm(r_heute_24h, r_mean, r_max, r_min, ex_mean, ex_max, ex_min, ex_q90, ex_sd, df_tag)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # SPEICHERN NACH JEDEM JAHR
  print(paste("Speichere Jahr", jahr, "..."))
  df_jahr_komplett <- do.call(rbind, jahres_ergebnisse)
  write.csv(df_jahr_komplett, file_name, row.names = FALSE)
  
  # RAM räumen
  rm(jahres_ergebnisse, df_jahr_komplett)
  gc()
}
