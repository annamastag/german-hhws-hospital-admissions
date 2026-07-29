# ==============================================================================
# Script: 01_05_prepare_daily_district_downwelling_radiation.R
# Original file: 00_radiation_downwelling_vorbereitung.R
#
# Purpose:
#   Aggregates hourly HOSTRADA downwelling-radiation raster data and derives
#   daily district-level radiation metrics and sunshine hours.
#
# Inputs:
#   - Hourly HOSTRADA NetCDF files containing downwelling radiation
#   - Harmonized district boundary shapefile containing AGS identifiers
#
# Outputs:
#   - One annual CSV file containing AGS, date, district-level mean radiation,
#     maximum radiation, and mean sunshine hours based on the applied
#     120 W/m² threshold
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   terra, sf, tidyverse
#
# Repository note:
#   Paths are project-relative and assume that the working directory is set to
#   the repository root. Only path strings and explanatory comments were
#   standardized for repository publication.
# ==============================================================================

# Load required packages

library(terra)
library(sf)
library(tidyverse) 

####################
## SETUP & CONFIG ##
####################

ordner_pfad_nc <- "data/raw/hostrada/downwelling_radiation"
pfad_shp       <- "data/raw/spatial/vg2500/vg2500"
output_ordner  <- "data/intermediate/exposure/daily_downwelling_radiation"

if (!dir.exists(output_ordner)) dir.create(output_ordner)

# Load input data

vg_shp <- st_read(pfad_shp, quiet = TRUE)
kreise_terra <- vect(vg_shp)

dateien_liste <- list.files(path = ordner_pfad_nc, pattern = "\\.nc$", full.names = TRUE)
if(length(dateien_liste) == 0) stop("Keine NetCDF Dateien gefunden!") 

r_stunden <- rast(dateien_liste)

if (crs(kreise_terra) != crs(r_stunden)) {
  print("Projeziere Shapefile auf Raster-CRS...")
  kreise_terra <- project(kreise_terra, crs(r_stunden))
}

# ZEITRAUM

zeit_stempel_raw <- time(r_stunden)

# Datum-Vektoren erstellen
alle_tage_vektor <- as.Date(zeit_stempel_raw) 
alle_jahre_vektor <- as.numeric(format(alle_tage_vektor, "%Y"))
alle_jahre_unique <- unique(alle_jahre_vektor)

# Zeitraum festlegen (2000 - 2022)
jahre_zu_berechnen <- alle_jahre_unique[alle_jahre_unique >= 2000 & alle_jahre_unique <= 2022] 

print(paste("Bearbeite Jahre:", paste(jahre_zu_berechnen, collapse=", ")))

##########
## LOOP ##
##########

for (jahr in jahre_zu_berechnen) {
  
  file_name <- file.path(output_ordner, paste0("Radiation_Kreise_Daily_", jahr, ".csv"))
  
  if (file.exists(file_name)) {
    print(paste("Jahr", jahr, "schon fertig. Next..."))
    next
  }
  
  print(paste(">>> Starte Analyse für Jahr:", jahr, "<<<"))
  
  tage_im_jahr <- unique(alle_tage_vektor[alle_jahre_vektor == jahr])
  
  jahres_ergebnisse <- list()
  pb <- txtProgressBar(min = 0, max = length(tage_im_jahr), style = 3)
  
  # LOOP ÜBER TAGE
  for (i in seq_along(tage_im_jahr)) {
    
    tag_datum <- tage_im_jahr[i]
    
    # Indizes für die Stunden des aktuellen Tages finden (meist 24)
    idx_tag <- which(alle_tage_vektor == tag_datum)
    
    # Raster-Slice für diesen Tag laden (Stack mit z.B. 24 Layern)
    r_heute <- r_stunden[[idx_tag]]
    
    # --- BERECHNUNG DER TAGES-METRIKEN (Pixel-Basis) ---
    
    # 1. Mean Radiation (Durchschnittliche Intensität)
    r_pixel_mean <- mean(r_heute, na.rm = TRUE)
    
    # 2. Max Radiation (Wann war die Belastungsspitze?)
    # WICHTIG: Erst Max pro Pixel berechnen (zeitlich), DANN räumlich mitteln
    r_pixel_max <- max(r_heute, na.rm = TRUE)
    
    # 3. Sunshine Hours (Summe der Stunden > 120 W/m²)
    # Erzeugt 0/1 Raster und summiert es auf
    r_pixel_sun <- sum(r_heute > 120, na.rm = TRUE)
    
    # --- RÄUMLICHE EXTRAKTION (Zonal Statistics auf Kreisebene) ---
    # Es wird 'fun = mean' genutzt, um den Durchschnittswert für den ganzen Landkreis zu bekommen
    
    # Extraktion (ID = FALSE liefert Dataframe)
    ex_mean <- terra::extract(r_pixel_mean, kreise_terra, fun = mean, na.rm = TRUE, ID = FALSE)
    ex_max  <- terra::extract(r_pixel_max,  kreise_terra, fun = mean, na.rm = TRUE, ID = FALSE)
    ex_sun  <- terra::extract(r_pixel_sun,  kreise_terra, fun = mean, na.rm = TRUE, ID = FALSE)
    
    # --- TABELLE BAUEN ---
    
    # AGS sicherstellen
    if("AGS" %in% names(kreise_terra)) { 
      ags_vec <- values(kreise_terra)[["AGS"]] 
    } else { 
      ags_vec <- values(kreise_terra)[,1] # Fallback: Erste Spalte
    }
    
    df_tag <- data.frame(
      AGS = ags_vec,
      Datum = rep(tag_datum, length(ags_vec)),
      Radiation_Mean = round(as.vector(ex_mean[,1]), 2), 
      Radiation_Max  = round(as.vector(ex_max[,1]), 2),
      Sunshine_Hours = round(as.vector(ex_sun[,1]), 2)
    )
    
    jahres_ergebnisse[[i]] <- df_tag
    
    # Speicher freiräumen (Wichtig in Loops!)
    rm(r_heute, r_pixel_mean, r_pixel_max, r_pixel_sun, ex_mean, ex_max, ex_sun, df_tag)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # --- SPEICHERN ---
  print(paste("Speichere CSV für Jahr", jahr, "..."))
  df_jahr_komplett <- do.call(rbind, jahres_ergebnisse)
  write.csv(df_jahr_komplett, file_name, row.names = FALSE)
  
  # Aufräumen für nächstes Jahr
  rm(jahres_ergebnisse, df_jahr_komplett)
  gc() 
}
