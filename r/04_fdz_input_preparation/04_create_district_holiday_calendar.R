# ==============================================================================
# Script: 04_create_district_holiday_calendar.R
#
# Purpose:
#   Creates a complete district-level daily calendar containing a binary public
#   holiday indicator. State-specific public holidays are assigned to German
#   administrative districts using the federal-state identifier contained in
#   the harmonized district boundary data.
#
# Inputs:
#   - CSV table of state-specific public holidays from 2000 through 2022
#   - Harmonized German district boundary shapefile containing AGS and
#     federal-state identifiers
#
# Outputs:
#   - CSV table containing one row per district and date with a binary public
#     holiday indicator
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   data.table, sf
#
# Dependencies:
#   The output is used by
#   04create_fdz_exposure_and_alert_master_dataset.R.
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

# ---------------------------------------------------------
# PAKETE LADEN
# ---------------------------------------------------------

library(data.table)
library(sf)

# ---------------------------------------------------------
# FEIERTAGE LADEN UND VORBEREITEN
# ---------------------------------------------------------

feiertage_2000_2022 <- fread(
  "data/raw/calendar/public_holidays/feiertage_2000_2022.csv"
)
feiertage_2000_2022[, Datum := as.Date(Datum)]

feiertage_flag <- feiertage_2000_2022[, .(Datum, Bundesland, Feiertag = 1)]

# ---------------------------------------------------------
# SHAPEFILE LADEN (AGS & ZUORDNUNG)
# ---------------------------------------------------------

pfad_shp <- "data/raw/spatial/vg2500/vg2500.shp"
shp_data <- st_read(pfad_shp, quiet = TRUE)

# AGS und LKZ extrahieren
ags_zuordnung <- as.data.table(shp_data)[, .(AGS, LKZ)]
ags_zuordnung <- unique(ags_zuordnung)

# ---------------------------------------------------------
# TIMELINE ALLE TAGE X ALLE AGS
# ---------------------------------------------------------

# Alle Tage von 2000 bis 2022
alle_tage <- seq(as.Date("2000-01-01"), as.Date("2022-12-31"), by = "day")

data <- CJ(AGS = ags_zuordnung$AGS, date = alle_tage)

# ---------------------------------------------------------
# MERGE 1: BL AGS ZUORDNEN
# ---------------------------------------------------------

data <- merge(data, ags_zuordnung, by = "AGS", all.x = TRUE)

# ---------------------------------------------------------
# MERGE 2: FEIERTAG JA ODER NEIN
# ---------------------------------------------------------

data <- merge(
  data,
  feiertage_flag,
  by.x = c("date", "LKZ"),
  by.y = c("Datum", "Bundesland"),
  all.x = TRUE
)

# ---------------------------------------------------------
# AUFRÄUMEN
# ---------------------------------------------------------

# Spalte LKZ umbenennen 
setnames(data, "LKZ", "Bundesland")

# NAs werden zu 0
data[is.na(Feiertag), Feiertag := 0]

# ---------------------------------------------------------
# CHECK
# ---------------------------------------------------------

print(paste("Gefundene Feiertage:", sum(data$Feiertag)))
head(data)

# ---------------------------------------------------------
# SPEICHERN
# ---------------------------------------------------------

write.csv(
  data,
  "data/intermediate/covariates/holidays/District_Holiday_Calendar_2000_2022.csv",
  row.names = FALSE
)
