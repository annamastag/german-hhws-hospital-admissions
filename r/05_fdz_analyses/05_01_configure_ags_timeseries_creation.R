# ==============================================================================
# Script: 05_01_configure_ags_timeseries_creation.R
#
# Purpose:
#   Defines paths, aggregation variables, demographic strata, logging settings,
#   and reference data before creating district-specific daily hospital
#   admission time series.
#
# Paired analysis script:
#   05_02_create_ags_timeseries.R
#
# Inputs:
#   - Annual restricted hospital-admission data files
#   - Environmental and heat-alert master dataset
#   - Crosswalk for harmonizing historical district identifiers
#
# Outputs:
#   - One district-specific CSV time series per harmonized AGS
#   - Execution log
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   data.table, dplyr, stringr
#
# Data access:
#   Hospital data are restricted and are not included in the repository.
#
# Repository note:
#   File names, path strings, and explanatory comments were standardized.
# ==============================================================================

# =====================================================================
# 1. LADEN VON PACKAGES
# =====================================================================

library(data.table)
library(dplyr)    
library(stringr)

# =====================================================================
# 2. SPEICHER LEEREN
# =====================================================================

rm(list=ls())

# =====================================================================
# 3. DEFINITION ARBEITSUMGEBUNG
# =====================================================================

# Define the computing environment

FDZ <- 1

if (FDZ == 0) {
  
  data_path <- "data/restricted/local/annual_hospital_data"
  
  output_path <- "data/restricted/local/ags_timeseries"
  
  syntax_path <- "R/05_FDZ_Analyses"
  
  master_path <- paste0(
    "data/external/fdz_input/",
    "Masterfile_Final_for_FDZ.csv"
  )
  
  crosswalk_path <- paste0(
    "data/external/fdz_input/",
    "KDFV_Plausi_Kreis-AGSaufbereitetR.csv"
  )
}

if (FDZ == 1) {
  
  data_path <- "data/restricted/fdz/annual_hospital_data"
  
  output_path <- "data/restricted/fdz/ags_timeseries"
  
  syntax_path <- "R/04_FDZ_Analyses"
  
  master_path <- paste0(
    "data/external/fdz_input/",
    "Masterfile_Final_for_FDZ.csv"
  )
  
  crosswalk_path <- paste0(
    "data/external/fdz_input/",
    "KDFV_Plausi_Kreis-AGSaufbereitetR.csv"
  )
}

# =====================================================================
# 4. MAKROS & STRATA
# =====================================================================

data_file_prefix <- "output_taegl_einweisungen_"
syntax_name <- "05_02_create_ags_timeseries"
log_file_name <- "Log_Create_AGS_Timeseries.txt"

target_genders <- c("maennlich", "weiblich")
target_ages    <- c("0-64", "65+") 

vars_to_aggregate <- c(
  "Infect_par", "Resp_0", "Resp_1", "Resp_2",
  "Card_0", "Card_1", "Card_2", "Ment_0", "Ment_1",
  "Uri_0",  "Uri_1",  "Uri_2", "Other"
)

# START DER AUFZEICHNUNG
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
sink(file.path(output_path, log_file_name), type = c("output", "message"), split = TRUE)

# =======================================================================
# 5. REFERENZDATEN LADEN
# =======================================================================

if(file.exists(crosswalk_path)) {
  crosswalk <- fread(crosswalk_path)
  setnames(crosswalk, old = c("ehemaliger Kreis", "neuer Kreis"), new = c("ags_alt", "ags_neu"), skip_absent = TRUE)
  crosswalk[, ags_alt := sprintf("%05d", as.integer(ags_alt))]
  crosswalk[, ags_neu := sprintf("%05d", as.integer(ags_neu))]
  print("Config: Crosswalk erfolgreich geladen.")
} else {
  crosswalk <- NULL
  print("Config WARNUNG: Crosswalk nicht gefunden.")
}

if(file.exists(master_path)) {
  print("Lade Umweltdaten...")
  
  env_raw <- read.csv(master_path, sep = ",", stringsAsFactors = FALSE)
  cleaned_names <- gsub('"', '', names(env_raw))
  names(env_raw) <- cleaned_names
  
  if (!"AGS" %in% cleaned_names || !"Datum" %in% cleaned_names) {
    stop("FEHLER: Spalte 'AGS' oder 'Datum' in den Umweltdaten nicht gefunden!")
  }
  
  setDT(env_raw)
  env_raw[, AGS := gsub('"', '', as.character(AGS))]
  env_raw[, AGS_num := suppressWarnings(as.integer(AGS))]
  
  env_data <- env_raw[!is.na(AGS_num)] 
  env_data[, AGS := sprintf("%05d", AGS_num)]
  env_data[, AGS_num := NULL]
  
  env_data[, date := as.Date(Datum, format = "%Y-%m-%d")]
  
  valid_ags_list <- unique(env_data$AGS)
  print(paste("Config: Umweltdaten für", length(valid_ags_list), "AGS geladen."))
  rm(env_raw) 
} else {
  valid_ags_list <- NULL
  env_data <- NULL
  stop(paste("FEHLER: Umweltdaten nicht gefunden unter:", master_path))
}

# =======================================================================
# 6. AUSFÜHRUNG ANALYSE-SKRIPT
# =======================================================================

print("Starte Analyse-Skript...")
source(file.path(syntax_path, paste0(syntax_name, ".R")), echo = TRUE, max.deparse.length = 99999)

print("Skript beendet. Schließe Log.")

sink()
