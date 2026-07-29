# ==============================================================================
# Script: 04_create_fdz_exposure_and_alert_master_dataset.R
#
# Purpose:
#   Creates the environmental and heat-alert master dataset used in the secure
#   FDZ analyses. The script combines district-level meteorological exposures,
#   public holiday indicators, official DWD heat-alert status, and reconstructed
#   pre-implementation alert classifications from five Random Forest model
#   variants.
#
# Inputs:
#   - Random Forest input dataset containing daily district-level
#     meteorological exposures and official heat-alert status
#   - District-level daily public holiday calendar
#   - Historical Random Forest prediction files for 2000–2004
#
# Outputs:
#   - Masterfile_Final_for_FDZ.rds
#   - Masterfile_Final_for_FDZ.csv
#
# Alert-period definition:
#   - 2000–2004: reconstructed alert classifications from the respective
#     Random Forest models
#   - From 2005 onward: official DWD heat-alert status
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Required packages:
#   tidyverse, lubridate, stringr
#
# Dependencies:
#   Requires outputs from:
#   - 02_02_prepare_random_forest_input_dataset.R
#   - 03_01_fit_over_and_undersampled_random_forest_models.R
#   - 03_02_fit_random_forest_without_population_weighted_predictors.R
#   - 03_03_fit_weighted_random_forest_and_evaluate_thresholds.R
#   - 04_00_01_create_district_holiday_calendar.R
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

library(tidyverse)
library(lubridate)
library(stringr)

# ---------------------------------------------------------
# 1. DATEN LADEN
# ---------------------------------------------------------

data_input <- read_rds(
  "data/processed/random_forest_input/RFM_Input_Data.rds"
) %>%
  mutate(AGS = as.character(AGS)) %>%
  mutate(AGS = str_pad(AGS, width = 5, pad = "0")) # 5 Stellen 

holiday_input <- read.csv(
  "data/intermediate/covariates/holidays/District_Holiday_Calendar_2000_2022.csv"
) %>%
  rename(Datum = date) %>%            
  mutate(Datum = as.Date(Datum)) %>%  
  mutate(AGS = as.character(AGS)) %>%     
  mutate(AGS = str_pad(AGS, width = 5, pad = "0")) %>% # 5 Stellen
  rename(holiday = Feiertag)

# Check nach Umwandlung
print(head(holiday_input))

# ---------------------------------------------------------
# 2. ZUSAMMENBAUEN 
# ---------------------------------------------------------

master_file <- data_input %>%
  select(Datum, AGS, temp_pop_w, dewpoint_pop_w, HW_Status) %>%
  
  left_join(holiday_input, by = c("Datum", "AGS")) %>%
  
 mutate(holiday = replace_na(holiday, 0)) %>%

  select(Datum, AGS, temp_pop_w, dewpoint_pop_w, holiday, HW_Status) %>%
  arrange(AGS, Datum)

# Check
print(head(master_file))

# ---------------------------------------------------------
# 3. AUSGEWÄHLTE RF MODELLE EINFÜGEN 
# ---------------------------------------------------------

model_files <- list(
  
  "RFM_1_Down_07" = paste0(
    "outputs/random_forest_models/03_01_sampling_comparison/",
    "RFM_Sim_Downsampling_Th0.7.csv"
  ),
  
  "RFM_2" = paste0(
    "outputs/random_forest_models/03_01_sampling_comparison/",
    "RFM_Sim_Oversampling_Th0.95.csv"
  ),
  
  "RFM_3_07" = paste0(
    "outputs/random_forest_models/03_03_weighted_threshold_comparison/",
    "RFM_Sim_2000-2004_Th0.7.csv"
  ),
  
  "RFM_4_no_pop_w" = paste0(
    "outputs/random_forest_models/",
    "03_02_without_population_weighted_predictors/",
    "RFM_Sim_NoPop_2000-2004.csv"
  ),
  
  "RFM_5_08" = paste0(
    "outputs/random_forest_models/03_03_weighted_threshold_comparison/",
    "RFM_Sim_2000-2004_Th0.8.csv"
  )
)

# -------------------------------------------------------------------------
# 4. FUNKTION ZUM EINFÜGEN DER SPALTEN
# -------------------------------------------------------------------------

add_simulation_column <- function(master_df, col_name, file_path) {
  
  cat("Verarbeite:", col_name, "aus Datei:", basename(file_path), "...\n")
  
  # 1. Simulations-Datei laden (enthält nur 2000-2004)
  sim_data <- read.csv(file_path) %>%
    mutate(Datum = as.Date(Datum)) %>%
    mutate(AGS = as.character(AGS)) %>%
    mutate(AGS = str_pad(AGS, width = 5, pad = "0")) %>%
    
    select(Datum, AGS, Pred) %>%
    # Umwandlung: yes -> 1, no -> 0
    mutate(Pred_Numeric = ifelse(Pred == "yes", 1, 0)) %>%
    select(Datum, AGS, Pred_Numeric)
  
    # 2. An das Masterfile hängen
    master_df <- master_df %>%
    left_join(sim_data, by = c("Datum", "AGS")) %>%
    
    # 3. Die Logik anwenden:
    # Wenn Jahr < 2005  -> Nimm Wert aus Simulation (Pred_Numeric)
    # Wenn Jahr >= 2005 -> Nimm Wert aus HW_Status (Original)
    mutate(
      !!sym(col_name) := case_when(
        year(Datum) < 2005 ~ Pred_Numeric, 
        TRUE ~ HW_Status                            
      )
    ) %>%
    
    # 4. Aufräumen (Hilfsspalte vom Join löschen)
    select(-Pred_Numeric)
  
  return(master_df)
}

# -------------------------------------------------------------------------
# 3. SCHLEIFE DURCH ALLE 5 MODELLE
# -------------------------------------------------------------------------

master_final <- master_file

for (name in names(model_files)) {
  path <- model_files[[name]]
  
  if (file.exists(path)) {
    master_final <- add_simulation_column(master_final, name, path)
  } else {
    cat("ACHTUNG: Datei nicht gefunden:", path, "\n")
  }
}

# -------------------------------------------------------------------------
# 4. VARIABLEN UMBENENNEN & SPEICHERN
# -------------------------------------------------------------------------

master_final <- master_final %>%
  rename(
    T_mean_popw      = temp_pop_w,
    T_dew_popw       = dewpoint_pop_w,
    heat_alerts_dwd  = HW_Status,
    heat_alerts_rfm1 = RFM_1_Down_07,
    heat_alerts_rfm2 = RFM_2,
    heat_alerts_rfm3 = RFM_3_07,
    heat_alerts_rfm4 = RFM_4_no_pop_w,
    heat_alerts_rfm5 = RFM_5_08
  )

# Check der neuen Namen
print(names(master_final))

# Speichern
write_rds(
  master_final,
  "data/external/fdz_input/Masterfile_Final_for_FDZ.rds"
)
write.csv(
  master_final,
  "data/external/fdz_input/Masterfile_Final_for_FDZ.csv",
  row.names = FALSE
)
