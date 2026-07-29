# ==============================================================================
# Script: 05_02_create_ags_timeseries.R
#
# Purpose:
#   Processes annual hospital-admission data, harmonizes historical district
#   identifiers, aggregates admissions by date, sex, and age group, and merges
#   the resulting counts with district-level environmental and heat-alert data.
#   Complete daily district-specific time series are created by adding
#   zero-valued admission counts for combinations without recorded admissions.
#
# Inputs:
#   - Annual restricted hospital-admission data files
#   - Environmental and heat-alert data loaded by the paired configuration
#     script
#   - District crosswalk loaded by the paired configuration script
#
# Outputs:
#   - One CSV file named AGS_<district identifier>.csv for each harmonized
#     German administrative district
#
# Spatial unit:
#   Harmonized German administrative districts identified by AGS
#
# Temporal coverage:
#   2000–2022
#
# Demographic strata:
#   Sex and age groups defined in the paired configuration script
#
# Dependencies:
#   This script is sourced by
#   05_01_configure_ags_timeseries_creation.R. Required packages, file paths,
#   aggregation variables, demographic strata, and reference data are defined
#   in that configuration script.
#
# Data access:
#   Hospital data and district-specific time series are restricted and are not
#   included in the repository.
#
# Repository note:
#   The file name and explanatory comments were standardized.
# ==============================================================================

hosp_list <- list()

# 1. KRANKENHAUS-DATEN JAHRESWEISE VERARBEITEN

for (year in 2000:2022) {
  
  # DATEINAME LADEN
  file_name <- file.path(data_path, paste0(data_file_prefix, year, ".csv"))
  
  if (!file.exists(file_name)) {
    message("ACHTUNG: Datei existiert nicht: ", file_name)
    next
  }
  
  # EINLESEN
  df <- fread(file_name)
  vars_current <- intersect(vars_to_aggregate, names(df))
  if (length(vars_current) == 0) next
  
  # BEREINIGEN & FORMATIEREN 
  df <- df[!is.na(kreis_2022) & kreis_2022 != "" & kreis_2022 != "NA"]
  df <- df[age_group %in% target_ages & geschlecht %in% target_genders]
  
  df[, kreis_2022_num := suppressWarnings(as.integer(kreis_2022))]
  df <- df[!is.na(kreis_2022_num)] 
  df[, kreis_2022 := sprintf("%05d", kreis_2022_num)]
  df[, kreis_2022_num := NULL]
  
  # KREISHARMONISIERUNG ANWENDEN 
  if (exists("crosswalk") && !is.null(crosswalk)) {
    df[crosswalk, on = .(kreis_2022 = ags_alt), kreis_2022_neu := i.ags_neu]
    df[, kreis_2022 := fcoalesce(kreis_2022_neu, kreis_2022)]
    df[, kreis_2022_neu := NULL]
  } else {
    warning(paste("Jahr", year, "- Harmonisierung übersprungen: Crosswalk fehlt!"))
  }
  
  # MASTERFILE-ABGLEICH FÜR AGS AUF STAND 2022
  if (!exists("valid_ags_list") || is.null(valid_ags_list) || length(valid_ags_list) == 0) {
    stop("FEHLER: valid_ags_list fehlt. Pfad zum Masterfile im Config-Skript prüfen!")
  }
  
  # Zählen 
  ags_before <- uniqueN(df$kreis_2022)
  
  # Filter
  df <- df[kreis_2022 %in% valid_ags_list]
  
  ags_after <- uniqueN(df$kreis_2022)
  if (ags_before != ags_after) {
    print(paste("Jahr", year, ":", ags_before - ags_after, "falsche/veraltete AGS entfernt."))
  }
  
  ###################
  
  # Datum formatieren 
  df[, date := as.Date(datum)]
  
  # AGGREGIEREN (nach Datum, Geschlecht und Alter)
  agg_data <- df[, lapply(.SD, sum, na.rm = TRUE), 
                 by = .(date, AGS = kreis_2022, geschlecht, age_group), 
                 .SDcols = vars_current]
  
  hosp_list[[as.character(year)]] <- agg_data
}

# -------------------------------------------------------------------------
# JAHRE ZUSAMMENFÜHREN & NAs AUFFÜLLEN

full_hosp <- rbindlist(hosp_list, fill = TRUE)

for (j in vars_to_aggregate) {
  if (j %in% names(full_hosp)) {
    set(full_hosp, which(is.na(full_hosp[[j]])), j, 0)
  }
}

# Umwelt-Variablen definieren, die mit in die CSV sollen
env_vars_to_keep <- c(
  "T_mean_popw", "T_dew_popw", "holiday", 
  "heat_alerts_dwd", "heat_alerts_rfm1", "heat_alerts_rfm2", 
  "heat_alerts_rfm3", "heat_alerts_rfm4", "heat_alerts_rfm5"
)

# ------------------------------------------------------------------------
# CSVs PRO KREIS ERSTELLEN 

counter <- 0

for (current_ags in valid_ags_list) {
  
  counter <- counter + 1
  if(counter %% 20 == 0) cat(paste0("... ", counter, " AGS verarbeitet\n")) 
  
  # Umweltdaten für den speziellen Kreis holen
  cols_to_select <- c("date", intersect(env_vars_to_keep, names(env_data)))
  env_sub <- env_data[AGS == current_ags, ..cols_to_select]
  
  if(nrow(env_sub) == 0) next 
  
  # Grid bauen: Jeder Tag x Jedes Alter x Jedes Geschlecht
  grid <- CJ(date = env_sub$date, 
             geschlecht = target_genders, 
             age_group = target_ages)
  
  # Wetter/Feiertage ans Grid heften
  grid_env <- merge(grid, env_sub, by = "date", all.x = TRUE)
  
  # Krankenhaus-Daten für den Kreis filtern und anheften
  hosp_sub <- full_hosp[AGS == current_ags]
  hosp_sub[, AGS := NULL] 
  
  final_df <- merge(grid_env, hosp_sub, 
                    by = c("date", "geschlecht", "age_group"),
                    all.x = TRUE)
  
  # Nullen auffüllen für Tage ohne KH-Einweisungen
  for (col in vars_to_aggregate) {
    if (col %in% names(final_df)) {
      set(final_df, i = which(is.na(final_df[[col]])), j = col, value = 0)
    } else {
      # Falls Diagnose in diesem Kreis komplett fehlt
      final_df[, (col) := 0]
    }
  }
  
  # Sortieren und Exportieren
  setorder(final_df, date, geschlecht, age_group)
  
  file_out <- file.path(output_path, paste0("AGS_", current_ags, ".csv"))
  fwrite(final_df, file_out, sep = ";")
}
