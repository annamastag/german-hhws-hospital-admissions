# ==============================================================================
# Script: 06_07_01_compare_rfm3_and_rfm5_alert_definitions.R
#
# Purpose:
#   Evaluates the robustness of the study results to the alternative RFM3
#   definition of reconstructed heat-alert days. The script first pools
#   diagnosis-specific pre-implementation alert-day, Difference-in-Differences,
#   and post-implementation associations for RFM3. A second analysis block
#   compares the number and distribution of alert days classified by RFM3 and
#   RFM5.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Released cause-specific first-stage FDZ output table including
#     cardiovascular admissions
#   - Released reference-switch first-stage FDZ output table
#   - Environmental and heat-alert master dataset
#
# Outputs:
#   - CSV table containing pooled RFM3 pre-implementation, DiD, and
#     post-implementation estimates by diagnosis
#   - CSV table comparing RFM3 and RFM5 alert-day classifications
#   - Numerical summaries printed to the console
#
# Outcomes:
#   - All-cause admissions
#   - Cardiovascular admissions
#   - Infectious and parasitic disease admissions
#   - Respiratory admissions
#   - Urogenital admissions
#
# Analysis period:
#   May–September, 2000–2009
#
# Required packages:
#   mixmeta, dplyr, readxl, data.table, lubridate
#
# Repository note:
#   The script contains two sequential analysis blocks and clears the R
#   environment between them. The file name, path strings, and explanatory
#   comments were standardized. 
# ==============================================================================

library(mixmeta)
library(dplyr)
library(readxl)
library(data.table)
library(lubridate)

rm(list = ls())

# =====================================================================
# 1. DATEIPFADE DEFINIEREN
# =====================================================================

# A) HAUPT-ANALYSE
file_main_all <-
  "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
file_main_strat <-
  "data/processed/fdz_outputs/02_Analysis_FirstStage_default_allCause_and_Strat(mit_Card)_und_sens_Res_gh.xlsx"

# B) REFERENCE SWITCH (Beinhaltet AllAdmissions + Causes)
file_ref_combined <-
  "data/processed/fdz_outputs/Results_FirstStage_ReferenceSwitch_gh.xlsx"

# RFM Parameter
current_rfm <- "heat_alerts_rfm3"

# =====================================================================
# 2. HILFSFUNKTIONEN ZUM BEREINIGEN (Getrennt für Main & Ref-Switch)
# =====================================================================

clean_main_data <- function(df, rfm_filter) {
  df %>%
    filter(RFM == rfm_filter, Model == "model1") %>%
    mutate(
      Diagnosis = case_when(
        Diagnosis == "AllAdmissions" ~ "All-Cause",
        Diagnosis == "Card_0" ~ "Cardiovascular",
        Diagnosis == "Infect_par" ~ "Infectious & Parasitic",
        Diagnosis == "Resp_0" ~ "Respiratory",
        Diagnosis == "Uri_0" ~ "Urogenital",
        TRUE ~ Diagnosis
      ),
      RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
      CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
      CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
      RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
      CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
      CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD)))
    ) %>%
    filter(!is.na(RR_Alert))
}

clean_ref_data <- function(df, rfm_filter) {
  df %>%
    filter(RFM == rfm_filter, Model == "model1") %>%
    mutate(
      Diagnosis = case_when(
        Diagnosis == "AllAdmissions" ~ "All-Cause",
        Diagnosis == "Card_0" ~ "Cardiovascular",
        Diagnosis == "Infect_par" ~ "Infectious & Parasitic",
        Diagnosis == "Resp_0" ~ "Respiratory",
        Diagnosis == "Uri_0" ~ "Urogenital",
        TRUE ~ Diagnosis
      ),
      RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert_Post2005))),
      CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert_Post2005))),
      CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert_Post2005)))
    ) %>%
    filter(!is.na(RR_Alert))
}

# =====================================================================
# 3. DATEN EINLESEN & ZUSAMMENFÜHREN
# =====================================================================

print(paste("Lese Daten für", current_rfm, "ein..."))

df_main_raw_all <- read_excel(file_main_all)
df_main_raw_strat <- read_excel(file_main_strat)
df_main <- bind_rows(df_main_raw_all, df_main_raw_strat) %>% 
  clean_main_data(current_rfm)

df_ref_raw <- read_excel(file_ref_combined)
df_ref <- df_ref_raw %>% 
  clean_ref_data(current_rfm)

# =====================================================================
# 4. HILFSFUNKTION FÜR META-ANALYSE
# =====================================================================

get_pooled_values <- function(data, target) {
  if (nrow(data) < 2) return(list(est = NA, ci_l = NA, ci_u = NA))
  
  if (target %in% c("Alert", "Post")) {
    yi <- log(data$RR_Alert)
    sei <- (log(data$CI_Up_Alert) - log(data$CI_Low_Alert)) / (2 * 1.96)
  } else if (target == "DiD") {
    yi <- log(data$RR_DiD)
    sei <- (log(data$CI_Up_DiD) - log(data$CI_Low_DiD)) / (2 * 1.96)
  }
  
  vi <- sei^2
  valid <- is.finite(yi) & is.finite(vi) & vi > 0
  
  if(sum(valid) < 2) return(list(est = NA, ci_l = NA, ci_u = NA))
  
  m <- mixmeta(yi[valid] ~ 1, S = vi[valid], method = "reml")
  p <- predict(m, ci = TRUE)
  
  list(est = exp(as.numeric(p[1, "fit"])), 
       ci_l = exp(as.numeric(p[1, "ci.lb"])), 
       ci_u = exp(as.numeric(p[1, "ci.ub"])))
}

format_res <- function(res) {
  if (is.na(res$est)) return("N/A")
  sprintf("%.3f (%.3f - %.3f)", res$est, res$ci_l, res$ci_u)
}

# =====================================================================
# 5. EFFEKTE BERECHNEN UND TABELLE BAUEN
# =====================================================================

diagnoses <- c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital")

results_list <- list()

print("Berechne gepoolte Schätzer per Diagnose...")

for (diag in diagnoses) {
  
  sub_main <- df_main %>% filter(Diagnosis == diag)
  sub_ref  <- df_ref %>% filter(Diagnosis == diag)
  
  # 1. PRE-Effekt (aus Main)
  pre_res <- get_pooled_values(sub_main, "Alert")
  pre_text <- format_res(pre_res)
  
  # 2. DiD-Effekt (aus Main)
  did_res <- get_pooled_values(sub_main, "DiD")
  did_text <- format_res(did_res)
  
  # 3. POST-Effekt
  if (nrow(sub_ref) > 0) {
    post_res <- get_pooled_values(sub_ref, "Post")
    post_text <- format_res(post_res)
  } else {
    if (!is.na(pre_res$est) && !is.na(did_res$est)) {
      calc_post <- pre_res$est * did_res$est
      post_text <- sprintf("%.3f (Calculated, no CI)", calc_post)
    } else {
      post_text <- "N/A"
    }
  }
  
  results_list[[length(results_list) + 1]] <- data.frame(
    Diagnosis = diag,
    RR_Pre = pre_text,
    RR_DiD = did_text,
    RR_Post = post_text
  )
}

final_table <- bind_rows(results_list)

cat(paste0("\n============================================================================\n"))
cat(paste0("ERGEBNISSE SENSITIVITÄTSANALYSE: ", current_rfm, " (Model 1)\n"))
cat(paste0("============================================================================\n\n"))
print(final_table, row.names = FALSE)
cat(paste0("\n============================================================================\n"))

write.csv(
  final_table,
  paste0(
    "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
    "01_rfm3_vs_rfm5/Sensitivity_Results_Table_",
    current_rfm,
    ".csv"
  ),
  row.names = FALSE
)

#######################################################################

# =====================================================================
# 1. DATEN EINLESEN UND AUF SOMMER (2000-2009) FILTERN
# =====================================================================

print("Lese Masterfile ein...")

df <- fread(
  "data/external/fdz_input/Masterfile_Final_for_FDZ.csv"
)

# Datums-Formatierungen und Filter
df[, Datum := as.Date(Datum)]
df[, Jahr := year(Datum)]
df[, Monat := month(Datum)]

# Nur die Sommer-Monate der Analyse-Jahre behalten
df <- df[Jahr %in% 2000:2009 & Monat %in% 5:9]

# Neue Spalte für die Periode (Pre vs. Post)
df[, Periode := fifelse(Jahr < 2005, "1. Pre (2000-2004)", "2. Post (2005-2009)")]

# =====================================================================
# 2. BERECHNUNG DER WARNTAGE (RFM3 vs. RFM5)
# =====================================================================

print("Berechne Metriken für RFM3 und RFM5...")

# Schritt A: Aggregation pro Landkreis (AGS) und Periode
ags_summary <- df[, .(
  Tage_RFM5 = sum(heat_alerts_rfm5 == 1, na.rm = TRUE),
  Tage_RFM3 = sum(heat_alerts_rfm3 == 1, na.rm = TRUE)
), by = .(AGS, Periode)]

# Differenz berechnen (Wie viele Tage hat RFM3 *zusätzlich* gewarnt?)
ags_summary[, Extra_Tage_RFM3 := Tage_RFM3 - Tage_RFM5]

# Schritt B: Nationale Zusammenfassung (als saubere Tabelle aufbereitet)
national_summary <- data.frame(
  Periode = rep(c("1. Pre (2000-2004)", "2. Post (2005-2009)"), each = 4),
  Metric = rep(c(
    "1. Total Alert Days (Summed across all districts)",
    "2. Median Alert Days per District",
    "3. Number of Districts with >= 1 Alert",
    "4. Additional Alert Days (RFM3 vs RFM5)"
  ), 2),
  
  RFM5 = c(
    sum(ags_summary[Periode == "1. Pre (2000-2004)"]$Tage_RFM5),
    median(ags_summary[Periode == "1. Pre (2000-2004)"]$Tage_RFM5),
    sum(ags_summary[Periode == "1. Pre (2000-2004)"]$Tage_RFM5 > 0),
    0, # Referenzwert
    
    sum(ags_summary[Periode == "2. Post (2005-2009)"]$Tage_RFM5),
    median(ags_summary[Periode == "2. Post (2005-2009)"]$Tage_RFM5),
    sum(ags_summary[Periode == "2. Post (2005-2009)"]$Tage_RFM5 > 0),
    0  # Referenzwert
  ),
  
  RFM3 = c(
    sum(ags_summary[Periode == "1. Pre (2000-2004)"]$Tage_RFM3),
    median(ags_summary[Periode == "1. Pre (2000-2004)"]$Tage_RFM3),
    sum(ags_summary[Periode == "1. Pre (2000-2004)"]$Tage_RFM3 > 0),
    sum(ags_summary[Periode == "1. Pre (2000-2004)"]$Extra_Tage_RFM3),
    
    sum(ags_summary[Periode == "2. Post (2005-2009)"]$Tage_RFM3),
    median(ags_summary[Periode == "2. Post (2005-2009)"]$Tage_RFM3),
    sum(ags_summary[Periode == "2. Post (2005-2009)"]$Tage_RFM3 > 0),
    sum(ags_summary[Periode == "2. Post (2005-2009)"]$Extra_Tage_RFM3)
  )
)

# =====================================================================
# 3. OUTPUT
# =====================================================================

cat("\n=====================================================================================\n")
cat("DESKRIPTIVER VERGLEICH DER WARNTAGE: RFM5 vs. RFM3\n")
cat("=====================================================================================\n\n")
print(national_summary, row.names = FALSE)
cat("\n=====================================================================================\n")

# Als CSV speichern
write.csv(
  national_summary,
  paste0(
    "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
    "01_rfm3_vs_rfm5/Alert_Days_Comparison_RFM5_vs_RFM3.csv"
  ),
  row.names = FALSE
)
