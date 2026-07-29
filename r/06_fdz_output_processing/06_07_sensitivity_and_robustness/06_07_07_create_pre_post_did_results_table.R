# ==============================================================================
# Script: 06_07_07_create_pre_post_did_results_table.R
#
# Purpose:
#   Creates a pooled results table containing pre-implementation alert-day,
#   post-implementation alert-day, and Difference-in-Differences estimates
#   across diagnosis groups and model specifications.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Released cause-specific first-stage FDZ output table including
#     cardiovascular admissions
#   - Released reference-switch first-stage FDZ output table
#
# Outputs:
#   - CSV table containing pooled RRs, 95% confidence intervals, I-squared
#     statistics, and Cochran's Q values for pre-implementation,
#     post-implementation, and DiD estimates
#
# Outcomes:
#   - All-cause admissions
#   - Cardiovascular admissions
#   - Infectious and parasitic disease admissions
#   - Respiratory admissions
#   - Urogenital admissions
#
# Model specifications:
#   Models 0, 1, and 2 using the RFM5 alert definition
#
# Post-period estimates:
#   Reference-switch estimates are used when available. For cardiovascular
#   admissions, the preserved script approximates the post-period association
#   by combining the district-specific pre-period and DiD log estimates and
#   their variances.
#
# Required packages:
#   dplyr, readxl, mixmeta
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized.
# ==============================================================================

rm(list=ls())

library(dplyr)
library(readxl)
library(mixmeta)

print("Lese Standard- und Reference-Switch-Dateien ein...")

# =====================================================================
# 1. Pfade definieren
# =====================================================================

file_reg_all <-
  "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"

file_reg_strat <- paste0(
  "data/processed/fdz_outputs/",
  "02_Analysis_FirstStage_default_allCause_and_Strat",
  "(mit_Card)_und_sens_Res_gh.xlsx"
)

file_ref_switch <-
  "data/processed/fdz_outputs/Results_FirstStage_ReferenceSwitch_gh.xlsx"

# =====================================================================
# 2A. PRE- UND DiD-EFFEKTE (Aus Standard-Analysen für ALLE Diagnosen)
# =====================================================================

df_pre_did <- bind_rows(read_excel(file_reg_all), read_excel(file_reg_strat)) %>%
  filter(RFM == "heat_alerts_rfm5") %>%
  rename(Diagnosis = any_of(c("Diagnosis", "Diagnose", "Diagnosegruppe", "diagnosis"))) %>%
  mutate(
    Diagnosis = case_when(
      Diagnosis == "AllAdmissions" ~ "All-Cause",
      Diagnosis == "Card_0" ~ "Cardiovascular",
      Diagnosis == "Infect_par" ~ "Infectious & Parasitic",
      Diagnosis == "Resp_0" ~ "Respiratory",
      Diagnosis == "Uri_0" ~ "Urogenital",
      TRUE ~ Diagnosis
    ),
    AGS_Clean = sprintf("%05d", as.numeric(gsub("AGS_|\\.csv", "", AGS_File))),
    
    RR_Alert   = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Up_Alert  = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
    
    RR_DiD     = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Up_DiD  = as.numeric(gsub(",", ".", as.character(CI_Up_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD)))
  ) %>%
  select(AGS_Clean, Diagnosis, Model, RR_Alert, CI_Low_Alert, CI_Up_Alert, RR_DiD, CI_Low_DiD, CI_Up_DiD)

# =====================================================================
# 2B. POST-EFFEKTE AUS DEM ECHTEN REFERENCE SWITCH (Außer Cardio)
# =====================================================================

df_post_switch <- read_excel(file_ref_switch) %>%
  filter(RFM == "heat_alerts_rfm5") %>%
  rename(Diagnosis = any_of(c("Diagnosis", "Diagnose", "Diagnosegruppe", "diagnosis"))) %>%
  mutate(
    Diagnosis = case_when(
      Diagnosis == "AllAdmissions" ~ "All-Cause",
      Diagnosis == "Card_0" ~ "Cardiovascular", 
      Diagnosis == "Infect_par" ~ "Infectious & Parasitic",
      Diagnosis == "Resp_0" ~ "Respiratory",
      Diagnosis == "Uri_0" ~ "Urogenital",
      TRUE ~ Diagnosis
    ),
    AGS_Clean = sprintf("%05d", as.numeric(gsub("AGS_|\\.csv", "", AGS_File))),
    
    RR_Post     = as.numeric(gsub(",", ".", as.character(RR_Alert_Post2005))),
    CI_Low_Post = as.numeric(gsub(",", ".", as.character(CI_Low_Alert_Post2005))),
    CI_Up_Post  = as.numeric(gsub(",", ".", as.character(CI_Up_Alert_Post2005)))
  ) %>%
  select(AGS_Clean, Diagnosis, Model, RR_Post, CI_Low_Post, CI_Up_Post)

# =====================================================================
# 3. HILFSFUNKTION FÜR DAS EXAKTE POOLING
# =====================================================================

get_pooled_results <- function(rr_val, ci_l, ci_u, model_name, effect_name, diag_name) {
  yi <- log(rr_val)
  sei <- (log(ci_u) - log(ci_l)) / (2 * 1.96)
  vi <- sei^2
  
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  if (length(valid_idx) < 5) return(NULL)
  
  meta_model <- mixmeta(yi[valid_idx] ~ 1, S = vi[valid_idx], method = "reml")
  
  RRpooled <- exp(predict(meta_model, ci=TRUE))
  summ <- summary(meta_model)
  
  i2_wert <- paste0(round(max(summ$i2stat), 1), " %")
  q_wert  <- round(summ$qstat$Q, 1)
  
  rr_str <- sprintf("%.3f (%.3f, %.3f)", round(RRpooled[1, "fit"], 3), round(RRpooled[1, "ci.lb"], 3), round(RRpooled[1, "ci.ub"], 3))
  
  data.frame(
    Diagnosis = diag_name,
    Model = model_name,
    Estimate = effect_name,
    `RR (95% CI)` = rr_str,
    `I² (%)` = i2_wert,
    `Cochran's Q` = as.character(q_wert),
    check.names = FALSE
  )
}

# =====================================================================
# 4. SCHLEIFE ÜBER DIAGNOSEN UND MODELLE
# =====================================================================

print("Starte Meta-Analysen für Pre, Post und DiD...")
results_list <- list()
diagnosen_liste <- c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital")
modelle_liste   <- c("model0", "model1", "model2")

for (diag in diagnosen_liste) {
  for (m in modelle_liste) {
    
    df_sub_pre_did <- df_pre_did %>% filter(Diagnosis == diag, Model == m)
    if(nrow(df_sub_pre_did) < 5) next
    
    # 1. PRE-EFFEKT (Standard)
    
    res_pre <- get_pooled_results(
      df_sub_pre_did$RR_Alert, df_sub_pre_did$CI_Low_Alert, df_sub_pre_did$CI_Up_Alert, 
      m, "Pre-implementation alert-day association", diag
    )
    if(!is.null(res_pre)) results_list[[length(results_list) + 1]] <- res_pre
    
    # 2. DiD-EFFEKT (Standard)
    
    res_did <- get_pooled_results(
      df_sub_pre_did$RR_DiD, df_sub_pre_did$CI_Low_DiD, df_sub_pre_did$CI_Up_DiD, 
      m, "Difference-in-differences", diag
    )
    if(!is.null(res_did)) results_list[[length(results_list) + 1]] <- res_did
    
    # 3. POST-EFFEKT (Unterscheidung Cardio-Hack vs. Rest)
    
    if (diag == "Cardiovascular") {
      
      # --- OPTION B HACK NUR FÜR CARDIO (Pre + DiD) ---
      
      yi_pre <- log(df_sub_pre_did$RR_Alert)
      vi_pre <- ((log(df_sub_pre_did$CI_Up_Alert) - log(df_sub_pre_did$CI_Low_Alert)) / (2 * 1.96))^2
      
      yi_did <- log(df_sub_pre_did$RR_DiD)
      vi_did <- ((log(df_sub_pre_did$CI_Up_DiD) - log(df_sub_pre_did$CI_Low_DiD)) / (2 * 1.96))^2
      
      yi_post_hack <- yi_pre + yi_did
      vi_post_hack <- vi_pre + vi_did
      
      rr_hack <- exp(yi_post_hack)
      ci_l_hack <- exp(yi_post_hack - 1.96 * sqrt(vi_post_hack))
      ci_u_hack <- exp(yi_post_hack + 1.96 * sqrt(vi_post_hack))
      
      res_post <- get_pooled_results(
        rr_hack, ci_l_hack, ci_u_hack, 
        m, "Post-implementation alert-day association (Hack)", diag
      )
      if(!is.null(res_post)) results_list[[length(results_list) + 1]] <- res_post
      
    } else {
      
      # --- ECHTER POST EFFEKT AUS DEM SWITCH FÜR ALLE ANDEREN ---
      
      df_sub_post <- df_post_switch %>% filter(Diagnosis == diag, Model == m)
      
      if(nrow(df_sub_post) > 0) {
        df_matched <- inner_join(df_sub_pre_did %>% select(AGS_Clean), df_sub_post, by = "AGS_Clean")
        
        res_post <- get_pooled_results(
          df_matched$RR_Post, df_matched$CI_Low_Post, df_matched$CI_Up_Post, 
          m, "Post-implementation alert-day association", diag
        )
        if(!is.null(res_post)) results_list[[length(results_list) + 1]] <- res_post
      }
    }
    
  }
}

# =====================================================================
# 5. EXPORT DER TABELLE
# =====================================================================

df_final <- bind_rows(results_list) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = diagnosen_liste),
    Estimate = factor(Estimate, levels = c(
      "Pre-implementation alert-day association",
      "Post-implementation alert-day association",
      "Post-implementation alert-day association (Hack)",
      "Difference-in-differences"
    ))
  ) %>%
  arrange(Diagnosis, Model, Estimate)

print(as.data.frame(df_final))

csv_filename <- paste0(
  "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
  "07_pre_post_did_table/Table_Final_Results_ALL_DIAGNOSES_PrePostDiD.csv"
)

write.csv2(
  df_final,
  csv_filename,
  row.names = FALSE
)

print(paste("--> Tabelle erfolgreich gespeichert unter:", csv_filename))
