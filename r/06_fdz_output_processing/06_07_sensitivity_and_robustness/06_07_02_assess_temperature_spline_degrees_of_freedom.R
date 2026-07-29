# ==============================================================================
# Script: 06_07_02_assess_temperature_spline_degrees_of_freedom.R
#
# Purpose:
#   Assesses the robustness of pooled alert-day and Difference-in-Differences
#   estimates to alternative degrees of freedom for the temperature-admission
#   response function. District-specific estimates based on 3, 4, and 5 degrees
#   of freedom are pooled separately by diagnosis.
#
# Inputs:
#   - Released first-stage FDZ sensitivity output containing estimates for
#     temperature functions with 3, 4, and 5 degrees of freedom
#
# Outputs:
#   - PDF plot comparing pooled estimates across degrees-of-freedom settings
#   - CSV table containing pooled RRs and 95% confidence intervals
#   - Numerical summaries printed to the console
#
# Outcomes:
#   Diagnosis groups contained in the first-stage sensitivity output
#
# Model specification:
#   RFM5 with separate pooling for temperature-function DF values 3, 4, and 5
#
# Required packages:
#   dplyr, readxl, mixmeta, ggplot2, tidyr
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

rm(list=ls())

library(dplyr)
library(readxl)
library(mixmeta)
library(ggplot2)
library(tidyr)

# =====================================================================
# 1. DATEN EINLESEN UND VARIANZEN BERECHNEN
# =====================================================================

print("Lese Daten für DF-Sensitivitätsanalyse ein...")

# Datei einlesen (Pfad ggf. anpassen)
file_df <- paste0(
  "data/processed/fdz_outputs/",
  "02_Analysis_FirstStage_default_allCause_and_Strat",
  "(mit_Card)_und_sensi(df)_gh.xlsx"
)
df_raw <- read_excel(file_df)

# Bereinigen und Varianzen aus den CIs zurückrechnen
df_clean <- df_raw %>%
  filter(RFM == "heat_alerts_rfm5") %>% # 
  mutate(
    Diagnosis = case_when(
      Diagnosis == "AllAdmissions" ~ "All-Cause",
      Diagnosis == "Card_0" ~ "Cardiovascular Diseases",
      Diagnosis == "Infect_par" ~ "Infectious & Parasitic Diseases",
      Diagnosis == "Resp_0" ~ "Respiratory Diseases",
      Diagnosis == "Uri_0" ~ "Urogenital Diseases",
      TRUE ~ Diagnosis
    ),
    Tested_DF = as.character(Tested_DF),
    DF_Label = paste0("DF = ", Tested_DF),
    
    RR_Alert = as.numeric(RR_Alert),
    CI_Low_Alert = as.numeric(CI_Low_Alert),
    CI_Up_Alert = as.numeric(CI_Up_Alert),
    
    RR_DiD = as.numeric(RR_DiD),
    CI_Low_DiD = as.numeric(CI_Low_DiD),
    CI_Up_DiD = as.numeric(CI_Up_DiD),
    
    # Rückrechnung der Varianzen für ALERT
    log_RR_Alert = log(RR_Alert),
    SE_Alert = (log(CI_Up_Alert) - log(CI_Low_Alert)) / (2 * 1.96),
    Var_Alert = SE_Alert^2,
    
    # Rückrechnung der Varianzen für DiD
    log_RR_DiD = log(RR_DiD),
    SE_DiD = (log(CI_Up_DiD) - log(CI_Low_DiD)) / (2 * 1.96),
    Var_DiD = SE_DiD^2
  )

diagnosen_liste <- unique(df_clean$Diagnosis)
df_liste <- c("3", "4", "5")

# =====================================================================
# 2. META-ANALYSE ÜBER ALLE DF UND DIAGNOSEN SCHLEIFEN
# =====================================================================

print("Starte Meta-Analysen für DF 3, 4 und 5...")

results_list <- list()

for (diag in diagnosen_liste) {
  for (aktuelles_df in df_liste) {
    
    df_sub <- df_clean %>% 
      filter(Diagnosis == diag, Tested_DF == aktuelles_df) %>%
      filter(
        !is.na(log_RR_Alert) & !is.infinite(log_RR_Alert),
        !is.na(Var_Alert) & !is.infinite(Var_Alert) & Var_Alert > 0,
        !is.na(log_RR_DiD) & !is.infinite(log_RR_DiD),
        !is.na(Var_DiD) & !is.infinite(Var_DiD) & Var_DiD > 0
      )
    
    if(nrow(df_sub) < 5) {
      warning(paste("Überspringe", diag, "mit DF", aktuelles_df, "- zu wenig Daten!"))
      next
    }
    
    # Pooling ALERT
    meta_alert <- mixmeta(log_RR_Alert ~ 1, S = Var_Alert, data = df_sub, method = "reml")
    alert_coef <- coef(meta_alert)[1]
    alert_se <- sqrt(vcov(meta_alert)[1,1])
    
    # Pooling DID
    meta_did <- mixmeta(log_RR_DiD ~ 1, S = Var_DiD, data = df_sub, method = "reml")
    did_coef <- coef(meta_did)[1]
    did_se <- sqrt(vcov(meta_did)[1,1])
    
    results_list[[paste0(diag, "_DF", aktuelles_df, "_Alert")]] <- data.frame(
      Diagnosis = diag,
      DF_Label = paste0("DF = ", aktuelles_df),
      Effect_Type = "Pre-Implementation Alert-Day Association", 
      RR = exp(alert_coef),
      CI_low = exp(alert_coef - 1.96 * alert_se),
      CI_high = exp(alert_coef + 1.96 * alert_se)
    )
    
    results_list[[paste0(diag, "_DF", aktuelles_df, "_DiD")]] <- data.frame(
      Diagnosis = diag,
      DF_Label = paste0("DF = ", aktuelles_df),
      Effect_Type = "Change in Alert-Day Association (DiD)",
      RR = exp(did_coef),
      CI_low = exp(did_coef - 1.96 * did_se),
      CI_high = exp(did_coef + 1.96 * did_se)
    )
  }
}

plot_data <- bind_rows(results_list)

# Precision-Status festlegen (für die Rauten)
plot_data <- plot_data %>%
  mutate(
    Precision_Status = case_when(
      is.na(CI_low) | is.na(CI_high) ~ "Unknown", 
      CI_low > 1 | CI_high < 1 ~ "Excludes 1.0",
      TRUE ~ "Includes 1.0"
    ),
    Precision_Status = factor(Precision_Status, levels = c("Excludes 1.0", "Includes 1.0", "Unknown")),
    
    # Y-Achse sortieren (DF 5 oben, DF 3 unten)
    DF_Label = factor(DF_Label, levels = rev(c("DF = 3", "DF = 4", "DF = 5"))),
    
    # METHODISCH SAUBERE FAKTOREN
    Effect_Type = factor(Effect_Type, levels = c("Pre-Implementation Alert-Day Association", "Change in Alert-Day Association (DiD)"))
  )

# =====================================================================
# 3. PUBLICATION PLOT FÜR DF SENSITIVITÄT
# =====================================================================

print("Generiere Sensitivitäts-Plot...")

p_df_sens <- ggplot(plot_data, aes(x = RR, y = DF_Label, color = Effect_Type, group = Effect_Type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.3, linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = Precision_Status), size = 4, stroke = 1.2, position = position_dodge(width = 0.5)) + 
  
  facet_wrap(~ Diagnosis, ncol = 2) +
  
  scale_color_manual(
    values = c("Pre-Implementation Alert-Day Association" = "#1f78b4", "Change in Alert-Day Association (DiD)" = "#e31a1c")
  ) +
  scale_shape_manual(values = c("Excludes 1.0" = 18, "Includes 1.0" = 5, "Unknown" = 4)) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Sensitivity to Temperature–Response Spline Flexibility",
    subtitle = "National pooled cause-specific estimates using 3–5 degrees of freedom for the temperature–admission function",
    x = "Relative Risk (RR) with 95% Confidence Interval",
    y = "Degrees of Freedom",
    color = "Estimate",
    shape = "95% CI"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical", 
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

pdf_filename <- paste0(
  "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
  "02_temperature_df/Sensitivity_DF_3_4_5_English.pdf"
)
# BILD-DIMENSIONEN ANGEPASST FÜR 2-SPALTIGES LAYOUT:
ggsave(pdf_filename, plot = p_df_sens, width = 10, height = 10)
print(paste("--> Plot erfolgreich gespeichert als:", pdf_filename))

# =====================================================================
# 4. EXAKTE WERTE (RR UND 95%-KI) FÜR DEN TEXT EXPORTIEREN
# =====================================================================

print("Erstelle Tabelle mit exakten RRs und CIs für den Ergebnisteil...")

df_export <- plot_data %>%
  mutate(
    RR_rounded = sprintf("%.2f", RR),
    CI_low_rounded = sprintf("%.2f", CI_low),
    CI_high_rounded = sprintf("%.2f", CI_high),
    Estimate_String = paste0(RR_rounded, " [", CI_low_rounded, ", ", CI_high_rounded, "]"),
    
    Term = case_when(
      Effect_Type == "Pre-Implementation Alert-Day Association" ~ "Pre-Implementation Association",
      Effect_Type == "Change in Alert-Day Association (DiD)" ~ "DiD Shift",
      TRUE ~ as.character(Effect_Type)
    )
  ) %>%
  select(Diagnosis, DF_Label, Term, Estimate_String)

df_export_wide <- df_export %>%
  pivot_wider(
    names_from = Term, 
    values_from = Estimate_String
  ) %>%
  arrange(Diagnosis, DF_Label)

csv_filename <- paste0(
  "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
  "02_temperature_df/Table_Sensitivity_DF_3_4_5_Estimates.csv"
)
write.csv2(
  df_export_wide,
  csv_filename,
  row.names = FALSE
)
print(paste("--> Tabelle erfolgreich gespeichert als:", csv_filename))

# =====================================================================
# 4. EXAKTE WERTE (RR UND 95%-KI) FÜR DEN TEXT EXPORTIEREN
# =====================================================================

print("Erstelle Tabelle mit exakten RRs und CIs für den Ergebnisteil...")

# 1. Daten formatieren: RR [95% CI] in einem schönen String zusammenfassen
df_export <- plot_data %>%
  mutate(
    # Werte auf 2 Nachkommastellen runden
    RR_rounded = sprintf("%.2f", RR),
    CI_low_rounded = sprintf("%.2f", CI_low),
    CI_high_rounded = sprintf("%.2f", CI_high),
    
    # Publikationsreifer String: "1.05 [0.98, 1.12]"
    Estimate_String = paste0(RR_rounded, " [", CI_low_rounded, ", ", CI_high_rounded, "]"),
    
    # Die neuen sauberen Labels für den Text
    Term = case_when(
      Effect_Type == "Baseline Risk (Alert-Qualifying Days)" ~ "Pre-Implementation Association",
      Effect_Type == "Added System Effect (DiD)" ~ "DiD Shift",
      TRUE ~ as.character(Effect_Type)
    )
  ) %>%
  # FIX: Wir greifen auf DF_Label zu, weil Tested_DF nicht im plot_data existiert
  select(Diagnosis, DF_Label, Term, Estimate_String)

# 2. Die Tabelle ins "Wide Format" bringen
df_export_wide <- df_export %>%
  pivot_wider(
    names_from = Term, 
    values_from = Estimate_String
  ) %>%
  arrange(Diagnosis, DF_Label)

# 3. CSV exportieren
csv_filename <- paste0(
  "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
  "02_temperature_df/Table_Sensitivity_DF_3_4_5_Estimates.csv"
)
write.csv2(
  df_export_wide,
  csv_filename,
  row.names = FALSE
)
print(paste("--> Tabelle erfolgreich gespeichert als:", csv_filename))

# 4. Vorschau für die Console
print("Vorschau:")
print(df_export_wide)
