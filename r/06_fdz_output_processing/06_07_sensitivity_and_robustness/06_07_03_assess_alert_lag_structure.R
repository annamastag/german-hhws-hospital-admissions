# ==============================================================================
# Script: 06_07_03_assess_alert_lag_structure.R
#
# Purpose:
#   Describes the temporal lag structure of district-specific heat-alert
#   associations from the alert day through three subsequent days. The median
#   lag-specific relative risk is calculated across districts separately for
#   each diagnosis group.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Released cause-specific first-stage FDZ output table including
#     cardiovascular admissions
#
# Outputs:
#   - PDF plot of median alert-day relative risks for lags 0–3 by diagnosis
#   - Lag-specific numerical summaries printed to the console
#
# Outcomes:
#   - All-cause admissions
#   - Cardiovascular admissions
#   - Infectious and parasitic disease admissions
#   - Respiratory admissions
#   - Urogenital admissions
#
# Model specification:
#   RFM5 and Model 1
#
# Required packages:
#   dplyr, tidyr, ggplot2, readxl
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

print("Lese Daten für Lag-Analyse ein...")

# 1. Daten einlesen 

df_all_cause <- read_excel(
  "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
)
df_diagnosen <- read_excel(
  paste0(
    "data/processed/fdz_outputs/",
    "02_Analysis_FirstStage_default_allCause_and_Strat",
    "(mit_Card)_und_sens_Res_gh.xlsx"
  )
)
df_combined <- bind_rows(df_all_cause, df_diagnosen)

# 2. Daten bereinigen und Lags extrahieren

df_lags <- df_combined %>%
  filter(Model == "model1", RFM == "heat_alerts_rfm5") %>%
  mutate(
    Diagnosis = case_when(
      Diagnosis == "AllAdmissions" ~ "All-Cause",
      Diagnosis == "Card_0" ~ "Cardiovascular",
      Diagnosis == "Infect_par" ~ "Infectious & Parasitic",
      Diagnosis == "Resp_0" ~ "Respiratory",
      Diagnosis == "Uri_0" ~ "Urogenital",
      TRUE ~ Diagnosis
    )
  ) %>%
  
  # Nur die relevanten Lag-Spalten behalten
  select(AGS_File, Diagnosis, RR_Alert_Lag0, RR_Alert_Lag1, RR_Alert_Lag2, RR_Alert_Lag3) %>%
  
  # Kommas in Punkte umwandeln und zu numerisch machen
  mutate(across(starts_with("RR_Alert_Lag"), ~ as.numeric(gsub(",", ".", as.character(.))))) %>%
  
  # Pivot Longer: Tabelle umbauen zum plotten
  pivot_longer(
    cols = starts_with("RR_Alert_Lag"),
    names_to = "Lag",
    values_to = "RR"
  ) %>%
  mutate(
    # Namen für die X-Achse vergeben
    Lag = case_when(
      Lag == "RR_Alert_Lag0" ~ "Lag 0 (Same Day)",
      Lag == "RR_Alert_Lag1" ~ "Lag 1 (Day +1)",
      Lag == "RR_Alert_Lag2" ~ "Lag 2 (Day +2)",
      Lag == "RR_Alert_Lag3" ~ "Lag 3 (Day +3)"
    ),
    Lag = factor(Lag, levels = c("Lag 0 (Same Day)", "Lag 1 (Day +1)", "Lag 2 (Day +2)", "Lag 3 (Day +3)")),
    Diagnosis = factor(Diagnosis, levels = c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital"))
  ) %>%
  filter(!is.na(RR))

# 3. Aggregation: Den Median pro Diagnose und Lag-Tag berechnen

df_lag_summary <- df_lags %>%
  group_by(Diagnosis, Lag) %>%
  summarise(
    Median_RR = median(RR, na.rm = TRUE),
    .groups = "drop"
  )

# 4. ZAHLEN IN DIE KONSOLE AUSGEBEN

cat("\n=== ZEITLICHER VERLAUF: MEDIAN ALERT EFFECT (RR) NACH LAGS ===\n")
print(df_lag_summary %>% 
        pivot_wider(names_from = Lag, values_from = Median_RR) %>%
        mutate(across(where(is.numeric), ~ round(., 3))))
cat("================================================================\n")

# 5. PUBLIKATIONS-PLOT ERSTELLEN

p_lags <- ggplot(df_lag_summary, aes(x = Lag, y = Median_RR, group = Diagnosis, color = Diagnosis)) +
  
  # Null-Effekt-Linie einzeichnen
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  
  # Linien und Punkte
  geom_line(linewidth = 1.2) +
  geom_point(size = 4, shape = 18) +
  
  # Facetten für jede Diagnose (Matrix)
  facet_wrap(~ Diagnosis, scales = "fixed", ncol = 3) +
  
  scale_color_manual(values = c("All-Cause" = "black", 
                                "Cardiovascular" = "#e41a1c", 
                                "Infectious & Parasitic" = "#4daf4a", 
                                "Respiratory" = "#377eb8", 
                                "Urogenital" = "#984ea3")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Temporal Dynamics of Heat Alerts (Lag Structure)",
    subtitle = "Median Relative Risk across 400 districts following an alert (Lag 0 to 3)",
    x = "Time Since Alert",
    y = "Median Relative Risk (RR)"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

pdf_filename <- paste0(
  "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
  "03_alert_lags/Plot_Alert_Lag_Structure_Dynamics_EN.pdf"
)
ggsave(
  pdf_filename,
  plot = p_lags,
  width = 12,
  height = 7
)
print(paste("--> Lag-Plot erfolgreich gespeichert als:", pdf_filename))
