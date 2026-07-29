# ==============================================================================
# Script: 06_07_04_assess_temperature_alert_overlap.R
#
# Purpose:
#   Assesses temperature overlap and positivity for alert and non-alert
#   district-days in the pre-implementation and post-implementation periods.
#   Alert probabilities are summarized within 1-degree temperature bins, and
#   overlap coefficients are calculated nationally and separately by district.
#
# Inputs:
#   - Environmental and heat-alert master dataset containing district-level
#     daily temperature and RFM5 alert classifications
#
# Outputs:
#   - PDF positivity plot showing alert probability and district-day density
#     across temperature bins
#   - National overlap coefficients printed to the console
#   - District-specific overlap summaries printed to the console
#
# Analysis period:
#   May–September, 2000–2009
#
# Temperature grouping:
#   One-degree Celsius bins
#
# Required packages:
#   dplyr, ggplot2, lubridate, tidyr
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

rm(list=ls())

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

print("Lese Masterfile für Overlap-Diagnostik ein...")

# =====================================================================
# 1. Masterfile einlesen
# =====================================================================

df_master <- read.csv(
  "data/external/fdz_input/Masterfile_Final_for_FDZ.csv",
  sep = ","
)

# =====================================================================
# 2. Daten bereinigen, filtern und Perioden zuweisen
# =====================================================================

df_overlap <- df_master %>%
  mutate(
    Datum = as.Date(Datum, format="%Y-%m-%d"),
    Jahr = year(Datum),
    Monat = month(Datum),
    
    # Alert Status numerisch für einfache Anteilsberechnung
    Alert_Num = case_when(
      heat_alerts_rfm5 == 1 ~ 1,
      heat_alerts_rfm5 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Pre vs. Post Periode (inkl. Labeling für den Plot)
    Period = case_when(
      Jahr >= 2000 & Jahr <= 2004 ~ "Pre-Implementation (2000-2004)\nReconstructed Alerts",
      Jahr >= 2005 & Jahr <= 2009 ~ "Post-Implementation (2005-2009)\nOfficial Alerts",
      TRUE ~ NA_character_
    ),
    
    # 1-Grad Temperatur Bins (Floor: 15.6 wird zu 15)
    Temp_Bin = floor(T_mean_popw)
  ) %>%
  filter(Jahr >= 2000 & Jahr <= 2009, Monat >= 5 & Monat <= 9) %>%
  filter(!is.na(T_mean_popw) & !is.na(Alert_Num) & !is.na(Period)) %>%
  mutate(Period = factor(Period, levels = c(
    "Pre-Implementation (2000-2004)\nReconstructed Alerts",
    "Post-Implementation (2005-2009)\nOfficial Alerts"
  )))

# =====================================================================
# 3. DATENAGGREGATION FÜR DEN POSITIVITY PLOT
# =====================================================================

# Aggregation pro Periode und Temp-Bin über ALLE Landkreise
df_plot <- df_overlap %>%
  group_by(Period, Temp_Bin) %>%
  summarise(
    n_total = n(),
    n_alert = sum(Alert_Num),
    prop_alert = n_alert / n_total,
    .groups = "drop"
  ) %>%
  # Nur Bins anzeigen, in denen überhaupt nennenswerte Daten vorhanden sind
  filter(n_total >= 10)

# =====================================================================
# 4. DER POSITIVITY PLOT
# =====================================================================

# Kombinierter Plot: Anteil Alert-Tage (Punkte/Linie) + Gesamtmasse (Balken im Hintergrund)

max_n <- max(df_plot$n_total)

p_positivity <- ggplot(df_plot, aes(x = Temp_Bin)) +
  # Hintergrund: Histogramm der absoluten Fallzahlen (Grau)
  geom_col(aes(y = n_total / max_n), fill = "grey80", width = 0.8, alpha = 0.6) +
  
  # Vordergrund: Propensity of Alert (Rote Linie/Punkte)
  geom_line(aes(y = prop_alert), color = "#d73027", linewidth = 1) +
  geom_point(aes(y = prop_alert, size = n_total), color = "#a50026", alpha = 0.8) +
  
  # Faceting nach Pre/Post
  facet_wrap(~ Period) +
  
  # Doppel-Achse: Links Probability, Rechts Count
  scale_y_continuous(
    name = "Probability of Alert Day",
    labels = scales::percent_format(accuracy = 1),
    sec.axis = sec_axis(~ . * max_n, name = "Total District-Days (Count)")
  ) +
  scale_x_continuous(breaks = seq(min(df_plot$Temp_Bin), max(df_plot$Temp_Bin), by = 2)) +
  
  labs(
    title = "Positivity Analysis: Alert Probability by Temperature",
    subtitle = "Points show proportion of alerts; grey bars and point size indicate common support density",
    x = "Mean Temperature Bin (°C)",
    size = "Total Days\nin Bin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold"),
    axis.title.y.left = element_text(color = "#a50026", face = "bold"),
    axis.text.y.left = element_text(color = "#a50026"),
    axis.title.y.right = element_text(color = "grey50", face = "bold"),
    axis.text.y.right = element_text(color = "grey50")
  )

ggsave(
  paste0(
    "outputs/06_fdz_output_processing/06_07_sensitivity_and_robustness/",
    "04_overlap/Overlap_Diagnostic_Positivity.pdf"
  ),
  plot = p_positivity,
  width = 12,
  height = 6
)
print("--> Positivity Plot gespeichert unter: Overlap_Diagnostic_Positivity.pdf")

# =====================================================================
# 5. BERECHNUNG DES OVERLAP KOEFFIZIENTEN (OVL)
# =====================================================================

print("Berechne Overlap-Koeffizienten...")

# Hilfsfunktion zur OVL Berechnung
calculate_ovl <- function(df) {
  
  # Gesamtanzahl Alerts und Non-Alerts im DF
  total_alerts <- sum(df$Alert_Num == 1)
  total_non_alerts <- sum(df$Alert_Num == 0)
  
  # Fallback, falls keine Alerts/Non-Alerts existieren (OVL = 0)
  if(total_alerts == 0 | total_non_alerts == 0) return(0)
  
  # Verteilungen (p_b) pro Bin berechnen
  bin_dist <- df %>%
    group_by(Temp_Bin) %>%
    summarise(
      n_alert = sum(Alert_Num == 1),
      n_non = sum(Alert_Num == 0),
      .groups = "drop"
    ) %>%
    mutate(
      p_alert = n_alert / total_alerts,
      p_non = n_non / total_non_alerts,
      # Minimum aus beiden Verteilungen (Gemeinsamer Anteil)
      min_p = pmin(p_alert, p_non)
    )
  
  # Summe der Minima ist der OVL
  sum(bin_dist$min_p)
}

# --- A) G gepoolter OVL (Gesamtdeutschland, Pre vs Post) ---

ovl_pooled <- df_overlap %>%
  group_by(Period) %>%
  group_modify(~ data.frame(OVL = calculate_ovl(.x))) %>%
  ungroup()

print("--- G gepoolter Overlap-Koeffizient (OVL) nach Periode ---")
print(ovl_pooled)

# --- B) Kreisspezifischer OVL (Median & IQR) ---

ovl_kreis <- df_overlap %>%
  group_by(Period, AGS) %>%
  group_modify(~ data.frame(OVL = calculate_ovl(.x))) %>%
  ungroup()

ovl_kreis_summary <- ovl_kreis %>%
  group_by(Period) %>%
  summarise(
    Median_OVL = median(OVL, na.rm = TRUE),
    Q1_OVL = quantile(OVL, 0.25, na.rm = TRUE),
    Q3_OVL = quantile(OVL, 0.75, na.rm = TRUE),
    IQR_OVL = IQR(OVL, na.rm = TRUE),
    Min_OVL = min(OVL, na.rm = TRUE),
    Max_OVL = max(OVL, na.rm = TRUE),
    .groups = "drop"
  )

print("--- Kreisspezifischer Overlap-Koeffizient (Zusammenfassung) ---")
print(ovl_kreis_summary)
