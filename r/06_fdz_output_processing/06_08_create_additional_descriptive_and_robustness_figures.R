# ==============================================================================
# Script: 06_08_create_additional_descriptive_and_robustness_figures.R
#
# Purpose:
#   Creates additional descriptive, model-evaluation, regional-heterogeneity,
#   subgroup, and temporal-comparison figures for the study. The preserved
#   script contains several sequential reporting and exploratory analysis
#   blocks using Random Forest evaluation data, released FDZ first-stage
#   estimates, district-level environmental data, population data, and
#   harmonized district boundaries.
#
# Inputs:
#   - Environmental and heat-alert master dataset
#   - Released all-cause first-stage FDZ output table
#   - Released cause-specific first-stage FDZ output table
#   - Released reference-switch first-stage FDZ output table
#   - Random Forest validation output for 2011–2022
#   - List of districts selected in the exploratory heat-sensitivity analysis
#   - Population raster derived from the 2011 German Census
#   - Harmonized German district boundary shapefile
#
# Outputs:
#   - Three-panel descriptive map of temperature, alert days, and admissions
#   - Spatial Random Forest evaluation figures
#   - National and cause-specific pooled-effect figures
#   - Subgroup figures for selected heat-sensitive districts
#   - Federal-state pooled estimates and spatial DiD figures
#   - Pre-implementation versus post-implementation comparison figures
#   - Pooled temporal temperature-response comparison figure
#
# Analysis scope:
#   The script contains multiple sequential analysis and figure-generation
#   blocks. Some blocks use all-cause admissions, while others also include
#   infectious and parasitic, respiratory, and urinary admissions.
#
# Required packages:
#   sf, terra, exactextractr, dplyr, ggplot2, readxl, lubridate, patchwork,
#   scales, tidyverse, viridis, mixmeta, dlnm
#
# Repository note:
#   Outputs are written to topic-specific subdirectories. The file name,
#   path strings, and explanatory comments were standardized.
# ==============================================================================

library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(ggplot2)
library(readxl)
library(lubridate)
library(patchwork) 
library(scales) 
library(tidyverse)
library(viridis)
library(mixmeta)

# ---------------------------------------------------------------------
# 3-PANEL DEUTSCHLANDKARTE (Temperatur | Heiße Tage | Total Cases Rate)
# ---------------------------------------------------------------------

shapefile_path   <- "data/raw/spatial/vg2500/vg2500.shp"
file_temp_alerts <- "data/external/fdz_input/Masterfile_Final_for_FDZ.csv"
file_cases       <- "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
path_pop         <- "data/raw/population/census_2011/population_1km.tif"

# 1. SHAPEFILE & BEVÖLKERUNG LADEN UND BERECHNEN

map_de <- st_read(shapefile_path, quiet = TRUE)
map_de$AGS <- sprintf("%05d", as.numeric(map_de$AGS))

pop_raster <- rast(path_pop)

# CRS-Check 
if (st_crs(map_de)$wkt != crs(pop_raster)) {
  map_de <- st_transform(map_de, crs(pop_raster))
}

# Absolute Einwohnerzahl pro AGS aufsummieren
map_de$Population_2011 <- exact_extract(pop_raster, map_de, fun = "sum", progress = FALSE)

# 2. DATEN VORBEREITEN: TEMPERATUR & ALERTS (2000 - 2009)

df_csv <- read.csv(file_temp_alerts, sep=",") %>%
  mutate(
    AGS = sprintf("%05d", as.numeric(AGS)),
    Datum = as.Date(Datum, format="%Y-%m-%d"),
    Jahr = year(Datum)
  )

df_agg_csv <- df_csv %>%
  filter(month(Datum) >= 5 & month(Datum) <= 9) %>%
  filter(Jahr >= 2000 & Jahr <= 2009) %>%
  mutate(
    is_alert = ifelse(Jahr < 2005, heat_alerts_rfm5, heat_alerts_dwd)
  ) %>%
  group_by(AGS) %>%
  summarise(
    mean_temp = mean(T_mean_popw, na.rm = TRUE),
    total_alerts = sum(is_alert, na.rm = TRUE)
  )

# 3. DATEN VORBEREITEN: KRANKENHAUSEINWEISUNGEN 

print("Bereite Einweisungs-Daten vor...")
df_cases <- read_excel(file_cases) %>%
  filter(RFM == "heat_alerts_rfm5") %>% 
  mutate(
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = gsub(".csv", "", AGS_Zahlen, fixed = TRUE)
  ) %>%
  distinct(AGS_Zahlen, .keep_all = TRUE) %>%
  select(AGS_Zahlen, Total_Cases) %>%
  mutate(Total_Cases = as.numeric(Total_Cases))

# 4.  ZUSAMMENFÜHREN (Absolute Fallzahlen)

df_all <- df_agg_csv %>%
  left_join(df_cases, by = c("AGS" = "AGS_Zahlen"))

map_merged <- map_de %>%
  left_join(df_all, by = "AGS")

# 5. KARTEN ZEICHNEN (GGPLOT)

theme_map <- theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.width = unit(1.2, "cm"), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# --- KARTE A: Temperatur 
p1 <- ggplot(map_merged) +
  geom_sf(aes(fill = mean_temp), color = "white", linewidth = 0.05) +
  scale_fill_distiller(palette = "Reds", direction = 1, name = "Temp. (°C)") + 
  labs(title = "Average temperature (°C)\n2000-2009 (May-Sept)") +
  theme_map

# --- KARTE B: Heiße Tage / Alerts 
p2 <- ggplot(map_merged) +
  geom_sf(aes(fill = total_alerts), color = "white", linewidth = 0.05) +
  scale_fill_distiller(palette = "YlOrBr", direction = 1, name = "Count") + 
  labs(title = "No. of heat days\n(Pre-2005: RFM, Post-2005: DWD)") +
  theme_map

# --- OPTIMIERTE KARTE C: Absolute Hospital Admissions (Log10 transformiert)
p3 <- ggplot(map_merged) +
  geom_sf(aes(fill = Total_Cases), color = "white", linewidth = 0.05) +
  scale_fill_viridis_c(
    option = "mako", 
    direction = -1, 
    name = "Total Cases",
    trans = "log10", 
    breaks = scales::log_breaks(n = 5), 
    labels = scales::comma_format(big.mark = ",") 
  ) +
  labs(title = "Total hospital admissions\n(Absolute counts, 2000-2009)") +
  theme_map

# 6. PLOTS VERBINDEN & PDF EXPORT

print("Exportiere PDF...")
pdf(
  paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "descriptive_maps/3Panel_Deutschlandkarte_Deskriptiv_EN_absolut_log.pdf"
  ),
  width = 16,
  height = 7
)
combined_plot <- p1 + p2 + p3
print(combined_plot)
dev.off()

# ---------------------------------------------------------
# RÄUMLICHER VERGLEICH: REALE VS. PREDICTETE WARNUNGEN
# ---------------------------------------------------------

# 1. Shapefile der Landkreise laden

lk_shape <- st_read("data/raw/spatial/vg2500/vg2500.shp") 

# 2. Evaluierungsdaten einlesen (Testzeitraum 2011-2022 mit Threshold 0.8)

eval_data <- read.csv(
  paste0(
    "outputs/random_forest_models/03_03_weighted_threshold_comparison/",
    "RFM_Eval_2011-2022_Th0.8.csv"
  ),
  stringsAsFactors = FALSE
)

# 3. Daten aggregieren: Summe der Warntage pro AGS berechnen

warntage_agg <- eval_data %>%
  mutate(
    Real_Numeric = ifelse(Real == "yes", 1, 0),
    Pred_Numeric = ifelse(Pred == "yes", 1, 0)
  ) %>%
  group_by(AGS) %>%
  summarize(
    Actual_DWD = sum(Real_Numeric, na.rm = TRUE),
    Predicted_RFM = sum(Pred_Numeric, na.rm = TRUE)
  ) %>%
  mutate(AGS = as.character(AGS)) 

# CHECK: führende Nullen bei  AGS hinzufügen
warntage_agg <- warntage_agg %>% mutate(AGS = str_pad(AGS, width = 5, side = "left", pad = "0"))

# 4. Shapefile mit den aggregierten Daten verknüpfen

map_sf <- lk_shape %>%
  left_join(warntage_agg, by = c("AGS" = "AGS"))

# NEUER ABSCHNITT: EINHEITLICHE SKALA ERZWINGEN

# 1. Das absolute Maximum über beide Variablen hinweg ermitteln

max_gesamt <- max(c(map_sf$Actual_DWD, map_sf$Predicted_RFM), na.rm = TRUE)

# 2. Plots mit identischen Limits und identischem Skalen-Namen erstellen

# Karte A: Echte Warnungen
plot_actual <- ggplot(map_sf) +
  geom_sf(aes(fill = Actual_DWD), color = NA) +
  scale_fill_viridis_c(option = "plasma", 
                       name = "Total Days", 
                       limits = c(0, max_gesamt)) + 
  theme_void() +
  labs(title = "A) Actual DWD Heat Warnings",
       subtitle = "Total warning days (2011–2022)") +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right")

# Karte B: Vorhergesagte Warnungen
plot_predicted <- ggplot(map_sf) +
  geom_sf(aes(fill = Predicted_RFM), color = NA) +
  scale_fill_viridis_c(option = "plasma", 
                       name = "Total Days", 
                       limits = c(0, max_gesamt)) + 
  theme_void() +
  labs(title = "B) Reconstructed RFM Warnings",
       subtitle = "Total predicted days (Threshold = 0.8)") +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right")

# 3. Beide Karten nebeneinander kombinieren 

kombinierter_plot <- plot_actual + plot_predicted + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right") # gemeinsame Legende rechts

# 4. Speichern

ggsave(
  paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "random_forest_evaluation/RFM_Spatial_Comparison_2011_2022.png"
  ),
  plot = kombinierter_plot,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)


# ---------------------------------------------------------
# BUNDESWEITES POOLING: BASELINE VS. DID EFFEKT FÜR MODEL 0, 1 UND 2
# INKLUSIVE AGGREGIERTEM FOREST PLOT
# ---------------------------------------------------------

rm(list=ls())

# 1. DATEN EINLESEN UND BEREINIGUNG

# FDZ-Dateien
file_all_cause <- "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx" 
file_diagnosen <- "data/processed/fdz_outputs/2026-05-19_Results_FirstStage_Diagnosen_gh.xlsx"

df_all_cause <- read_excel(file_all_cause) 
df_diagnosen <- read_excel(file_diagnosen)

df_combined <- bind_rows(df_all_cause, df_diagnosen)

modelle_auswahl <- c("model0", "model1", "model2")
diagnosen_liste <- c("All-Cause", "Infectious and Parasitic Diseases", "Respiratory Diseases", "Urinary Diseases")

# 2. SCHLEIFE ÜBER DIE MODELLE

for (aktuelles_modell in modelle_auswahl) {
  
  print(paste("=================================================="))
  print(paste("Start Meta-Analyse und Plot für:", aktuelles_modell))
  
  # A) Daten filtern für das aktuelle Modell
  df_clean <- df_combined %>%
    filter(Model == aktuelles_modell, RFM == "heat_alerts_rfm5") %>% 
    mutate(
      Diagnosis = case_when(
        Diagnosis == "AllAdmissions" ~ "All-Cause",
        Diagnosis == "Infect_par" ~ "Infectious and Parasitic Diseases",
        Diagnosis == "Resp_0" ~ "Respiratory Diseases",
        Diagnosis == "Uri_0" ~ "Urinary Diseases",
        TRUE ~ Diagnosis
      ),
      RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
      CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
      CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
      
      log_RR_Alert = log(RR_Alert),
      SE_Alert = (log(CI_Up_Alert) - log(CI_Low_Alert)) / (2 * 1.96),
      Var_Alert = SE_Alert^2,
      
      DiD_Coef = as.numeric(gsub(",", ".", as.character(DiD_Coef))),
      DiD_Variance = as.numeric(gsub(",", ".", as.character(DiD_Variance)))
    )
  
  # B) Meta-Analyse für dieses Modell
  results_list <- list()
  
  for (diag in diagnosen_liste) {
    
    df_sub <- df_clean %>% 
      filter(Diagnosis == diag) %>%
      filter(
        !is.na(log_RR_Alert) & !is.infinite(log_RR_Alert),
        !is.na(Var_Alert) & !is.infinite(Var_Alert) & Var_Alert > 0,
        !is.na(DiD_Coef) & !is.infinite(DiD_Coef),
        !is.na(DiD_Variance) & !is.infinite(DiD_Variance) & DiD_Variance > 0
      )
    
    if(nrow(df_sub) < 5) {
      print(paste("   Überspringe", diag, "- zu wenig gültige Daten!"))
      next
    }
    
    # Pooling ALERT
    meta_alert <- mixmeta(log_RR_Alert ~ 1, S = Var_Alert, data = df_sub, method = "reml")
    alert_coef <- coef(meta_alert)[1]
    alert_se <- sqrt(vcov(meta_alert)[1,1])
    
    # Pooling DID
    meta_did <- mixmeta(DiD_Coef ~ 1, S = DiD_Variance, data = df_sub, method = "reml")
    did_coef <- coef(meta_did)[1]
    did_se <- sqrt(vcov(meta_did)[1,1])
    
    results_list[[paste0(diag, "_Alert")]] <- data.frame(
      Diagnosis = diag,
      Effect_Type = "Alert", 
      RR = exp(alert_coef),
      CI_low = exp(alert_coef - 1.96 * alert_se),
      CI_high = exp(alert_coef + 1.96 * alert_se)
    )
    
    results_list[[paste0(diag, "_DiD")]] <- data.frame(
      Diagnosis = diag,
      Effect_Type = "DiD",
      RR = exp(did_coef),
      CI_low = exp(did_coef - 1.96 * did_se),
      CI_high = exp(did_coef + 1.96 * did_se)
    )
  }
  
  plot_data <- bind_rows(results_list)
  
  # C) Forest Plot für dieses Modell zeichnen
  plot_data <- plot_data %>%
    mutate(
      Precision_Status = case_when(
        is.na(CI_low) | is.na(CI_high) ~ "Unknown", 
        CI_low > 1 | CI_high < 1 ~ "Excludes 1.0",
        TRUE ~ "Includes 1.0"
      ),
      Precision_Status = factor(Precision_Status, levels = c("Excludes 1.0", "Includes 1.0", "Unknown"))
    )
  
  plot_data$Diagnosis <- factor(plot_data$Diagnosis, levels = rev(diagnosen_liste))
  plot_data$Effect_Type <- factor(plot_data$Effect_Type, levels = c("Alert", "DiD"))
  
  titel_anzeige <- toupper(aktuelles_modell) 
  subtitle_anzeige <- paste("Pooled estimates from", titel_anzeige, "| National Average")
  
  p_forest <- ggplot(plot_data, aes(x = RR, y = Diagnosis, color = Effect_Type, group = Effect_Type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.3, linewidth = 1, position = position_dodge(width = 0.5)) +
    
    geom_point(aes(shape = Precision_Status), size = 4, stroke = 1.2, position = position_dodge(width = 0.5)) + 
    
    scale_color_manual(
      values = c("Alert" = "#1f78b4", "DiD" = "#e31a1c"),
      labels = c("Alert" = "Pre-Implementation Alert-Day Association", "DiD" = "Change in Alert-Day Association (DiD)")
    ) +
    
    scale_shape_manual(values = c("Excludes 1.0" = 18, "Includes 1.0" = 5, "Unknown" = 4)) +
    
    theme_bw() +
    labs(
      title = "Cause-Specific Alert-Day Associations and DiD Estimates",
      subtitle = subtitle_anzeige,
      x = "Relative Risk (RR) with 95% Confidence Interval",
      y = "",
      color = "Estimate",
      shape = "95% CI"
    ) +
    theme(
      legend.position = "top",
      legend.box = "vertical", 
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      plot.title.position = "plot"
    )
  
  pdf_filename <- paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "pooled_effects/National_Pooled_Alert_DiD_",
    titel_anzeige,
    "_Final_English.pdf"
  )
  ggsave(pdf_filename, plot = p_forest, width = 11, height = 6)
  print(paste("--> Plot erfolgreich gespeichert als:", pdf_filename))
}

# ---------------------------------------------------------
# POOLING: ALERT VS. DID EFFEKT FÜR MODEL 0, 1 UND 2
# NUR FÜR LANDKREISE (AGS) MIT SIGNIFIKANTEM HITZEEFFEKT
# ---------------------------------------------------------

# 1a. DATEN EINLESEN 

df_all_cause <- read_excel("data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx") 
df_diagnosen <- read_excel("data/processed/fdz_outputs/2026-05-19_Results_FirstStage_Diagnosen_gh.xlsx")

df_combined <- bind_rows(df_all_cause, df_diagnosen)

# 1b. SIGNIFIKANTE AGS EINLESEN UND FILTERN 

df_sig_ags <- read.csv(
  paste0(
    "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/",
    "Signifikante_Kreise_mit_Fallzahlen_99Perzentil_def.csv"
  ),
  stringsAsFactors = FALSE
)

signifikante_ags_liste <- as.numeric(df_sig_ags$AGS) 

df_combined_filtered <- df_combined %>%
  mutate(
    AGS_Clean = as.numeric(gsub("AGS_|\\.csv", "", AGS_File))
  ) %>%
  filter(AGS_Clean %in% signifikante_ags_liste) 

print(paste("Anzahl Zeilen VOR Filter:", nrow(df_combined)))
print(paste("Anzahl Zeilen NACH Filter:", nrow(df_combined_filtered)))

modelle_auswahl <- c("model0", "model1", "model2")
diagnosen_liste <- c("All-Cause", "Infectious and Parasitic Diseases", "Respiratory Diseases", "Urinary Diseases")

# 2. SCHLEIFE ÜBER DIE MODELLE

for (aktuelles_modell in modelle_auswahl) {
  
  print(paste("=================================================="))
  print(paste("Start Meta-Analyse und Plot für:", aktuelles_modell))
  
  # A) Daten filtern für das aktuelle Modell
  df_clean <- df_combined_filtered %>%
    filter(Model == aktuelles_modell, RFM == "heat_alerts_rfm5") %>% 
    mutate(
      Diagnosis = case_when(
        Diagnosis == "AllAdmissions" ~ "All-Cause",
        Diagnosis == "Infect_par" ~ "Infectious and Parasitic Diseases",
        Diagnosis == "Resp_0" ~ "Respiratory Diseases",
        Diagnosis == "Uri_0" ~ "Urinary Diseases",
        TRUE ~ Diagnosis
      ),
      RR_Alert = as.numeric(RR_Alert),
      CI_Low_Alert = as.numeric(CI_Low_Alert),
      CI_Up_Alert = as.numeric(CI_Up_Alert),
      
      log_RR_Alert = log(RR_Alert),
      SE_Alert = (log(CI_Up_Alert) - log(CI_Low_Alert)) / (2 * 1.96),
      Var_Alert = SE_Alert^2,
      
      DiD_Coef = as.numeric(DiD_Coef),
      DiD_Variance = as.numeric(DiD_Variance)
    )
  
  # B) Meta-Analyse für dieses Modell
  results_list <- list()
  
  for (diag in diagnosen_liste) {
    
    df_sub <- df_clean %>% 
      filter(Diagnosis == diag) %>%
      filter(
        !is.na(log_RR_Alert) & !is.infinite(log_RR_Alert),
        !is.na(Var_Alert) & !is.infinite(Var_Alert) & Var_Alert > 0,
        !is.na(DiD_Coef) & !is.infinite(DiD_Coef),
        !is.na(DiD_Variance) & !is.infinite(DiD_Variance) & DiD_Variance > 0
      )
    
    if(nrow(df_sub) < 5) {
      print(paste("   Überspringe", diag, "- zu wenig gültige Daten!"))
      next
    }
    
    # Pooling ALERT
    meta_alert <- mixmeta(log_RR_Alert ~ 1, S = Var_Alert, data = df_sub, method = "reml")
    alert_coef <- coef(meta_alert)[1]
    alert_se <- sqrt(vcov(meta_alert)[1,1])
    
    # Pooling DID
    meta_did <- mixmeta(DiD_Coef ~ 1, S = DiD_Variance, data = df_sub, method = "reml")
    did_coef <- coef(meta_did)[1]
    did_se <- sqrt(vcov(meta_did)[1,1])
    
    results_list[[paste0(diag, "_Alert")]] <- data.frame(
      Diagnosis = diag,
      Effect_Type = "Alert", 
      RR = exp(alert_coef),
      CI_low = exp(alert_coef - 1.96 * alert_se),
      CI_high = exp(alert_coef + 1.96 * alert_se)
    )
    
    results_list[[paste0(diag, "_DiD")]] <- data.frame(
      Diagnosis = diag,
      Effect_Type = "DiD",
      RR = exp(did_coef),
      CI_low = exp(did_coef - 1.96 * did_se),
      CI_high = exp(did_coef + 1.96 * did_se)
    )
  }
  
  plot_data <- bind_rows(results_list)
  
  # C) Forest Plot für dieses Modell zeichnen
  plot_data <- plot_data %>%
    mutate(
      Precision_Status = case_when(
        is.na(CI_low) | is.na(CI_high) ~ "Unknown", 
        CI_low > 1 | CI_high < 1 ~ "95% CI excludes 1",
        TRUE ~ "95% CI includes 1"
      ),
      Precision_Status = factor(Precision_Status, levels = c("95% CI excludes 1", "95% CI includes 1", "Unknown"))
    )
  
  plot_data$Diagnosis <- factor(plot_data$Diagnosis, levels = rev(diagnosen_liste))
  plot_data$Effect_Type <- factor(plot_data$Effect_Type, levels = c("Alert", "DiD"))
  
  titel_anzeige <- toupper(aktuelles_modell) 
  
  subtitle_anzeige <- paste("Pooled estimates from", titel_anzeige, "| Subgroup: Counties with specific heat effect")
  
  p_forest <- ggplot(plot_data, aes(x = RR, y = Diagnosis, color = Effect_Type, group = Effect_Type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.3, linewidth = 1, position = position_dodge(width = 0.5)) +
    
    geom_point(aes(shape = Precision_Status), size = 4, stroke = 1.2, position = position_dodge(width = 0.5)) + 
    
    scale_color_manual(
      values = c("Alert" = "#1f78b4", "DiD" = "#e31a1c"),
      labels = c("Alert" = "Warning Effect (Pre-Implementation)", "DiD" = "Added System Effect (DiD)")
    ) +
    
    scale_shape_manual(values = c("95% CI excludes 1" = 18, "95% CI includes 1" = 5, "Unknown" = 4)) +
    
    theme_bw() +
    labs(
      title = "Heat Warning Effects by Diagnosis (Vulnerable AGS)",
      subtitle = subtitle_anzeige,
      x = "Relative Risk (RR) with 95% Confidence Interval",
      y = "",
      color = "Model Term",
      shape = "Precision"
    ) +
    theme(
      legend.position = "top",
      legend.box = "vertical", 
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    )
  
  pdf_filename <- paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "subgroup_analyses/SigAGS_Pooled_Alert_DiD_",
    titel_anzeige,
    "_Final_English.pdf"
  )
  ggsave(pdf_filename, plot = p_forest, width = 10, height = 6)
  print(paste("--> Plot erfolgreich gespeichert als:", pdf_filename))
}

print("==================================================")
print("Alle Plots für die Subgruppe wurden erstellt.")

# ---------------------------------------------------------
# ALL-IN-ONE SCRIPT: Figure Forest Plot & Figure Spatial Map
# Model 1 | All-Cause Admissions | National & State-Level
# ---------------------------------------------------------

print("Starte All-In-One Pipeline...")

# 1. PFADE DEFINIEREN

file_cases <- "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
path_shape <- "data/raw/spatial/vg2500/vg2500.shp"

# 2. DATEN EINLESEN UND BEREINIGEN 

print("Lese FDZ-Daten ein und bereite AGS/Bundesländer vor...")

# Bundesländer auf ENGLISCH
bundesland_dict <- c(
  "01" = "Schleswig-Holstein", "02" = "Hamburg", "03" = "Lower Saxony",
  "04" = "Bremen", "05" = "North Rhine-Westphalia", "06" = "Hesse",
  "07" = "Rhineland-Palatinate", "08" = "Baden-Wuerttemberg", "09" = "Bavaria",
  "10" = "Saarland", "11" = "Berlin", "12" = "Brandenburg",
  "13" = "Mecklenburg-Western Pomerania", "14" = "Saxony",
  "15" = "Saxony-Anhalt", "16" = "Thuringia"
)

df_clean <- read_excel(file_cases) %>%
  filter(Model == "model1", RFM == "heat_alerts_rfm5", Diagnosis == "AllAdmissions") %>%
  mutate(
    AGS_Clean = sprintf("%05d", as.numeric(gsub("AGS_|\\.csv", "", AGS_File))),
    BL_Code = substr(AGS_Clean, 1, 2),
    Bundesland = bundesland_dict[BL_Code],
    DiD_Coef = as.numeric(DiD_Coef),
    DiD_Variance = as.numeric(DiD_Variance)
  ) %>%
  filter(!is.na(DiD_Coef) & !is.na(DiD_Variance) & DiD_Variance > 0)

# 3. META-ANALYSE FÜR DEN FOREST PLOT

print("Berechne gepoolte Effekte pro Bundesland und für ganz Deutschland...")

results_list <- list()
states <- unique(df_clean$Bundesland)

for (state in states) {
  df_state <- df_clean %>% filter(Bundesland == state)
  
  if (nrow(df_state) > 1) {
    meta_state <- mixmeta(DiD_Coef ~ 1, S = DiD_Variance, data = df_state, method = "reml")
    estimate <- coef(meta_state)[1]
    standard_error <- sqrt(vcov(meta_state)[1,1])
  } else if (nrow(df_state) == 1) {
    estimate <- df_state$DiD_Coef[1]
    standard_error <- sqrt(df_state$DiD_Variance[1])
  } else {
    next
  }
  
  results_list[[state]] <- data.frame(
    Bundesland = state,
    Type = "State",
    RR = exp(estimate),
    CI_low = exp(estimate - 1.96 * standard_error),
    CI_high = exp(estimate + 1.96 * standard_error)
  )
}

meta_overall <- mixmeta(DiD_Coef ~ 1, S = DiD_Variance, data = df_clean, method = "reml")
results_list[["Overall"]] <- data.frame(
  Bundesland = "Overall",
  Type = "Overall",
  RR = exp(coef(meta_overall)[1]),
  CI_low = exp(coef(meta_overall)[1] - 1.96 * sqrt(vcov(meta_overall)[1,1])),
  CI_high = exp(coef(meta_overall)[1] + 1.96 * sqrt(vcov(meta_overall)[1,1]))
)

df_forest <- bind_rows(results_list)

df_forest$Bundesland <- factor(df_forest$Bundesland, 
                               levels = c("Overall", sort(states, decreasing = TRUE)))

# 4. FOREST PLOT ZEICHNEN (FIGURE 2A)

print("Erstelle Figure 2A (Forest Plot)...")

p_forest <- ggplot(df_forest, aes(x = RR, y = Bundesland, color = Type, shape = Type, size = Type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 0.8) +
  geom_point() +
  scale_color_manual(values = c("State" = "steelblue", "Overall" = "firebrick")) +
  scale_shape_manual(values = c("State" = 16, "Overall" = 18)) +
  scale_size_manual(values = c("State" = 3, "Overall" = 4.5)) +
  labs(
    title = "Regional Heterogeneity of the DiD Estimates",
    subtitle = "State-level Difference-in-Differences estimates for all-cause admissions",
    x = "Relative Risk (RR) with 95% Confidence Interval",
    y = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", margin = margin(b = 15)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
    axis.text.y = element_text(face = ifelse(levels(df_forest$Bundesland) == "Overall", "bold", "plain"), size = 11)
  )

ggsave(
  paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "regional_heterogeneity/Figure2A_ForestPlot_States_EN.pdf"
  ),
  plot = p_forest,
  width = 10,
  height = 7
)

# 5. SPATIAL MAP VORBEREITEN & ZEICHNEN (FIGURE 2B)

print("Lese Shapefile ein und erstelle Figure 2B (Spatial Map mit RR)...")

# 6. Shapefile einlesen

map_de <- read_sf(path_shape)

# Check
if("RS" %in% names(map_de)) {
  map_de <- map_de %>% rename(AGS = RS)
}

# 7. Daten mergen und RR vorbereiten
map_sf_did <- map_de %>%
  left_join(df_clean, by = c("AGS" = "AGS_Clean")) %>%
  mutate(
    RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
    CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD))),
    
    RR_DiD_Clipped = scales::oob_squish(RR_DiD, range = c(0.8, 1.2)),
    
    Is_Sig = ifelse(CI_Low_DiD > 1 | CI_Up_DiD < 1, TRUE, FALSE)
  )

# 8. Centroids (Mittelpunkte) NUR für die signifikanten Landkreise berechnen

sig_centroids <- map_sf_did %>% 
  filter(Is_Sig == TRUE) %>% 
  st_centroid()

# 9. Plot bauen

p_map <- ggplot() +
  geom_sf(data = map_sf_did, aes(fill = RR_DiD_Clipped), color = "white", linewidth = 0.1) +
  geom_sf(data = sig_centroids, color = "black", size = 0.6, shape = 20) +
  scale_fill_gradient2(
    low = "dodgerblue4", mid = "white", high = "firebrick", 
    midpoint = 1.0, 
    limits = c(0.8, 1.2), 
    name = "Relative Risk (RR)", 
    breaks = c(0.8, 0.9, 1.0, 1.1, 1.2),
    labels = c("<= 0.8", "0.9", "1.0", "1.1", ">= 1.2") 
  ) +
  labs(
    title = "Spatial Distribution of District-Specific DiD Estimates",
    subtitle = "District-specific relative risk for all-cause admissions (2000-2009)",
    caption = "Note: Black dots indicate statistical significance (95% CI does not cross 1.0).\nColor scale is truncated at RR values of 0.8 and 1.2 for enhanced contrast.\nBlue indicates a downward shift in the alert-day association, red indicates an upward shift."
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0, color = "darkgray", size = 9, margin = margin(t = 10))
  )

ggsave(
  paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "regional_heterogeneity/Figure2B_SpatialMap_RR_Sig_EN.pdf"
  ),
  plot = p_map,
  width = 8,
  height = 9,
  bg = "white"
)

# ---------------------------------------------------------
# ALL-IN-ONE PIPELINE: PRE vs. POST (GETRENNTE EXCEL-DATEIEN)
# Methode: Pre aus Standard-Modell, Post aus Reference-Switch-Modell
# ---------------------------------------------------------

rm(list=ls()) 

# 1. Pfade definieren 

file_standard <- "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
file_switch   <- "data/processed/fdz_outputs/Results_FirstStage_ReferenceSwitch_gh.xlsx"

bundesland_dict_de <- c(
  "01" = "Schleswig-Holstein", "02" = "Hamburg", "03" = "Niedersachsen",
  "04" = "Bremen", "05" = "Nordrhein-Westfalen", "06" = "Hessen",
  "07" = "Rheinland-Pfalz", "08" = "Baden-Württemberg", "09" = "Bayern",
  "10" = "Saarland", "11" = "Berlin", "12" = "Brandenburg",
  "13" = "Mecklenburg-Vorpommern", "14" = "Sachsen",
  "15" = "Sachsen-Anhalt", "16" = "Thüringen"
)

# 2a. PRE-EFFEKTE (aus der Standard-Analyse)

df_pre <- read_excel(file_standard) %>%
  filter(Model == "model1", RFM == "heat_alerts_rfm5", Diagnosis == "AllAdmissions") %>%
  mutate(
    AGS_Clean = sprintf("%05d", as.numeric(gsub("AGS_|\\.csv", "", AGS_File))),
    BL_Code = substr(AGS_Clean, 1, 2),
    Region_DE = bundesland_dict_de[BL_Code],
    
    # Pre-Effekt: RR ableiten
    RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    
    Pre_Coef = log(RR_Alert),
    Pre_Variance = ((log(CI_Up_Alert) - log(RR_Alert)) / 1.96)^2
  ) %>%
  select(AGS_Clean, Region_DE, Pre_Coef, Pre_Variance) %>%
  filter(!is.na(Pre_Coef) & Pre_Variance > 0)

# 2b. POST-EFFEKTE (aus der Reference-Switch-Analyse)

df_post <- read_excel(file_switch) %>%
  filter(Model == "model1", RFM == "heat_alerts_rfm5", Diagnosis == "AllAdmissions") %>%
  mutate(
    AGS_Clean = sprintf("%05d", as.numeric(gsub("AGS_|\\.csv", "", AGS_File))),
    
    # Post-Effekt: Im Switch ist Alert_Coef die Post-Periode
    Post_Coef = as.numeric(gsub(",", ".", as.character(Alert_Coef))),
    Post_Variance = as.numeric(gsub(",", ".", as.character(Alert_Variance)))
  ) %>%
  select(AGS_Clean, Post_Coef, Post_Variance) %>%
  filter(!is.na(Post_Coef) & Post_Variance > 0)

# 2c. MERGEN

df_merged <- inner_join(df_pre, df_post, by = "AGS_Clean")

results_list <- list()
states <- unique(df_merged$Region_DE)

for(state in states) {
  df_state <- df_merged %>% filter(Region_DE == state)
  
  # --- PRE-EFFEKT POOLEN (Standard) ---
  
  if(nrow(df_state) > 1) {
    meta_pre <- mixmeta(Pre_Coef ~ 1, S = Pre_Variance, data = df_state, method = "reml")
    pre_est <- coef(meta_pre)[1]
    pre_se  <- sqrt(vcov(meta_pre)[1,1])
  } else {
    pre_est <- df_state$Pre_Coef[1]
    pre_se  <- sqrt(df_state$Pre_Variance[1])
  }
  
  results_list[[paste0(state, "_Pre")]] <- data.frame(
    Region = state,
    Period = "Pre-Intervention (2000-2004)",
    RR = exp(pre_est),
    Low = exp(pre_est - 1.96 * pre_se),
    Up = exp(pre_est + 1.96 * pre_se)
  )
  
  # --- POST-EFFEKT POOLEN (Switch) ---
  
  if(nrow(df_state) > 1) {
    meta_post <- mixmeta(Post_Coef ~ 1, S = Post_Variance, data = df_state, method = "reml")
    post_est <- coef(meta_post)[1]
    post_se  <- sqrt(vcov(meta_post)[1,1])
  } else {
    post_est <- df_state$Post_Coef[1]
    post_se  <- sqrt(df_state$Post_Variance[1])
  }
  
  results_list[[paste0(state, "_Post")]] <- data.frame(
    Region = state,
    Period = "Post-Intervention (2005-2009)",
    RR = exp(post_est),
    Low = exp(post_est - 1.96 * post_se),
    Up = exp(post_est + 1.96 * post_se)
  )
}

df_plot_bl <- bind_rows(results_list) %>%
  mutate(
    Significance = ifelse(Low > 1 | Up < 1, "Significant", "Non-Significant")
  )

# 3. Übersetzung und absolut bombensichere Facet-Reihenfolge

df_plot_bl_en <- df_plot_bl %>%
  mutate(Region = as.character(Region)) %>%
  mutate(Region = recode(Region,
                         "Baden-Württemberg" = "Baden-Wuerttemberg",
                         "Bayern" = "Bavaria",
                         "Brandenburg" = "Brandenburg",
                         "Hessen" = "Hesse",
                         "Mecklenburg-Vorpommern" = "Mecklenburg-Western Pomerania",
                         "Niedersachsen" = "Lower Saxony",
                         "Nordrhein-Westfalen" = "North Rhine-Westphalia",
                         "Rheinland-Pfalz" = "Rhineland-Palatinate",
                         "Sachsen" = "Saxony",
                         "Sachsen-Anhalt" = "Saxony-Anhalt",
                         "Schleswig-Holstein" = "Schleswig-Holstein",
                         "Thüringen" = "Thuringia"
  )) %>%
  mutate(
    # Regionen alphabetisch sortieren (von unten nach oben im Plot)
    Region = factor(Region, levels = rev(sort(unique(Region)))),
    
    # Den kompletten Text (inkl. Zeilenumbruch \n) direkt in die Spalte schreiben
    Period = case_when(
      grepl("Pre", Period) ~ "Pre-Implementation Period (2000–2004):\nReconstructed Alert-Qualifying Days",
      grepl("Post", Period) ~ "Post-Implementation Period (2005–2009):\nOfficial DWD Alert Days"
    ),
    
    # Genau diese neuen Texte als feste Faktoren sortieren (Pre muss als erstes stehen)
    Period = factor(Period, levels = c(
      "Pre-Implementation Period (2000–2004):\nReconstructed Alert-Qualifying Days",
      "Post-Implementation Period (2005–2009):\nOfficial DWD Alert Days"
    ))
  )

# 4. Plot

p_forest_compare <- ggplot(df_plot_bl_en, aes(x = RR, y = Region, color = Significance)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Low, xmax = Up), height = 0.3, linewidth = 0.8) +
  geom_point(size = 3.5) +
  
  facet_wrap(~ Period) + 
  
  scale_color_manual(values = c("Significant" = "dodgerblue4", "Non-Significant" = "grey60")) +
  labs(
    title = "Alert-Day Associations Before and After HHWS Implementation by Federal State",
    subtitle = "Relative risks comparing alert with non-alert days within each implementation period",
    x = "Relative Risk (RR) with 95% CI",
    y = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "white"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
    axis.text.y = element_text(color = "black"),
    
    # Zentrierung
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    
    plot.title.position = "plot", 
    plot.margin = margin(t = 15, r = 15, b = 15, l = 15) 
  )

ggsave(
  paste0(
    "outputs/06_fdz_output_processing/06_08_additional_figures/",
    "temporal_comparisons/Figure_ForestPlot_Alert_PrePost_EN_Final.pdf"
  ),
  plot = p_forest_compare,
  width = 12,
  height = 6
)
