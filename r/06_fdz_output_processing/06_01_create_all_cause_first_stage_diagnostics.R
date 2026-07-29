# ==============================================================================
# Script: 06_01_create_all_cause_first_stage_diagnostics.R
#
# Purpose:
#   Processes district-specific all-cause first-stage estimates and creates
#   diagnostic visualizations for alert-day associations,
#   Difference-in-Differences estimates, and district-specific temperature-
#   response functions.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Environmental and heat-alert master dataset
#
# Outputs:
#   - Multi-model forest plots for selected and all districts
#   - Forest plots grouped by federal state
#   - District-specific time-invariant temperature-response plots
#   - District-specific time-varying temperature-response plots
#   - CSV file containing district-specific Model 1 estimates
#
# Analysis scope:
#   All-cause hospital admissions and Random Forest alert definition RFM5
#
# Required packages:
#   dplyr, ggplot2, readxl, dlnm, lubridate
#
# Dependencies:
#   Requires all-cause first-stage estimates generated within the secure FDZ
#   analysis environment and the environmental master dataset.
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

library(ggplot2)
library(dplyr)
library(dlnm)
library(lubridate) 
library(readxl) 

# =====================================================================
# 1. DATEN VORBEREITEN
# =====================================================================

# FDZ Output einlesen
df <- read_excel("data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx") %>%
  
  # Nur RFM5 
  filter(RFM == "heat_alerts_rfm5") %>%
  
  # SICHERHEITSCHECK: Kommas in Punkte umwandeln
  mutate(
    RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
    CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
    CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD)))
  ) %>%
  
  # AGS Zahlen extrahieren für saubere Namen und Bundesländer
  mutate(
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = gsub(".csv", "", AGS_Zahlen, fixed = TRUE),
    BL_Code = substr(AGS_Zahlen, 1, 2),
    Bundesland = case_when(
      BL_Code == "01" ~ "Schleswig-Holstein", BL_Code == "02" ~ "Hamburg",
      BL_Code == "03" ~ "Niedersachsen", BL_Code == "04" ~ "Bremen",
      BL_Code == "05" ~ "Nordrhein-Westfalen", BL_Code == "06" ~ "Hessen",
      BL_Code == "07" ~ "Rheinland-Pfalz", BL_Code == "08" ~ "Baden-Württemberg",
      BL_Code == "09" ~ "Bayern", BL_Code == "10" ~ "Saarland",
      BL_Code == "11" ~ "Berlin", BL_Code == "12" ~ "Brandenburg",
      BL_Code == "13" ~ "Mecklenburg-Vorpommern", BL_Code == "14" ~ "Sachsen",
      BL_Code == "15" ~ "Sachsen-Anhalt", BL_Code == "16" ~ "Thüringen",
      TRUE ~ "Unbekannt"
    )
  )

df_plot <- df %>% 
  filter(Model %in% c("model0", "model1", "model2")) %>%
  mutate(
    Model_Name = case_when(
      Model == "model0" ~ "Model0 (no temperature dlnm)",
      Model == "model1" ~ "Model1 (time-invariant temperature dlnm)",
      Model == "model2" ~ "Model2 (time-variant temperature dlnm)"
    ),
    # Legenden-Reihenfolge 
    Model_Name = factor(Model_Name, levels = c(
      "Model0 (no temperature dlnm)",
      "Model1 (time-invariant temperature dlnm)",
      "Model2 (time-variant temperature dlnm)"
    ))
  )

# =====================================================================
# 2. FUNKTION: MULTI-MODEL-FOREST-PLOT BAUEN
# =====================================================================

plot_forest_multi <- function(data, title, variable) {
  
  if(variable == "Alert") {
    data$x_val <- data$RR_Alert; data$ci_l <- data$CI_Low_Alert; data$ci_u <- data$CI_Up_Alert
    x_label <- "RR (alerts)"   # 1) HEAT ALERT EFFECT- Wie stark stiegen die EW an HW-Tagen in den Jahren 2000 bis 2004
    # (RRalert.pdf; RR_Alert in FDZ Output)
  } else {
    data$x_val <- data$RR_DiD; data$ci_l <- data$CI_Low_DiD; data$ci_u <- data$CI_Up_DiD
    x_label <- "RR (alerts*implementation)" # 2) INTERACTION ALERTS*IMPLEMENTATION (RR)- Veränderung des Grund-Risikos nach Einführung der echten Warnungen 
    # DiD; RRalertint; RR_DiD in FDZ Output
  }
  
  # Sortierung der y-Achse basierend auf Modell 2 
  sort_data <- data %>% filter(Model == "model2") %>% arrange(x_val)
  data$AGS_Zahlen <- factor(data$AGS_Zahlen, levels = sort_data$AGS_Zahlen)
  
  # 3 Balken pro AGS sauber übereinander
  pd <- position_dodge(width = 0.6)
  
  ggplot(data, aes(x = x_val, y = AGS_Zahlen, color = Model_Name)) +
    
    # Null-Linie 
    geom_vline(xintercept = 1, linetype = "solid", color = "black", linewidth = 0.5) +
    
    # Fehlerbalken und Punkte mit Dodging (Verschiebung)
    geom_errorbarh(aes(xmin = ci_l, xmax = ci_u), height = 0.3, position = pd, linewidth = 0.8) +
    geom_point(position = pd, size = 2.5) +
    
    scale_color_manual(values = c(
      "Model0 (no temperature dlnm)" = "blue",
      "Model1 (time-invariant temperature dlnm)" = "black",
      "Model2 (time-variant temperature dlnm)" = "orange"
    )) +
    
    theme_bw() +
    labs(title = title, x = x_label, y = "") +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "top", # Legende nach oben
      legend.title = element_blank(), # Legenden-Titel ausblenden
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey80") # Gepunktete Führungslinien
    )
}

# =====================================================================
# 3. PDF EXPORT
# =====================================================================

# 1. TOP 10 LANDKREISE

top10_ags <- df_plot %>% filter(Model == "model2") %>% arrange(desc(Total_Cases)) %>% head(10) %>% pull(AGS_Zahlen)
df_top10 <- df_plot %>% filter(AGS_Zahlen %in% top10_ags)

pdf(
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/Check_Top10_MultiModel_def.pdf",
  width = 8,
  height = 6
)

print(plot_forest_multi(df_top10, "Top 10 Landkreise (Alert Effekt)", "Alert"))
print(plot_forest_multi(df_top10, "Top 10 Landkreise (DiD Effekt)", "DiD"))
dev.off()

# 2. ALLE 400 AGS 

ags_list <- unique(df_plot$AGS_Zahlen)
chunk_size <- 25
total_chunks <- ceiling(length(ags_list) / chunk_size)

pdf(
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/Check_Alle_400_AGS_Alert_MultiModel_def.pdf",
  width = 8,
  height = 9
)

for(i in 1:total_chunks) {
  current_ags <- ags_list[((i-1)*chunk_size + 1):min(i*chunk_size, length(ags_list))]
  df_chunk <- df_plot %>% filter(AGS_Zahlen %in% current_ags)
  
  p <- plot_forest_multi(df_chunk, paste("Alle AGS (Alert) - Teil", i, "von", total_chunks), "Alert")
  print(p)
}
dev.off()

pdf(
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/Check_Alle_400_AGS_DiD_MultiModel_def.pdf",
  width = 8,
  height = 9
)

for(i in 1:total_chunks) {
  current_ags <- ags_list[((i-1)*chunk_size + 1):min(i*chunk_size, length(ags_list))]
  df_chunk <- df_plot %>% filter(AGS_Zahlen %in% current_ags)
  
  p <- plot_forest_multi(df_chunk, paste("Alle AGS (DiD) - Teil", i, "von", total_chunks), "DiD")
  print(p)
}
dev.off()

# =====================================================================
# 4. OUTPUT 3: ALLE AGS GRUPPIERT NACH BUNDESLAND (Multi-Model)
# =====================================================================

bundeslaender <- sort(unique(df_plot$Bundesland))

pdf(
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/Check_Alle_AGS_nach_Bundesland_MultiModel_def.pdf",
  width = 8,
  height = 12
)

for(bl in bundeslaender) {
  
  # 1. Finde alle Landkreise, die zu diesem Bundesland gehören
  bl_ags <- df %>% 
    filter(Bundesland == bl) %>% 
    pull(AGS_Zahlen) %>% 
    unique()
  
  # 2. Filtere unseren Plot-Datensatz (der Modell 0, 1 und 2 enthält) auf genau diese Kreise
  df_bl <- df_plot %>% filter(AGS_Zahlen %in% bl_ags)
  
  # 3. Nur plotten, wenn auch wirklich Daten für das BL da sind
  if(nrow(df_bl) > 0) {
    # Alert Effekt plotten
    p1 <- plot_forest_multi(df_bl, paste(bl, "- Alert Effekt"), "Alert")
    print(p1)
    # DiD Effekt plotten
    p2 <- plot_forest_multi(df_bl, paste(bl, "- DiD Effekt"), "DiD")
    print(p2)
  }
}
dev.off()

#######################################################################
# PLOT DLNM CURVES (FIXED & TIME-VARIANT) AUS FDZ-EXPORTEN
#######################################################################

# =====================================================================
# 1. PARAMETER DEINES MODELLS (wie First-Stage Config-Skript)
# =====================================================================

varfun    <- "bs"              # Basis-Funktion Temp
vardegree <- 2                 # Grad
varper    <- c(50, 90)         # Knots (Perzentile)
lag       <- 3                 # Lag in Tagen (Verzögerung des Effekts)
dfseas    <- 4                 # Freiheitsgrade Saison
dftrend   <- 1                 # Freiheitsgrade Trend

# =====================================================================
# 2. DATEN EINLESEN
# =====================================================================

# FDZ Daten einlesen
df_fdz <- readxl::read_excel("data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx") %>%
  filter(RFM == "heat_alerts_rfm5") %>%
  mutate(
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = gsub(".csv", "", AGS_Zahlen, fixed = TRUE)
  )

# Temp-Daten einlesen und auf Sommer filtern
df_temp <- read.csv("data/external/fdz_input/Masterfile_Final_for_FDZ.csv") %>%
  mutate(
    AGS = sprintf("%05d", as.numeric(AGS)), 
    Datum = as.Date(Datum, format="%Y-%m-%d"),
    Monat = month(Datum) 
  ) %>%
  # nur Mai (5) bis September (9)
  filter(Monat >= 5 & Monat <= 9) 

ags_list <- unique(df_fdz$AGS_Zahlen)

# =====================================================================
# 3. HILFSFUNKTIONEN ZUM ENTPACKEN DER FDZ-MATRIZEN
# =====================================================================

get_coef <- function(data, mod) {
  val <- data %>% filter(Model == mod) %>% pull(Temp_Coef_Vector)
  if(length(val) == 0 || is.na(val)) return(NULL)
  as.numeric(unlist(strsplit(as.character(val), "\\|")))
}

get_vcov <- function(data, mod) {
  val <- data %>% filter(Model == mod) %>% pull(Temp_Vcov_Matrix)
  if(length(val) == 0 || is.na(val)) return(NULL)
  vec <- as.numeric(unlist(strsplit(as.character(val), "\\|")))
  n <- sqrt(length(vec)) 
  matrix(vec, nrow=n, ncol=n)
}

# =====================================================================
# 4. PLOT: TIME INVARIANT DLNM (Model 1)
# =====================================================================

pdf(
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/DlnmFixCheck_AGS_def.pdf",
  width = 12,
  height = 14
)
layout(matrix(1:15, ncol=3, byrow=T))
par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.2, 0.8, 0), las=1)

for (i in seq_along(ags_list)) {
  current_ags <- ags_list[i]
  
  fdz_sub <- df_fdz %>% filter(AGS_Zahlen == current_ags)
  citydat <- df_temp %>% filter(AGS == current_ags)
  
  if(nrow(fdz_sub) == 0 || nrow(citydat) == 0) next
  
  coef1 <- get_coef(fdz_sub, "model1")
  vcov1 <- get_vcov(fdz_sub, "model1")
  
  if(is.null(coef1) || is.null(vcov1)) next
  
  # varfun, vardegree und varper, um die Kurve ("bvar") zu bauen
  argvar <- list(x = citydat$T_mean_popw, fun = varfun, degree = vardegree, 
                 knots = quantile(citydat$T_mean_popw, varper/100, na.rm=T))
  bvar <- do.call(onebasis, argvar)
  
  mean_temp <- mean(citydat$T_mean_popw, na.rm=T)
  
  cp <- crosspred(bvar, coef=coef1, vcov=vcov1, model.link="log", 
                  by=0.1, cen=mean_temp)
  
  plot(cp, ylim=c(0.9, max(cp$allRRhigh, na.rm=T)+0.1),
       xlim=c(min(citydat$T_mean_popw, na.rm=T)-1, max(citydat$T_mean_popw, na.rm=T)+1),
       xlab="Temperature", ylab="RR", lab=c(6,5,7), lwd=1, main=current_ags, col="blue")
  
  abline(v=mean_temp, lty=3, lwd=0.9, col=grey(0.8))
}
dev.off()

print("Time Invariant Plot (Model 1) erstellt!")

# =====================================================================
# 5. PLOT: TIME-VARIABLE DLNM (Model 2 vs Model 3)
# =====================================================================

pdf(
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/DlnmVarCheck_AGS_def.pdf",
  width = 12,
  height = 14
)
layout(matrix(1:15, ncol=3, byrow=T))
par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.2, 0.8, 0), las=1)

for (i in seq_along(ags_list)) {
  current_ags <- ags_list[i]
  
  fdz_sub <- df_fdz %>% filter(AGS_Zahlen == current_ags)
  citydat <- df_temp %>% filter(AGS == current_ags)
  
  if(nrow(fdz_sub) == 0 || nrow(citydat) == 0) next
  
  coef2 <- get_coef(fdz_sub, "model2"); vcov2 <- get_vcov(fdz_sub, "model2")
  coef3 <- get_coef(fdz_sub, "model3"); vcov3 <- get_vcov(fdz_sub, "model3")
  
  if(is.null(coef2) || is.null(coef3)) next
  
  argvar <- list(x = citydat$T_mean_popw, fun = varfun, degree = vardegree, 
                 knots = quantile(citydat$T_mean_popw, varper/100, na.rm=T))
  bvar <- do.call(onebasis, argvar)
  
  mean_temp <- mean(citydat$T_mean_popw, na.rm=T)
  
  cp2 <- crosspred(bvar, coef=coef2, vcov=vcov2, model.link="log", by=0.1, cen=mean_temp)
  cp3 <- crosspred(bvar, coef=coef3, vcov=vcov3, model.link="log", by=0.1, cen=mean_temp)
  
  plot(cp2, ylim=c(0.9, max(c(cp2$allRRhigh, cp3$allRRhigh), na.rm=T)+0.1),
       xlim=c(min(citydat$T_mean_popw, na.rm=T)-1, max(citydat$T_mean_popw, na.rm=T)+1),
       xlab="Temperature", ylab="RR", lab=c(6,5,7), lwd=1.5, main=current_ags, 
       col="blue", ci.arg=list(density=20, col="blue"))
  
  lines(cp3, lwd=1.5, col="black", ci="area", ci.arg=list(density=20, angle=-45, col="black"))
  
  abline(v=mean_temp, lty=3, lwd=0.9, col=grey(0.8))
  
  legend("topleft", c("Periode 1 (Potentiell)", "Periode 2 (Echt)"), 
         col=c("blue","black"), lwd=1.5, cex=0.8, inset=0.05, bty="n")
}
dev.off()

print("Time Variable Plot (Model 2 vs 3) erstellt.")

#######################################################################
# Zahlen für Paper extrahieren
#######################################################################

# =====================================================================
# EXTRA-SKRIPT: DIE WICHTIGSTEN ZAHLEN FÜR DEN TEXT EXTRAHIEREN
# =====================================================================

# 1. DATEN NEU EINLESEN 

df <- read_excel("data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx") %>%
  # Nur RFM5 
  filter(RFM == "heat_alerts_rfm5") %>%
  
  # SICHERHEITSCHECK: Kommas in Punkte umwandeln, damit R rechnen kann
  mutate(
    RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
    CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
    CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD)))
  ) %>%
  
  # AGS Zahlen extrahieren
  mutate(
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = gsub(".csv", "", AGS_Zahlen, fixed = TRUE),
    BL_Code = substr(AGS_Zahlen, 1, 2),
    Bundesland = case_when(
      BL_Code == "01" ~ "Schleswig-Holstein", BL_Code == "02" ~ "Hamburg",
      BL_Code == "03" ~ "Niedersachsen", BL_Code == "04" ~ "Bremen",
      BL_Code == "05" ~ "Nordrhein-Westfalen", BL_Code == "06" ~ "Hessen",
      BL_Code == "07" ~ "Rheinland-Pfalz", BL_Code == "08" ~ "Baden-Württemberg",
      BL_Code == "09" ~ "Bayern", BL_Code == "10" ~ "Saarland",
      BL_Code == "11" ~ "Berlin", BL_Code == "12" ~ "Brandenburg",
      BL_Code == "13" ~ "Mecklenburg-Vorpommern", BL_Code == "14" ~ "Sachsen",
      BL_Code == "15" ~ "Sachsen-Anhalt", BL_Code == "16" ~ "Thüringen",
      TRUE ~ "Unbekannt"
    )
  )

# 2. FINALES MODELL FILTERN (Model 1)

df_final <- df %>%
  filter(Model == "model1") %>%
  select(AGS_Zahlen, Bundesland, Total_Cases, 
         RR_Alert, CI_Low_Alert, CI_Up_Alert, 
         RR_DiD, CI_Low_DiD, CI_Up_DiD) %>%
  arrange(desc(Total_Cases))

# 3. ZUSAMMENFASSUNG BERECHNEN

summary_stats <- df_final %>%
  summarise(
    Median_RR_Alert = median(RR_Alert, na.rm = TRUE),
    Districts_Alert_Over_1 = sum(RR_Alert > 1, na.rm = TRUE),
    Districts_Alert_Under_1 = sum(RR_Alert < 1, na.rm = TRUE),
    
    Median_RR_DiD = median(RR_DiD, na.rm = TRUE),
    Districts_DiD_Over_1 = sum(RR_DiD > 1, na.rm = TRUE),
    Districts_DiD_Under_1 = sum(RR_DiD < 1, na.rm = TRUE)
  )

# 4. OUTPUT FÜR DEN TEXT ANZEIGEN & SPEICHERN

cat("\n=== SUMMARY STATS FÜR DEN FLIEßTEXT ===\n")
print(summary_stats)

cat("\n=== TOP 5 GRÖßTE LANDKREISE ===\n")
print(head(df_final, 5))

write.csv(
  df_final,
  "outputs/06_fdz_output_processing/06_01_all_cause_first_stage_diagnostics/FirstStage_ExactNumbers_Model1.csv",
  row.names = FALSE
)