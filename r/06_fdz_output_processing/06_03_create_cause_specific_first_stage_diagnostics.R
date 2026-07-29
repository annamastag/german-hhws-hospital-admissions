# ==============================================================================
# Script: 06_03_create_cause_specific_first_stage_diagnostics.R
#
# Purpose:
#   Combines all-cause and cause-specific first-stage estimates and creates
#   diagnosis-stratified diagnostic plots for alert-day associations,
#   Difference-in-Differences estimates, and temperature-response functions.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Released cause-specific first-stage FDZ output table including
#     cardiovascular admissions
#   - Environmental and heat-alert master dataset
#
# Outputs:
#   - Diagnosis-stratified multi-model forest plot for selected districts
#   - District-specific time-invariant temperature-response plots by diagnosis
#   - District-specific time-varying temperature-response plots by diagnosis
#
# Outcomes:
#   - All-cause admissions
#   - Cardiovascular admissions
#   - Infectious and parasitic disease admissions
#   - Respiratory admissions
#   - Urogenital admissions
#
# Required packages:
#   dplyr, ggplot2, readxl, dlnm, lubridate
#
# Dependencies:
#   Requires all-cause and cause-specific first-stage estimates generated in
#   the secure FDZ environment.
#
# Repository note:
#   Only the file name, path strings, and explanatory comments were standardized.
# ==============================================================================

rm(list=ls())

library(dplyr)
library(ggplot2)
library(readxl) 
library(dlnm)
library(lubridate) 

# =====================================================================
# 1. DATEN VORBEREITEN & ZUSAMMENFÜGEN (FRANKENSTEIN-MERGE)
# =====================================================================

# A) Das alte File für AllAdmissions laden
df_all_cause <- read_excel(
  "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
)

# B) Das NEUE File vom FDZ laden (enthält Resp_0, Uri_0, Infect_par UND Card_0)
df_diagnosen <- read_excel(
  "data/processed/fdz_outputs/02_Analysis_FirstStage_default_allCause_and_Strat(mit_Card)_und_sens_Res_gh.xlsx"
)

# C) Beide Tabellen verbinden
df_combined <- bind_rows(df_all_cause, df_diagnosen)

# D) Variablen bereinigen und übersetzen
df <- df_combined %>%
  filter(RFM == "heat_alerts_rfm5") %>%
  mutate(
    # Diagnosen in sprechende Namen übersetzen
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
    Model_Name = factor(Model_Name, levels = c(
      "Model0 (no temperature dlnm)",
      "Model1 (time-invariant temperature dlnm)",
      "Model2 (time-variant temperature dlnm)"
    )),
    # Wichtig: Die neue Reihenfolge der Facets festlegen
    Diagnosis = factor(Diagnosis, levels = c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital")) 
  )

print("Datenvorbereitung abgeschlossen. Es sind 5 Diagnosegruppen im Datensatz.")

# =====================================================================
# 2. FUNKTION: STRATIFIZIERTE FOREST-PLOTS (MIT FACET WRAP)
# =====================================================================

plot_forest_multi_strat <- function(data, title, variable) {
  
  if(variable == "Alert") {
    data$x_val <- data$RR_Alert; data$ci_l <- data$CI_Low_Alert; data$ci_u <- data$CI_Up_Alert
    x_label <- "RR (alerts)"   
  } else {
    data$x_val <- data$RR_DiD; data$ci_l <- data$CI_Low_DiD; data$ci_u <- data$CI_Up_DiD
    x_label <- "RR (alerts*implementation)" 
  }
  
  # Y-Achsen Sortierung (Referenz ist hier immer die erste Gruppe, also "All-Cause"): 
  ref_diag <- levels(data$Diagnosis)[1]
  sort_data <- data %>% filter(Model == "model2" & Diagnosis == ref_diag) %>% arrange(x_val)
  data$AGS_Zahlen <- factor(data$AGS_Zahlen, levels = sort_data$AGS_Zahlen)
  
  pd <- position_dodge(width = 0.6)
  
  x_limits <- c(min(data$ci_l, na.rm=TRUE)*0.95, min(max(data$ci_u, na.rm=TRUE)*1.05, 3.0))
  
  ggplot(data, aes(x = x_val, y = AGS_Zahlen, color = Model_Name)) +
    geom_vline(xintercept = 1, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_l, xmax = ci_u), height = 0.3, position = pd, linewidth = 0.8) +
    geom_point(position = pd, size = 2.5) +
    
    scale_color_manual(values = c(
      "Model0 (no temperature dlnm)" = "blue",
      "Model1 (time-invariant temperature dlnm)" = "black",
      "Model2 (time-variant temperature dlnm)" = "orange"
    )) +
    
    # Der Plot wird nach 5 Diagnosegruppen in Spalten zerlegt!
    facet_wrap(~ Diagnosis, ncol = length(levels(data$Diagnosis))) +
    
    coord_cartesian(xlim = x_limits) +
    theme_bw() +
    labs(title = title, x = x_label, y = "") +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "top", 
      legend.title = element_blank(), 
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

# =====================================================================
# 3. PDF EXPORT: FOREST PLOTS
# =====================================================================

ref_diag <- levels(df_plot$Diagnosis)[1]
top10_ags <- df_plot %>% filter(Model == "model2" & Diagnosis == ref_diag) %>% arrange(desc(Total_Cases)) %>% head(10) %>% pull(AGS_Zahlen)
df_top10 <- df_plot %>% filter(AGS_Zahlen %in% top10_ags)

print("Generiere stratifizierten Forest-Plot für Top 10 AGS...")
pdf(
  "outputs/06_fdz_output_processing/06_03_cause_specific_first_stage_diagnostics/Check_Top10_Stratified_MultiModel_inkl_Cardio.pdf",
  width = 16,
  height = 7
)
print(plot_forest_multi_strat(df_top10, "Top 10 Landkreise (Alert Effekt) - Nach Diagnose", "Alert"))
print(plot_forest_multi_strat(df_top10, "Top 10 Landkreise (DiD Effekt) - Nach Diagnose", "DiD"))
dev.off()

# =====================================================================
# 4. DLNM CURVES: PARAMETER & HILFSFUNKTIONEN
# =====================================================================

varfun    <- "bs"              
vardegree <- 2                 
varper    <- c(50, 90)         
lag       <- 3                 
dfseas    <- 4                 
dftrend   <- 1                 

df_temp <- read.csv(
  "data/external/fdz_input/Masterfile_Final_for_FDZ.csv",
  sep = ","
) %>%
  mutate(
    AGS = sprintf("%05d", as.numeric(AGS)), 
    Datum = as.Date(Datum, format="%Y-%m-%d"),
    Monat = month(Datum) 
  ) %>%
  filter(Monat >= 5 & Monat <= 9) 

get_coef <- function(data, mod) {
  val <- data %>% filter(Model == mod) %>% pull(Temp_Coef_Vector)
  if(length(val) == 0 || is.na(val[1])) return(NULL) 
  as.numeric(unlist(strsplit(as.character(val[1]), "\\|")))
}

get_vcov <- function(data, mod) {
  val <- data %>% filter(Model == mod) %>% pull(Temp_Vcov_Matrix)
  if(length(val) == 0 || is.na(val[1])) return(NULL)
  vec <- as.numeric(unlist(strsplit(as.character(val[1]), "\\|")))
  n <- sqrt(length(vec)) 
  matrix(vec, nrow=n, ncol=n)
}

# =====================================================================
# 5. DLNM CURVES ZEICHNEN (GETRENNT NACH DIAGNOSE)
# =====================================================================

ags_list <- unique(df$AGS_Zahlen)
diagnosen <- levels(df_plot$Diagnosis)

for (diag in diagnosen) {
  
  safe_diag_name <- gsub(" & | |&", "_", diag)
  
  pdf_name <- paste0(
    "outputs/06_fdz_output_processing/06_03_cause_specific_first_stage_diagnostics/",
    "DlnmFixCheck_",
    safe_diag_name,
    ".pdf"
  )
  pdf(pdf_name, width=12, height=14) 
  layout(matrix(1:15, ncol=3, byrow=T))
  par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.2, 0.8, 0), las=1)
  
  for (i in seq_along(ags_list)) {
    current_ags <- ags_list[i]
    
    fdz_sub <- df %>% filter(AGS_Zahlen == current_ags & Diagnosis == diag)
    citydat <- df_temp %>% filter(AGS == current_ags)
    
    if(nrow(fdz_sub) == 0 || nrow(citydat) == 0) next
    
    coef1 <- get_coef(fdz_sub, "model1")
    vcov1 <- get_vcov(fdz_sub, "model1")
    
    if(is.null(coef1) || is.null(vcov1)) next
    
    argvar <- list(x = citydat$T_mean_popw, fun = varfun, degree = vardegree, 
                   knots = quantile(citydat$T_mean_popw, varper/100, na.rm=T))
    bvar <- do.call(onebasis, argvar)
    
    mean_temp <- mean(citydat$T_mean_popw, na.rm=T)
    
    cp <- crosspred(bvar, coef=coef1, vcov=vcov1, model.link="log", by=0.1, cen=mean_temp)
    
    plot(cp, ylim=c(0.9, max(cp$allRRhigh, na.rm=T)+0.1),
         xlim=c(min(citydat$T_mean_popw, na.rm=T)-1, max(citydat$T_mean_popw, na.rm=T)+1),
         xlab="Temperature", ylab="RR", lab=c(6,5,7), lwd=1, 
         main=paste(current_ags, "-", diag), col="blue") 
    
    abline(v=mean_temp, lty=3, lwd=0.9, col=grey(0.8))
  }
  dev.off()
  print(paste("Time Invariant Plot für", diag, "erstellt!"))
}

# =====================================================================
# 6. DLNM VARIABLE CHECK (Model 2 vs Model 3) GETRENNT NACH DIAGNOSE
# =====================================================================

for (diag in diagnosen) {
  
  safe_diag_name <- gsub(" & | |&", "_", diag)
  
  pdf_name <- paste0(
    "outputs/06_fdz_output_processing/06_03_cause_specific_first_stage_diagnostics/",
    "DlnmVarCheck_",
    safe_diag_name,
    ".pdf"
  )
  pdf(pdf_name, width=12, height=14)
  layout(matrix(1:15, ncol=3, byrow=T))
  par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.2, 0.8, 0), las=1)
  
  for (i in seq_along(ags_list)) {
    current_ags <- ags_list[i]
    
    fdz_sub <- df %>% filter(AGS_Zahlen == current_ags & Diagnosis == diag)
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
         xlab="Temperature", ylab="RR", lab=c(6,5,7), lwd=1.5, 
         main=paste(current_ags, "-", diag), 
         col="blue", ci.arg=list(density=20, col="blue"))
    
    lines(cp3, lwd=1.5, col="black", ci="area", 
          ci.arg=list(density=20, angle=-45, col="black"))
    
    abline(v=mean_temp, lty=3, lwd=0.9, col=grey(0.8))
    
    legend("topleft", c("Period 1 (Pre-System)", "Period 2 (Post-System)"), 
           col=c("blue","black"), lwd=1.5, cex=0.7, inset=0.05, bty="n")
  }
  dev.off()
  print(paste("Time Variable Plot (Model 2 vs 3) für", diag, "erstellt!"))
}
