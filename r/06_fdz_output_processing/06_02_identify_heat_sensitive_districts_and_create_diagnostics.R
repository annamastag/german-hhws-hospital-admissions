# ==============================================================================
# Script: 06_02_identify_heat_sensitive_districts_and_create_diagnostics.R
#
# Purpose:
#   Identifies districts with statistically detectable temperature-admission
#   associations at the district-specific 99th temperature percentile and
#   creates exploratory diagnostic plots for the selected districts.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Environmental and heat-alert master dataset
#   - Harmonized German district boundary shapefile
#
# Outputs:
#   - CSV file listing selected districts and total admission counts
#   - Multi-model forest plots for selected districts
#   - District-specific temperature-response plots
#   - Choropleth map of statistically detectable Model 1 alert associations
#
# Analysis scope:
#   Exploratory, data-dependent subgroup analysis based on district-specific
#   Model 1 temperature associations at the 99th percentile
#
# Required packages:
#   dplyr, ggplot2, dlnm, lubridate, readxl, sf, viridis
#
# Dependencies:
#   Requires all-cause first-stage estimates, the environmental master dataset,
#   and harmonized district boundaries.
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(dlnm)
library(lubridate) 
library(readxl)   
library(sf)
library(viridis)

# =====================================================================
# 0. ZENTRALE PFADE DEFINIEREN
# =====================================================================

# Pfad zu FDZ Excel-Datei
file_first_stage <-
  "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"

# Pfad zu Temperatur-Daten
file_temp_data <-
  "data/external/fdz_input/Masterfile_Final_for_FDZ.csv"

# Pfad zum Shapefile für die Karte
shapefile_path <-
  "data/raw/spatial/vg2500/vg2500.shp"

# =====================================================================
# 1. DATEN EINLESEN UND AUFBEREITEN
# =====================================================================

df_fdz <- read_excel(file_first_stage) %>%
  
  # Nur RFM5 betrachten
  filter(RFM == "heat_alerts_rfm5") %>%
  
  # SICHERHEITSCHECK
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

ags_list <- unique(df_fdz$AGS_Zahlen)

# Temperaturdaten einlesen (nur Sommer)
df_temp <- read.csv(file_temp_data, sep=",") %>%
  mutate(
    AGS = sprintf("%05d", as.numeric(AGS)), 
    Datum = as.Date(Datum, format="%Y-%m-%d"),
    Monat = month(Datum) 
  ) %>%
  filter(Monat >= 5 & Monat <= 9) 


# Hilfsfunktionen für DLNM Matrizen
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
# 2. EXPLORATION: SIGNIFIKANTE HITZE-EFFEKTE (99. PERZENTIL)
# =====================================================================

varfun    <- "bs"              
vardegree <- 2                 
varper    <- c(50, 90)         

results_99 <- data.frame(
  AGS = character(), Temp_99 = numeric(), RR_99 = numeric(),
  CI_low = numeric(), CI_high = numeric(), Signifikanz = character(),
  stringsAsFactors = FALSE
)

for (current_ags in ags_list) {
  fdz_sub <- df_fdz %>% filter(AGS_Zahlen == current_ags)
  citydat <- df_temp %>% filter(AGS == current_ags)
  if(nrow(fdz_sub) == 0 || nrow(citydat) == 0) next
  
  cases_val <- as.numeric(fdz_sub$Total_Cases[1])
  coef1 <- get_coef(fdz_sub, "model1")
  vcov1 <- get_vcov(fdz_sub, "model1")
  if(is.null(coef1) || is.null(vcov1)) next
  
  temp_99_val <- quantile(citydat$T_mean_popw, 0.99, na.rm=TRUE)
  mean_temp <- mean(citydat$T_mean_popw, na.rm=TRUE)
  
  argvar <- list(x = citydat$T_mean_popw, fun = varfun, degree = vardegree, 
                 knots = quantile(citydat$T_mean_popw, varper/100, na.rm=T))
  bvar <- do.call(onebasis, argvar)
  
  cp <- crosspred(bvar, coef=coef1, vcov=vcov1, model.link="log", 
                  at=temp_99_val, cen=mean_temp)
  
  rr <- as.numeric(cp$allRRfit[1]); ci_l <- as.numeric(cp$allRRlow[1]); ci_u <- as.numeric(cp$allRRhigh[1])
  
  sig <- "Nicht signifikant"
  if(ci_l > 1) sig <- "Signifikant Positiv (RR > 1)"
  if(ci_u < 1) sig <- "Signifikant Negativ (RR < 1)"
  
  results_99 <- rbind(results_99, data.frame(
    AGS = current_ags, Temp_99 = round(temp_99_val, 1), RR_99 = round(rr, 3),
    CI_low = round(ci_l, 3), CI_high = round(ci_u, 3), Signifikanz = sig
  ))
}

# Extrahieren der Kreise, die jetzt signifikant sind!
sig_ags <- results_99 %>% filter(Signifikanz != "Nicht signifikant")
sig_ags_liste <- sig_ags$AGS

# Fallzahlen dazuspielen und sortieren
sig_ags_with_cases <- sig_ags %>%
  left_join(df_fdz %>% select(AGS_Zahlen, Total_Cases) %>% unique(), 
            by = c("AGS" = "AGS_Zahlen")) %>%
  arrange(desc(Total_Cases))

# Datei neu abspeichern (enthält alle 37 Kreise)
write.csv(
  sig_ags_with_cases,
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/Signifikante_Kreise_mit_Fallzahlen_99Perzentil_def.csv",
  row.names = FALSE
)

# =====================================================================
# 3. FOREST PLOTS (NUR FÜR DIE SIGNIFIKANTEN KREISE)
# =====================================================================

df_plot <- df_fdz %>% 
  filter(AGS_Zahlen %in% sig_ags_liste) %>%
  filter(Model %in% c("model0", "model1", "model2")) %>%
  mutate(
    Model_Name = case_when(
      Model == "model0" ~ "Model0 (no temperature dlnm)",
      Model == "model1" ~ "Model1 (time-invariant temperature dlnm)",
      Model == "model2" ~ "Model2 (time-variant temperature dlnm)"
    ),
    Model_Name = factor(Model_Name, levels = c(
      "Model0 (no temperature dlnm)", "Model1 (time-invariant temperature dlnm)", "Model2 (time-variant temperature dlnm)"
    ))
  )

plot_forest_multi <- function(data, title, variable) {
  if(variable == "Alert") {
    data$x_val <- data$RR_Alert; data$ci_l <- data$CI_Low_Alert; data$ci_u <- data$CI_Up_Alert
    x_label <- "RR (alerts)"
  } else {
    data$x_val <- data$RR_DiD; data$ci_l <- data$CI_Low_DiD; data$ci_u <- data$CI_Up_DiD
    x_label <- "RR (alerts*implementation)" 
  }
  sort_data <- data %>% filter(Model == "model2") %>% arrange(x_val)
  data$AGS_Zahlen <- factor(data$AGS_Zahlen, levels = sort_data$AGS_Zahlen)
  pd <- position_dodge(width = 0.6)
  
  ggplot(data, aes(x = x_val, y = AGS_Zahlen, color = Model_Name)) +
    geom_vline(xintercept = 1, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_l, xmax = ci_u), height = 0.3, position = pd, linewidth = 0.8) +
    geom_point(position = pd, size = 2.5) +
    scale_color_manual(values = c("blue", "black", "orange")) +
    theme_bw() + labs(title = title, x = x_label, y = "") +
    theme(axis.text.y = element_text(size = 7), legend.position = "top", legend.title = element_blank(), 
          legend.text = element_text(size = 9), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"))
}

# Top 10 Plots
top10_ags <- df_plot %>% filter(Model == "model2") %>% arrange(desc(Total_Cases)) %>% head(10) %>% pull(AGS_Zahlen)
df_top10 <- df_plot %>% filter(AGS_Zahlen %in% top10_ags)

pdf(
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/Check_Top10_SIGNIFIKANT_MultiModel.pdf",
  width = 8,
  height = 6
)
print(plot_forest_multi(df_top10, "Top 10 Signifikante Landkreise (Alert Effekt)", "Alert"))
print(plot_forest_multi(df_top10, "Top 10 Signifikante Landkreise (DiD Effekt)", "DiD"))
dev.off()

# Alle signifikanten Plots
pdf(
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/Check_Alle_SIGNIFIKANT_Alert_MultiModel.pdf",
  width = 8,
  height = 9
)
print(plot_forest_multi(df_plot, "Alle signifikanten AGS (Alert Effekt)", "Alert"))
dev.off()
pdf(
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/Check_Alle_SIGNIFIKANT_DiD_MultiModel.pdf",
  width = 8,
  height = 9
)
print(plot_forest_multi(df_plot, "Alle signifikanten AGS (DiD Effekt)", "DiD"))
dev.off()

# =====================================================================
# 4a. PLOT: TIME INVARIANT DLNM ZENTRIERT AUF 75. PERZENTIL (Model 1)
# =====================================================================

pdf(
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/DlnmFixCheck_SIGNIFIKANT_75Perc_def.pdf",
  width = 12,
  height = 14
)
layout(matrix(1:15, ncol=3, byrow=T))
par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.2, 0.8, 0), las=1)

for (current_ags in sig_ags_liste) {
  fdz_sub <- df_fdz %>% filter(AGS_Zahlen == current_ags)
  citydat <- df_temp %>% filter(AGS == current_ags)
  if(nrow(fdz_sub) == 0 || nrow(citydat) == 0) next
  
  coef1 <- get_coef(fdz_sub, "model1")
  vcov1 <- get_vcov(fdz_sub, "model1")
  if(is.null(coef1) || is.null(vcov1)) next
  
  argvar <- list(x = citydat$T_mean_popw, fun = varfun, degree = vardegree, 
                 knots = quantile(citydat$T_mean_popw, varper/100, na.rm=T))
  bvar <- do.call(onebasis, argvar)
  
  cen_temp <- quantile(citydat$T_mean_popw, 0.75, na.rm=T)
  
  cp <- crosspred(bvar, coef=coef1, vcov=vcov1, model.link="log", 
                  by=0.1, cen=cen_temp)
  
  plot(cp, ylim=c(0.9, max(cp$allRRhigh, na.rm=T)+0.1),
       xlim=c(min(citydat$T_mean_popw, na.rm=T)-1, max(citydat$T_mean_popw, na.rm=T)+1),
       xlab="Temperature", ylab="RR (Ref: 75th Perc)", lab=c(6,5,7), lwd=1, 
       main=paste("AGS:", current_ags), col="blue")
  
  abline(v=cen_temp, lty=3, lwd=0.9, col=grey(0.8))
}
dev.off()

# =====================================================================
# 4. DLNM PLOTS ZENTRIERT AUF 75. PERZENTIL
# =====================================================================

pdf(
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/DlnmVarCheck_SIGNIFIKANT_75Perc.pdf",
  width = 12,
  height = 14
)
layout(matrix(1:15, ncol=3, byrow=T))
par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.2, 0.8, 0), las=1)

for (current_ags in sig_ags_liste) {
  fdz_sub <- df_fdz %>% filter(AGS_Zahlen == current_ags)
  citydat <- df_temp %>% filter(AGS == current_ags)
  if(nrow(fdz_sub) == 0 || nrow(citydat) == 0) next
  
  coef2 <- get_coef(fdz_sub, "model2"); vcov2 <- get_vcov(fdz_sub, "model2")
  coef3 <- get_coef(fdz_sub, "model3"); vcov3 <- get_vcov(fdz_sub, "model3")
  if(is.null(coef2) || is.null(coef3)) next
  
  argvar <- list(x = citydat$T_mean_popw, fun = varfun, degree = vardegree, 
                 knots = quantile(citydat$T_mean_popw, varper/100, na.rm=T))
  bvar <- do.call(onebasis, argvar)
  
  cen_temp <- quantile(citydat$T_mean_popw, 0.75, na.rm=T)
  
  cp2 <- crosspred(bvar, coef=coef2, vcov=vcov2, model.link="log", by=0.1, cen=cen_temp)
  cp3 <- crosspred(bvar, coef=coef3, vcov=vcov3, model.link="log", by=0.1, cen=cen_temp)
  
  plot(cp2, ylim=c(0.9, max(c(cp2$allRRhigh, cp3$allRRhigh), na.rm=T)+0.1),
       xlim=c(min(citydat$T_mean_popw, na.rm=T)-1, max(citydat$T_mean_popw, na.rm=T)+1),
       xlab="Temperature", ylab="RR (Ref: 75th Perc)", lab=c(6,5,7), lwd=1.5, 
       main=paste("AGS:", current_ags), col="blue", ci.arg=list(density=20, col="blue"))
  lines(cp3, lwd=1.5, col="black", ci="area", ci.arg=list(density=20, angle=-45, col="black"))
  abline(v=cen_temp, lty=3, lwd=0.9, col=grey(0.8))
  legend("topleft", c("Periode 1 (Potentiell)", "Periode 2 (Echt)"), col=c("blue","black"), lwd=1.5, cex=0.8, inset=0.05, bty="n")
}
dev.off()

# =====================================================================
# 5. KARTENDARSTELLUNG (CHOROPLETH MAP) FÜR DEN ALERT-EFFEKT
# =====================================================================

print("Erstelle Deutschlandkarte...")

df_map_full <- df_fdz %>%
  filter(Model == "model1") %>%
  mutate(
    AGS_Zahlen = sprintf("%05d", as.numeric(AGS_Zahlen)), 
    ist_signifikant = ifelse(CI_Up_Alert < 1.0, TRUE, FALSE),
    Plot_RR = ifelse(ist_signifikant, RR_Alert, NA)
  )

deutschland_map <- st_read(shapefile_path)
map_merged <- left_join(deutschland_map, df_map_full, by = c("AGS" = "AGS_Zahlen")) 

pdf(
  "outputs/06_fdz_output_processing/06_02_heat_sensitive_districts/Deutschlandkarte_Alert_Model1_Neu.pdf",
  width = 8,
  height = 10
)
ggplot(data = map_merged) +
  geom_sf(aes(fill = Plot_RR), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "mako",       
    direction = -1,        
    name = "Relatives Risiko\n(Signifikant)",
    na.value = "grey85"   
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right", legend.title = element_text(face = "bold", size=10)) +
  labs(title = "Regionaler Effekt der Hitzewarnungen (Alert Effect)", subtitle = "Grau = Nicht signifikanter Effekt",
       caption = "Dargestellt ist Model 1 (kontrolliert für Temperatur)")
dev.off()
