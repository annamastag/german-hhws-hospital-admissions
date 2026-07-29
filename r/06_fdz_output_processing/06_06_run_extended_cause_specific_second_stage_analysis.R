# ==============================================================================
# Script: 06_06_run_extended_cause_specific_second_stage_analysis.R
#
# Purpose:
#   Conducts the extended cause-specific second-stage analysis including
#   cardiovascular admissions. District-specific first-stage estimates are
#   pooled across German administrative districts using random-effects
#   meta-analysis. The script reconstructs diagnosis-specific pointwise pooled
#   temperature-response curves and pools alert-day and Difference-in-
#   Differences estimates.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Released cause-specific first-stage FDZ output table including
#     cardiovascular admissions
#   - Environmental and heat-alert master dataset
#   - Population raster derived from the 2011 German Census
#   - Harmonized German district boundary shapefile
#
# Outputs:
#   - Diagnosis-specific pooled temperature-response curves
#   - District-level and pooled alert-day and Difference-in-Differences
#     forest plots
#   - CSV table of pooled cause-specific estimates
#   - Main cause-specific publication plot
#   - Supplementary model-comparison plot
#   - CSV table of admission counts and population-based admission rates
#   - Numerical results printed to the console
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
# Primary model:
#   Model 1 with a time-invariant temperature-response function and the RFM5
#   definition of reconstructed pre-implementation alert-qualifying days
#
# Required packages:
#   mixmeta, dplyr, ggplot2, readxl, dlnm, lubridate, terra, sf
#
# Dependencies:
#   Requires first-stage estimates generated within the secure FDZ environment
#   and the environmental exposure dataset prepared for the FDZ analysis.
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

rm(list=ls())

library(mixmeta)
library(dplyr)
library(ggplot2)
library(readxl) 
library(dlnm)
library(lubridate) 
library(terra)
library(sf)

# =====================================================================
# 1. PARAMETER & DATEN VORBEREITEN (FRANKENSTEIN-MERGE & TEMP-FIX)
# =====================================================================

varfun    <- "bs"              
vardegree <- 2                 
varper    <- c(50, 90)         
lag       <- 3                 
dfseas    <- 4                 
dftrend   <- 1                 
by_temperature <- 0.1

# A) Exporte lesen

df_all_cause <- read_excel("data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx") 
df_diagnosen <- read_excel("data/processed/fdz_outputs/02_Analysis_FirstStage_default_allCause_and_Strat(mit_Card)_und_sens_Res_gh.xlsx")

df_combined <- bind_rows(df_all_cause, df_diagnosen)

df <- df_combined %>%
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
    RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
    CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
    CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD))),
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = sprintf("%05d", as.numeric(gsub(".csv", "", AGS_Zahlen, fixed = TRUE)))
  )

# B) Masterfile lesen 

df_temp <- read.csv("data/external/fdz_input/Masterfile_Final_for_FDZ.csv", sep=",") %>%
  mutate(
    AGS = sprintf("%05d", as.numeric(AGS)), 
    Datum = as.Date(Datum, format="%Y-%m-%d"),
    Jahr = year(Datum),
    Monat = month(Datum) 
  ) %>%
  filter(Jahr >= 2000 & Jahr <= 2009, Monat >= 5 & Monat <= 9) 

ags_list <- unique(df$AGS_Zahlen)

df_plot <- df %>% 
  filter(Model %in% c("model0", "model1", "model2")) %>%
  mutate(
    Model_Name = factor(case_when(
      Model == "model0" ~ "Model0 (no temperature dlnm)",
      Model == "model1" ~ "Model1 (time-invariant temperature dlnm)",
      Model == "model2" ~ "Model2 (time-variant temperature dlnm)"
    ), levels = c("Model0 (no temperature dlnm)", "Model1 (time-invariant temperature dlnm)", "Model2 (time-variant temperature dlnm)")),
    Diagnosis = factor(Diagnosis, levels = c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital"))
  )

diagnosen <- levels(df_plot$Diagnosis)

# =====================================================================
# 2 & 3.  TEMPERATUR-POOLING (PUNKTWEISES POOLING PRO DIAGNOSE)
# =====================================================================

# Basis-Grid für alle Diagnosen berechnen (Temperaturverteilung ist ja gleich)
temp_all <- df_temp %>% filter(AGS %in% ags_list) %>% pull(T_mean_popw)
reference_temperature <- as.numeric(quantile(temp_all, probs = 0.75, na.rm = TRUE, names = FALSE))
national_p99 <- as.numeric(quantile(temp_all, probs = 0.99, na.rm = TRUE, names = FALSE))
plot_limits <- as.numeric(quantile(temp_all, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE))

temperature_grid <- sort(unique(c(
  seq(plot_limits[1], plot_limits[2], by = by_temperature),
  reference_temperature, national_p99
)))

# Extraktionsfunktion für Kontraste (angepasst auf spezifische Diagnose)
parse_pipe_numeric <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(as.character(x))) return(numeric(0))
  suppressWarnings(as.numeric(strsplit(gsub(",", ".", as.character(x), fixed = TRUE), "\\|")[[1]]))
}

get_diag_contrasts <- function(current_ags, cur_diag) {
  row_i <- df_plot %>% filter(AGS_Zahlen == current_ags, Diagnosis == cur_diag, Model == "model1")
  if(nrow(row_i) == 0) return(NULL)
  
  beta_i <- parse_pipe_numeric(row_i$Temp_Coef_Vector[[1]])
  vcov_vec_i <- parse_pipe_numeric(row_i$Temp_Vcov_Matrix[[1]])
  
  if (length(beta_i) == 0L || anyNA(beta_i)) return(NULL)
  
  V_i <- matrix(vcov_vec_i, nrow = length(beta_i), ncol = length(beta_i))
  temp_i <- df_temp %>% filter(AGS == current_ags) %>% arrange(Datum) %>% pull(T_mean_popw)
  
  knots_i <- as.numeric(quantile(temp_i, probs = varper / 100, na.rm = TRUE, names = FALSE))
  boundary_i <- range(temp_i, na.rm = TRUE)
  
  # Absolute Kurve
  if (reference_temperature < boundary_i[1] || reference_temperature > boundary_i[2]) {
    df_abs <- NULL
  } else {
    grid_supported_i <- temperature_grid[temperature_grid >= boundary_i[1] & temperature_grid <= boundary_i[2]]
    if (length(grid_supported_i) > 0L) {
      B_grid_i <- onebasis(x = grid_supported_i, fun = varfun, degree = vardegree, knots = knots_i, Boundary.knots = boundary_i, intercept = FALSE)
      B_ref_i <- onebasis(x = reference_temperature, fun = varfun, degree = vardegree, knots = knots_i, Boundary.knots = boundary_i, intercept = FALSE)
      contrast_matrix_abs <- sweep(B_grid_i, MARGIN = 2, STATS = as.numeric(B_ref_i[1, ]), FUN = "-")
      df_abs <- data.frame(
        AGS_Zahlen = current_ags, Temperature = grid_supported_i,
        Log_RR = as.numeric(contrast_matrix_abs %*% beta_i),
        Variance = rowSums((contrast_matrix_abs %*% V_i) * contrast_matrix_abs)
      )
    } else df_abs <- NULL
  }
  return(df_abs)
}

# Schleife über Diagnosen für Temperaturkurven
for (diag in diagnosen) {
  print(paste("Extrahiere und poole absolute Temperaturkurve für:", diag))
  safe_diag_name <- gsub(" & | |&", "_", diag)
  
  list_abs <- list()
  for (ags in ags_list) {
    res_abs <- tryCatch(get_diag_contrasts(ags, diag), error = function(e) NULL)
    if (!is.null(res_abs)) list_abs[[ags]] <- res_abs
  }
  
  df_abs_contrasts <- bind_rows(list_abs) %>% filter(is.finite(Log_RR), is.finite(Variance))
  
  # Punktweises Pooling
  pool_one_temp <- function(cur_temp) {
    dat_t <- df_abs_contrasts %>% filter(abs(Temperature - cur_temp) < 1e-8)
    if (abs(cur_temp - reference_temperature) < 1e-8) {
      return(data.frame(Temperature = cur_temp, Log_RR = 0, RR = 1, CI_Low = 1, CI_Up = 1))
    }
    dat_meta <- dat_t %>% filter(Variance > 0)
    if (nrow(dat_meta) < 2L) return(NULL)
    
    meta_t <- mixmeta(Log_RR ~ 1, S = Variance, data = dat_meta, method = "reml")
    pred_t <- predict(meta_t, ci = TRUE)
    data.frame(
      Temperature = cur_temp, Log_RR = as.numeric(pred_t[1, "fit"]), RR = exp(as.numeric(pred_t[1, "fit"])),
      CI_Low = exp(as.numeric(pred_t[1, "ci.lb"])), CI_Up = exp(as.numeric(pred_t[1, "ci.ub"]))
    )
  }
  
  pooled_curve <- bind_rows(lapply(temperature_grid, pool_one_temp)) %>% arrange(Temperature)
  
  # Absolute Kurve als PDF speichern
  pdf(
    paste0(
      "outputs/06_fdz_output_processing/06_06_extended_cause_specific_second_stage/figures/",
      "Temperature_Pooled_Curve_Model1_",
      safe_diag_name,
      ".pdf"
    ),
    width = 9,
    height = 5.5
  )
  p_curve <- ggplot(pooled_curve, aes(x = Temperature, y = RR)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_vline(xintercept = reference_temperature, linetype = "dotted", linewidth = 0.8, color = "grey40") +
    geom_ribbon(aes(ymin = CI_Low, ymax = CI_Up), alpha = 0.20, fill = "#d73027") +
    geom_line(linewidth = 1.1, color = "#a50026") +
    
    coord_cartesian(ylim = c(min(pooled_curve$CI_Low, na.rm = TRUE) * 0.99, 
                             max(pooled_curve$CI_Up, na.rm = TRUE) * 1.02)) + 
    
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Daily Mean Air Temperature (°C)",
      y = paste0("Relative Risk (Ref: ", round(reference_temperature, 1), " °C)"),
      title = paste("Pooled Temperature–Admission Association -", diag)
    )
  print(p_curve)
  dev.off()
  
  # Konsolen-Ausgabe für P99 (pro Diagnose)
  national_p99_result <- pooled_curve %>% filter(abs(Temperature - national_p99) < 1e-8)
  if(nrow(national_p99_result) > 0) {
    cat(sprintf("\nTemp-Effekt (P99 vs P75) für %s: RR %.3f [%.3f - %.3f]\n", 
                diag, national_p99_result$RR, national_p99_result$CI_Low, national_p99_result$CI_Up))
  }
}

# =====================================================================
# 4. FOREST PLOTS FÜR ALERT & DID (Getrennt nach Diagnose)
# =====================================================================

plot_pooled_effect <- function(data, effect_type, mod, diag) {
  dat_mod <- data %>% filter(Model == mod & Diagnosis == diag) %>% arrange(AGS_Zahlen)
  
  if(effect_type == "Alert") {
    rr_val <- dat_mod$RR_Alert; ci_l <- dat_mod$CI_Low_Alert; ci_u <- dat_mod$CI_Up_Alert
    xlab_text <- paste("RR (Alert Effect) -", mod)
  } else {
    rr_val <- dat_mod$RR_DiD; ci_l <- dat_mod$CI_Low_DiD; ci_u <- dat_mod$CI_Up_DiD
    xlab_text <- paste("RR (Alert*Implementation) -", mod)
  }
  
  safe_diag_name <- gsub(" & | |&", "_", diag)
  pdf_name <- paste0(
    "outputs/06_fdz_output_processing/06_06_extended_cause_specific_second_stage/figures/",
    "RR_",
    effect_type,
    "_blup_ALLE_AGS_",
    mod,
    "_",
    safe_diag_name,
    ".pdf"
  )
  
  yi <- log(rr_val)
  sei <- (log(ci_u) - log(ci_l)) / (2 * 1.96)
  vi <- sei^2
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  
  if(length(valid_idx) < 3) return(NULL) 
  
  yi <- yi[valid_idx]; vi <- vi[valid_idx]; cities_valid <- dat_mod$AGS_Zahlen[valid_idx]
  rr_val <- rr_val[valid_idx]; ci_l <- ci_l[valid_idx]; ci_u <- ci_u[valid_idx]
  
  meta_model <- mixmeta(yi ~ 1, S = vi, method = "reml")
  RRblup <- exp(blup(meta_model, pi=TRUE))
  RRpooled <- exp(predict(meta_model, ci=TRUE))
  
  nc <- length(cities_valid)
  yseq <- 1:(nc+2)
  ylab <- c("Pooled", "", rev(cities_valid))
  
  pdf(pdf_name, width=8, height=max(10, nc * 0.12))
  par(mar=c(5,6,3,2)) 
  
  plot(1, type="n", xlim=c(min(ci_l, na.rm=T)*0.95, max(ci_u, na.rm=T)*1.15), 
       ylim=c(0.5, nc+2.5), yaxs="i", yaxt="n", xlab=xlab_text, ylab="", mgp=c(2.5,1,0), main=paste("Diagnosis:", diag)) 
  axis(2, at=yseq, las=1, tick=TRUE, labels=ylab, cex.axis=0.5) 
  abline(v=1, lwd=1, col="black")
  legend("topright", c(paste0("First-stage ", mod), "BLUPs"), col=c("blue", "black"), pch=19, bty="o", bg="white", box.col="white", horiz=FALSE, cex=0.8)
  
  arrows(rev(ci_l),  yseq[c(-1,-2)]+0.1, rev(ci_u),  yseq[c(-1,-2)]+0.1, col="blue", code=3, angle=90, length=0.02, lwd=0.5)
  points(rev(rr_val), yseq[c(-1,-2)]+0.1, col="blue", pch=19, cex=0.6)
  arrows(rev(RRblup[,"pi.lb"]), yseq[c(-1,-2)]-0.1, rev(RRblup[,"pi.ub"]), yseq[c(-1,-2)]-0.1, col="black", code=3, angle=90, length=0.02, lwd=0.5)
  points(rev(RRblup[,1]), yseq[c(-1,-2)]-0.1, col="black", pch=19, cex=0.6) 
  arrows(RRpooled[1, "ci.lb"], 1, RRpooled[1, "ci.ub"], 1, col="orange", code=3, angle=90, length=0.05, lwd=2)
  points(RRpooled[1, "fit"], 1, col="orange", pch=19, cex=1.5) 
  dev.off()
}

# =====================================================================
# 5. POOLED ERGEBNISTABELLE 
# =====================================================================

get_pooled_numbers <- function(data, effect_type, mod, diag) {
  dat_mod <- data %>% filter(Model == mod & Diagnosis == diag) %>% arrange(AGS_Zahlen)
  
  if(effect_type == "Alert") {
    rr_val <- dat_mod$RR_Alert; ci_l <- dat_mod$CI_Low_Alert; ci_u <- dat_mod$CI_Up_Alert
  } else {
    rr_val <- dat_mod$RR_DiD; ci_l <- dat_mod$CI_Low_DiD; ci_u <- dat_mod$CI_Up_DiD
  }
  
  yi <- log(rr_val); sei <- (log(ci_u) - log(ci_l)) / (2 * 1.96); vi <- sei^2
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  
  if(length(valid_idx) < 3) {
    return(data.frame(Diagnosis = diag, Effect = effect_type, Model = mod, Pooled_RR = NA, Pooled_CI_Low = NA, Pooled_CI_Up = NA, I2_Heterogenity = NA))
  }
  
  meta_model <- mixmeta(yi[valid_idx] ~ 1, S = vi[valid_idx], method = "reml")
  RRpooled <- exp(predict(meta_model, ci=TRUE))
  i2_wert <- paste0(round(summary(meta_model)$i2stat[1], 1), " %")
  
  return(data.frame(Diagnosis = diag, Effect = effect_type, Model = mod, Pooled_RR = round(RRpooled[1, "fit"], 6), Pooled_CI_Low = round(RRpooled[1, "ci.lb"], 6), Pooled_CI_Up = round(RRpooled[1, "ci.ub"], 6), I2_Heterogenity = i2_wert))
}

print("Generiere gepoolte Tabelle...")
pooled_results_list <- list()
for(diag in diagnosen) {
  for(m in c("model0", "model1", "model2")) {
    pooled_results_list[[length(pooled_results_list) + 1]] <- get_pooled_numbers(df_plot, "Alert", m, diag)
    pooled_results_list[[length(pooled_results_list) + 1]] <- get_pooled_numbers(df_plot, "DiD", m, diag)
    plot_pooled_effect(df_plot, "Alert", m, diag)
    plot_pooled_effect(df_plot, "DiD", m, diag)
  }
}

Df_pooled_ALL <- bind_rows(pooled_results_list)
write.csv(
  Df_pooled_ALL,
  "outputs/06_fdz_output_processing/06_06_extended_cause_specific_second_stage/tables/Gepoolte_Ergebnisse_GanzDeutschland_Stratifiziert_inkl_Cardio.csv",
  row.names = FALSE
)

# =====================================================================
# 6. PUBLICATION PLOT (MAIN RESULTS - CAUSE-SPECIFIC, MODEL 1 ONLY)
# =====================================================================


# 1. Daten für den Plot vorbereiten

df_plot_paper <- Df_pooled_ALL %>%
  filter(Model == "model1") %>%
  filter(Diagnosis %in% c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital")) %>%
  mutate(
    
    # Signifikanz-Logik für die Rauten
    Precision_Status = case_when(
      is.na(Pooled_CI_Low) | is.na(Pooled_CI_Up) ~ "Unknown", 
      Pooled_CI_Low > 1 | Pooled_CI_Up < 1 ~ "Excludes 1.0",
      TRUE ~ "Includes 1.0"
    ),
    Precision_Status = factor(Precision_Status, levels = c("Excludes 1.0", "Includes 1.0", "Unknown")),
    
    # Einheitliche, schöne Y-Achsen Labels
    Diagnosis = case_when(
      Diagnosis == "Cardiovascular" ~ "Cardiovascular Diseases",
      Diagnosis == "Infectious & Parasitic" ~ "Infectious & Parasitic Diseases",
      Diagnosis == "Respiratory" ~ "Respiratory Diseases",
      Diagnosis == "Urogenital" ~ "Urogenital Diseases",
      TRUE ~ Diagnosis # "All-Cause" bleibt
    ),
    
    # Feste Sortierung (von oben nach unten im Plot: All-Cause ganz oben)
    Diagnosis = factor(Diagnosis, levels = rev(c(
      "All-Cause", 
      "Cardiovascular Diseases", 
      "Infectious & Parasitic Diseases", 
      "Respiratory Diseases", 
      "Urogenital Diseases"
    ))),

    Effect_Label = ifelse(Effect == "Alert", 
                          "Pre-Implementation Alert-Day Association", 
                          "Change in Alert-Day Association (DiD)"),
    Effect_Label = factor(Effect_Label, levels = c(
      "Pre-Implementation Alert-Day Association", 
      "Change in Alert-Day Association (DiD)"
    ))
  )

# 2. Den Plot im One-Panel Design zeichnen

pdf(
  "outputs/06_fdz_output_processing/06_06_extended_cause_specific_second_stage/figures/Plot_CauseSpecific_Model1_OnePanel.pdf",
  width = 11,
  height = 6
)

p_forest <- ggplot(df_plot_paper, aes(x = Pooled_RR, y = Diagnosis, color = Effect_Label, group = Effect_Label)) +
  
  # Null-Effekt Linie bei 1.0
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Fehlerbalken und Punkte zeichnen (position_dodge schiebt Blau und Rot leicht übereinander)
  geom_errorbarh(aes(xmin = Pooled_CI_Low, xmax = Pooled_CI_Up), height = 0.3, linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = Precision_Status), size = 4, stroke = 1.2, fill = "white", position = position_dodge(width = 0.5)) + 
  
  # Farben exakt matchen
  scale_color_manual(
    values = c("Pre-Implementation Alert-Day Association" = "#1f78b4", 
               "Change in Alert-Day Association (DiD)" = "#e31a1c")
  ) +
  
  # Rauten für Signifikanz (18 = gefüllt, 5 = hohl/durchsichtig)
  scale_shape_manual(
    values = c("Excludes 1.0" = 18, "Includes 1.0" = 5, "Unknown" = 4)
  ) +
  
  theme_bw() +
  labs(
    title = "Cause-Specific Alert-Day Associations and DiD Estimates",
    subtitle = "Pooled estimates from MODEL1 | National Average",
    x = "Relative Risk (RR) with 95% Confidence Interval",
    y = "",
    color = "Estimate",
    shape = "95% CI"
  ) +
  
  # Design-Anpassungen
  theme(
    legend.position = "top",
    legend.box = "vertical", 
    legend.margin = margin(t = 0, b = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    plot.title.position = "plot" 
  )

print(p_forest)
dev.off()

print("--> Publikations-Plot 'Plot_CauseSpecific_Model1_OnePanel.pdf' erfolgreich gespeichert!")

# =====================================================================
# 7. ZAHLEN-OUTPUT FÜR DAS MANUSKRIPT (CAUSE-SPECIFIC DiD - MODEL 1)
# =====================================================================

cat("\n\n=== ZAHLEN FÜR DAS MANUSKRIPT (CAUSE-SPECIFIC DiD - MODEL 1) ===\n")

for (diag in diagnosen) {
  dat_diag <- df_plot %>% filter(Model == "model1" & Diagnosis == diag & !is.na(RR_DiD))
  yi <- log(dat_diag$RR_DiD); sei <- (log(dat_diag$CI_Up_DiD) - log(dat_diag$CI_Low_DiD)) / (2 * 1.96); vi <- sei^2
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  
  if(length(valid_idx) < 3) {
    cat(sprintf("\n--- Diagnose: %s ---\nZu wenig Daten.\n", diag))
    next
  }
  
  meta_diag <- mixmeta(yi[valid_idx] ~ 1, S = vi[valid_idx], method = "reml")
  RRpooled <- exp(predict(meta_diag, ci=TRUE)); summ <- summary(meta_diag)
  
  cat(sprintf("\n--- Diagnose: %s ---\n", diag))
  cat(sprintf("Gepooltes RR (DiD): %.3f\n", RRpooled[1, "fit"]))
  cat(sprintf("95%% CI:             [%.3f - %.3f]\n", RRpooled[1, "ci.lb"], RRpooled[1, "ci.ub"]))
  cat(sprintf("I^2 (Heterog.):     %.1f %%\n", max(summ$i2stat)))
  cat(sprintf("Cochran's Q:        %.1f\n", summ$qstat$Q))
}
cat("\n=================================================================\n")

# =====================================================================
# 8. ABSOLUTE FALLZAHLEN & INZIDENZRATEN MIT ZENSUS-POPULATION (ALLE DIAGNOSEN)
# =====================================================================

cat("\n=== FALLZAHLEN & RATEN MIT ZENSUS 2011 (ALLE DIAGNOSEN) ===\n")

# 1. Zensus-TIF und AGS-Shapefile/Vektordaten laden

pop_raster <- rast(
  "data/raw/population/census_2011/population_1km.tif"
)
ags_polygons <- read_sf(
  "data/raw/spatial/vg2500/vg2500.shp"
)

# 2. Bevölkerung pro AGS extrahieren

ags_polygons <- ags_polygons %>% mutate(AGS_Zahlen = sprintf("%05d", as.numeric(AGS)))
pop_extracted <- terra::extract(pop_raster, ags_polygons, fun = sum, na.rm = TRUE)
ags_polygons$population <- pop_extracted[, 2] 

df_pop_lookup <- ags_polygons %>% st_drop_geometry() %>% select(AGS_Zahlen, population)

# 3. Datenbasis für Model 1 filtern und Population anfügen

df_all_analysis <- df_plot %>% 
  filter(Model == "model1") %>%
  mutate(AGS_Zahlen = sprintf("%05d", as.numeric(AGS_Zahlen))) %>%
  left_join(df_pop_lookup, by = "AGS_Zahlen")

# 4. Schleife über alle Diagnosen zur Berechnung der Raten

results_rates <- list()

for (diag in diagnosen) {
  
  df_diag <- df_all_analysis %>% filter(Diagnosis == diag)
  
  # Falls eine Diagnose komplett fehlt, überspringen
  if (nrow(df_diag) == 0) next
  
  # Fallzahlen aggregieren
  df_diag_cases <- df_diag %>%
    summarise(
      Sum_Total = sum(as.numeric(Total_Cases), na.rm = TRUE),
      Sum_Pre = sum(as.numeric(Total_Cases_Pre_2005), na.rm = TRUE),
      Sum_Post = sum(as.numeric(Total_Cases_Post_2005), na.rm = TRUE),
      Total_Study_Pop = sum(population, na.rm = TRUE) # Summe der Einwohner aller für DIESE Diagnose gültigen Kreise
    )
  
  pop_denominator <- df_diag_cases$Total_Study_Pop[1]
  
  # Jährliche Raten pro 100.000 Einwohner berechnen
  rate_total <- (df_diag_cases$Sum_Total[1] / 10) / (pop_denominator / 100000)
  rate_pre   <- (df_diag_cases$Sum_Pre[1] / 5) / (pop_denominator / 100000)
  rate_post  <- (df_diag_cases$Sum_Post[1] / 5) / (pop_denominator / 100000)
  
  # Konsolen-Output
  cat(sprintf("\n--- %s ---\n", toupper(diag)))
  cat(sprintf("Studienpopulation (Nenner): %s Einwohner\n", format(pop_denominator, big.mark=",", scientific=FALSE)))
  cat(sprintf("Fälle Total (2000-2009):    %s\n", format(df_diag_cases$Sum_Total[1], big.mark=",", scientific=FALSE)))
  cat(sprintf("Fälle Pre-2005:             %s\n", format(df_diag_cases$Sum_Pre[1], big.mark=",", scientific=FALSE)))
  cat(sprintf("Fälle Post-2005:            %s\n", format(df_diag_cases$Sum_Post[1], big.mark=",", scientific=FALSE)))
  cat(sprintf("Rate Total (p.a.):          %.1f / 100.000\n", rate_total))
  cat(sprintf("Rate Pre-2005 (p.a.):       %.1f / 100.000\n", rate_pre))
  cat(sprintf("Rate Post-2005 (p.a.):      %.1f / 100.000\n", rate_post))
  
  # Für die Export-Tabelle speichern
  results_rates[[diag]] <- data.frame(
    Diagnosis = diag,
    Study_Population = pop_denominator,
    Cases_Total = df_diag_cases$Sum_Total[1],
    Cases_Pre = df_diag_cases$Sum_Pre[1],
    Cases_Post = df_diag_cases$Sum_Post[1],
    Rate_Total = round(rate_total, 1),
    Rate_Pre = round(rate_pre, 1),
    Rate_Post = round(rate_post, 1)
  )
}

# 5. Export als CSV für das Paper
df_rates_export <- bind_rows(results_rates)
write.csv(
  df_rates_export,
  "outputs/06_fdz_output_processing/06_06_extended_cause_specific_second_stage/tables/Inzidenzraten_Alle_Diagnosen.csv",
  row.names = FALSE
)

cat("\n=====================================================\n")
print("--> Deskriptive Raten für alle Diagnosen berechnet und als CSV gespeichert!")

# =====================================================================
# SUPPLEMENTARY PLOT: MATRIX-LAYOUT (ALLE MODELLE X ALLE DIAGNOSEN)
# =====================================================================

df_plot_matrix <- Df_pooled_ALL %>%
  mutate(
    Model_Label = case_when(
      Model == "model0" ~ "Model 0\n(No temp. DLNM)",
      Model == "model1" ~ "Model 1\n(Time-invariant)",
      Model == "model2" ~ "Model 2\n(Time-varying)"
    ),
    Model_Label = factor(Model_Label, levels = rev(c(
      "Model 0\n(No temp. DLNM)", 
      "Model 1\n(Time-invariant)", 
      "Model 2\n(Time-varying)"
    ))),
    Effect_Label = ifelse(Effect == "Alert", "Baseline Risk (Alert-Qualifying Days)", "Added System Effect (DiD)"),
    Effect_Label = factor(Effect_Label, levels = c("Baseline Risk (Alert-Qualifying Days)", "Added System Effect (DiD)")),
    Diagnosis = factor(Diagnosis, levels = c("All-Cause", "Cardiovascular", "Infectious & Parasitic", "Respiratory", "Urogenital"))
  )

pdf(
  "outputs/06_fdz_output_processing/06_06_extended_cause_specific_second_stage/figures/Plot_Appendix_Matrix_AllModels_EN.pdf",
  width = 12,
  height = 9
)

ggplot(df_plot_matrix, aes(x = Pooled_RR, y = Model_Label, color = Effect_Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Pooled_CI_Low, xmax = Pooled_CI_Up), height = 0.2, linewidth = 1.2) +
  geom_point(size = 4, shape = 18) + 
  
  # Matrix-Design: Zeilen = Diagnose, Spalten = Effekt
  facet_grid(Diagnosis ~ Effect_Label, scales = "free_x") +
  
  scale_color_manual(values = c("Baseline Risk (Alert-Qualifying Days)" = "#2c7bb6", 
                                "Added System Effect (DiD)" = "#d7191c")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Sensitivity Analysis: Effect estimates across different temperature controls",
    x = "Relative Risk (RR) with 95% Confidence Interval", 
    y = ""
  ) +
  theme(
    legend.position = "none", 
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), 
    axis.text.y = element_text(face = "bold", color = "black")
  )

dev.off()
print("Matrix-Plot für das Supplement erfolgreich gespeichert!")

print("--> Alle stratifizierten Analysen, Tabellen und Plots (inklusive neuer Temperatur-Poolings) erfolgreich gespeichert!")
