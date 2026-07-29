# ==============================================================================
# Script: 06_05_run_cause_specific_second_stage_meta_analysis.R
#
# Purpose:
#   Conducts the earlier cause-specific second-stage analysis by pooling
#   district-specific temperature-response coefficients and alert-related
#   estimates across districts using random-effects meta-analysis.
#
# Inputs:
#   - Released cause-specific first-stage FDZ output table
#   - Environmental and heat-alert master dataset
#
# Outputs:
#   - Diagnosis-specific pooled temperature-response curves
#   - Diagnosis-specific BLUP spaghetti plots
#   - District-level and pooled alert-day and DiD forest plots
#   - CSV table of pooled cause-specific estimates
#   - Cause-specific publication plots
#   - Numerical meta-analysis results printed to the console
#
# Outcomes:
#   - Respiratory admissions
#   - Urogenital admissions
#   - Infectious and parasitic disease admissions
#
# Required packages:
#   mixmeta, dplyr, ggplot2, readxl, dlnm, lubridate
#
# Dependencies:
#   Requires the earlier cause-specific first-stage FDZ output and the
#   environmental master dataset.
#
# Repository note:
#   This script represents an earlier cause-specific analysis version without
#   the later cardiovascular extension. The file name, path strings, and
#   explanatory comments were standardized.
# ==============================================================================

library(mixmeta)
library(dplyr)
library(ggplot2)
library(readxl) 
library(dlnm)
library(lubridate) 

# =====================================================================
# 1. PARAMETER & DATEN VORBEREITEN
# =====================================================================

varfun    <- "bs"              
vardegree <- 2                 
varper    <- c(50, 90)         
lag       <- 3                 
dfseas    <- 4                 
dftrend   <- 1                 

# Daten einlesen (mit Diagnose-Spalte)
df <- read_excel(
  "data/processed/fdz_outputs/2026-05-19_Results_FirstStage_Diagnosen_gh.xlsx"
) %>%
  filter(RFM == "heat_alerts_rfm5") %>%
  rename(Diagnosis = any_of(c("Diagnosis", "Diagnose", "Diagnosegruppe", "diagnosis"))) %>%
  mutate(
    RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
    CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
    CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD)))
  ) %>%
  mutate(
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = gsub(".csv", "", AGS_Zahlen, fixed = TRUE)
  )

# Wetterdaten laden 
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

ags_list <- unique(df$AGS_Zahlen)

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
    Diagnosis = factor(Diagnosis)
  )

diagnosen <- levels(df_plot$Diagnosis)

# Hilfsfunktionen

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
# 2 & 3. META-ANALYSE, GEPOLTER EFFEKT & SPAGHETTI PLOT (PRO DIAGNOSE)
# =====================================================================

for(diag in diagnosen) {
  
  coef_list <- list()
  vcov_list <- list()
  gueltige_ags <- c()
  
  for (ags in ags_list) {
    # WICHTIG: Filter auf AGS *und* Diagnose!
    fdz_sub <- df %>% filter(AGS_Zahlen == ags & Diagnosis == diag)
    c1 <- get_coef(fdz_sub, "model1")
    v1 <- get_vcov(fdz_sub, "model1")
    
    if(!is.null(c1) && !is.null(v1)) {
      coef_list[[ags]] <- c1
      vcov_list[[ags]] <- v1
      gueltige_ags <- c(gueltige_ags, ags)
    }
  }
  
  coef_mat <- do.call(rbind, coef_list)
  meta_model <- mixmeta(coef_mat ~ 1, S = vcov_list, method = "reml")
  pooled_coef <- coef(meta_model)
  pooled_vcov <- vcov(meta_model)
  blup_mat <- blup(meta_model)
  
  # Temperatur-Setup
  temp_all <- df_temp %>% filter(AGS %in% gueltige_ags) %>% pull(T_mean_popw)
  perc75_temp_all <- quantile(temp_all, 0.75, na.rm=TRUE)
  perc99_temp_all <- quantile(temp_all, 0.99, na.rm=TRUE)
  
  argvar_pool <- list(x = temp_all, fun = varfun, degree = vardegree, knots = quantile(temp_all, varper/100, na.rm=TRUE))
  bvar_pool <- do.call(onebasis, argvar_pool)
  temp_grid <- seq(min(temp_all, na.rm=TRUE), max(temp_all, na.rm=TRUE), by=0.1)
  
  # Vorhersage
  cp_pool_75 <- crosspred(bvar_pool, coef=pooled_coef, vcov=pooled_vcov, 
                          model.link="log", at=temp_grid, cen=perc75_temp_all)
  
  # --- PLOT A: Einzelne gepoolte Kurve ---
  
  pdf(
    paste0(
      "outputs/06_fdz_output_processing/06_05_cause_specific_second_stage/",
      "Gesamteffekt_ALLE_AGS_75Perc_",
      diag,
      ".pdf"
    ),
    width = 8,
    height = 6
  )
  plot(cp_pool_75, ylim=c(0.85, max(cp_pool_75$allRRhigh, na.rm=T)+0.1),
       xlim=c(min(temp_all, na.rm=T), max(temp_all, na.rm=T)),
       xlab="Temperatur (°C)", ylab="Relative Risk (RR) - Ref: 75th Percentile", 
       main=paste("Pooled national effect -", diag), 
       col="darkblue", lwd=2, ci.arg=list(density=30, col="darkblue"))
  abline(h=1, lty=1, col="black"); abline(v=perc75_temp_all, lty=3, lwd=1.5, col="grey")
  legend("topleft", legend=paste("Reference:", round(perc75_temp_all, 1), "°C"), bty="n")
  dev.off()
  
  # --- PLOT B: Spaghetti Plot ---
  
  dummy_vcov <- matrix(0, nrow = ncol(blup_mat), ncol = ncol(blup_mat))
  pdf(
    paste0(
      "outputs/06_fdz_output_processing/06_05_cause_specific_second_stage/",
      "Appendix_SpaghettiPlot_Germany_75Perc_",
      diag,
      ".pdf"
    ),
    width = 10,
    height = 7
  )
  par(mar=c(4.5, 4.5, 3, 2), mgp=c(2.5, 0.8, 0), las=1)
  plot(1, type="n", xlim=c(min(temp_all, na.rm=T), max(temp_all, na.rm=T)), ylim=c(0.85, 1.30), yaxs="i",  
       xlab="Temperature (°C)", ylab="Relative Risk (RR) - Ref: 75th Percentile", 
       main=paste("Heterogeneity of Temperature-Response Curves -", diag), cex.main=1.1)
  
  for(i in 1:nrow(blup_mat)) {
    cp_local <- crosspred(bvar_pool, coef=as.numeric(blup_mat[i, ]), vcov=dummy_vcov, model.link="log", at=temp_grid, cen=perc75_temp_all)
    lines(cp_local$predvar, cp_local$allRRfit, col=adjustcolor("grey40", alpha.f=0.2), lwd=1)
  }
  abline(h=1, col="black"); abline(v=perc75_temp_all, lty=3, col="grey40") 
  lines(cp_pool_75$predvar, cp_pool_75$allRRfit, col="#d7191c", lwd=3.5)
  lines(cp_pool_75$predvar, cp_pool_75$allRRlow, col="#d7191c", lty=2, lwd=1.5)
  lines(cp_pool_75$predvar, cp_pool_75$allRRhigh, col="#d7191c", lty=2, lwd=1.5)
  legend("topleft", legend=c("Pooled Estimate", "95% CI", "Regional Estimates (400 AGS)"), col=c("#d7191c", "#d7191c", "grey50"), lty=c(1, 2, 1), lwd=c(3.5, 1.5, 1.2), bty="n", inset=0.02)
  dev.off()
  
  print(paste("Kurven & Spaghetti-Plots für", diag, "gespeichert."))
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
  
  pdf_name <- paste0(
    "outputs/06_fdz_output_processing/06_05_cause_specific_second_stage/",
    "RR_",
    effect_type,
    "_blup_ALLE_AGS_",
    mod,
    "_",
    diag,
    ".pdf"
  )
  
  yi <- log(rr_val)
  sei <- (log(ci_u) - log(ci_l)) / (2 * 1.96)
  vi <- sei^2
  
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  
  # === SICHERHEITSCHECK FÜR ZU WENIG DATEN ===
  
  if(length(valid_idx) < 3) {
    warning(paste("Überspringe", effect_type, "für", mod, "[", diag, "]: Nur", length(valid_idx), "gültige Landkreise vorhanden! Data gaps too large."))
    return(NULL) # Beendet die Funktion vorzeitig für dieses Modell, ohne abzustürzen
  }
  
  # ==========================================================
  
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
       ylim=c(0.5, nc+2.5), yaxs="i", yaxt="n", 
       xlab=xlab_text, ylab="", mgp=c(2.5,1,0), 
       main=paste("Diagnosis:", diag)) 
  
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
# 5. POOLED ERGEBNISTABELLE (mit Diagnose-Spalte & Sicherheits-Check)
# =====================================================================

get_pooled_numbers <- function(data, effect_type, mod, diag) {
  dat_mod <- data %>% filter(Model == mod & Diagnosis == diag) %>% arrange(AGS_Zahlen)
  
  if(effect_type == "Alert") {
    rr_val <- dat_mod$RR_Alert; ci_l <- dat_mod$CI_Low_Alert; ci_u <- dat_mod$CI_Up_Alert
  } else {
    rr_val <- dat_mod$RR_DiD; ci_l <- dat_mod$CI_Low_DiD; ci_u <- dat_mod$CI_Up_DiD
  }
  
  # Logarithmieren und auf unendliche/fehlende Werte prüfen
  
  yi <- log(rr_val)
  sei <- (log(ci_u) - log(ci_l)) / (2 * 1.96)
  vi <- sei^2
  
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  
  if(length(valid_idx) < 3) {
    warning(paste("Tabelle: Zu wenig Daten für", effect_type, "-", mod, "[", diag, "]. Fülle mit NA auf."))
    return(data.frame(
      Diagnosis = diag,
      Effect = effect_type, 
      Model = mod, 
      Pooled_RR = NA, 
      Pooled_CI_Low = NA, 
      Pooled_CI_Up = NA, 
      I2_Heterogenity = NA
    ))
  }
  
  meta_model <- mixmeta(yi[valid_idx] ~ 1, S = vi[valid_idx], method = "reml")
  RRpooled <- exp(predict(meta_model, ci=TRUE))
  i2_wert <- paste0(round(summary(meta_model)$i2stat[1], 1), " %")
  
  return(data.frame(
    Diagnosis = diag, 
    Effect = effect_type, 
    Model = mod, 
    Pooled_RR = round(RRpooled[1, "fit"], 6), 
    Pooled_CI_Low = round(RRpooled[1, "ci.lb"], 6), 
    Pooled_CI_Up = round(RRpooled[1, "ci.ub"], 6), 
    I2_Heterogenity = i2_wert
  ))
}

print("Generiere gepoolte Tabelle...")
pooled_results_list <- list()
for(diag in diagnosen) {
  for(m in c("model0", "model1", "model2")) {
    pooled_results_list[[length(pooled_results_list) + 1]] <- get_pooled_numbers(df_plot, "Alert", m, diag)
    pooled_results_list[[length(pooled_results_list) + 1]] <- get_pooled_numbers(df_plot, "DiD", m, diag)
  }
}

Df_pooled_ALL <- bind_rows(pooled_results_list)
print("+++ GEPOLTE ERGEBNISSE FÜR GANZ DEUTSCHLAND (STRATIFIZIERT) +++")
print(Df_pooled_ALL)
write.csv(
  Df_pooled_ALL,
  "outputs/06_fdz_output_processing/06_05_cause_specific_second_stage/Gepoolte_Ergebnisse_GanzDeutschland_Stratifiziert.csv",
  row.names = FALSE
)

# =====================================================================
# 6. PUBLICATION PLOT (MATRIX-LAYOUT FÜR PAPER)
# =====================================================================

df_plot_paper <- Df_pooled_ALL %>%
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
    Effect_Label = ifelse(Effect == "Alert", "Main Warning Effect (Alert)", "Added System Effect (DiD)")
  )

pdf(
  "outputs/06_fdz_output_processing/06_05_cause_specific_second_stage/Plot_Pooled_Effects_Germany_Stratified_EN.pdf",
  width = 10,
  height = 7
)

ggplot(df_plot_paper, aes(x = Pooled_RR, y = Model_Label, color = Effect_Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Pooled_CI_Low, xmax = Pooled_CI_Up), height = 0.2, linewidth = 1.2) +
  geom_point(size = 6, shape = 18) + 

  # Zeilen = Diagnosen, Spalten = Effekte (Alert / DiD)
  facet_grid(Diagnosis ~ Effect_Label, scales = "free_x") +
  
  scale_color_manual(values = c("Main Warning Effect (Alert)" = "#2c7bb6", 
                                "Added System Effect (DiD)" = "#d7191c")) +
  theme_bw(base_size = 14) +
  labs(x = "Relative Risk (RR) with 95% Confidence Interval", y = "") +
  theme(
    legend.position = "none", 
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), 
    axis.text.y = element_text(face = "bold", color = "black")
  )
dev.off()

#######################################################################

# 1. Diagnosenamen direkt aus dem rohen Datensatz ziehen

diagnosen <- unique(df$Diagnosis)
diagnosen <- diagnosen[!is.na(diagnosen)] # Leere Einträge rauswerfen

cat("\n=== ZAHLEN FÜR DEN TEXT IN KAPITEL 4.4.2 (CAUSE-SPECIFIC DiD) ===\n")

# 2. Schleife nur für die reinen Text-Zahlen

for (diag in diagnosen) {
  
  # Daten filtern auf Model 1
  dat_diag <- df %>% filter(Model == "model1" & Diagnosis == diag) %>% arrange(AGS_Zahlen)
  
  yi <- log(dat_diag$RR_DiD)
  sei <- (log(dat_diag$CI_Up_DiD) - log(dat_diag$CI_Low_DiD)) / (2 * 1.96)
  vi <- sei^2
  
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0 & !is.infinite(yi) & !is.infinite(vi))
  
  if(length(valid_idx) < 3) {
    cat(sprintf("\n--- Diagnose: %s ---\n", diag))
    cat("Zu wenig gültige Daten für Meta-Analyse.\n")
    next
  }
  
  # Meta-Modell rechnen
  meta_diag <- mixmeta(yi[valid_idx] ~ 1, S = vi[valid_idx], method = "reml")
  
  # Werte ziehen
  RRpooled <- exp(predict(meta_diag, ci=TRUE))
  summ <- summary(meta_diag)
  
  # Output für die Konsole
  cat(sprintf("\n--- Diagnose: %s ---\n", diag))
  cat("Gepooltes RR:      ", round(RRpooled[1, "fit"], 3), "\n")
  cat("95% CI Lower:      ", round(RRpooled[1, "ci.lb"], 3), "\n")
  cat("95% CI Upper:      ", round(RRpooled[1, "ci.ub"], 3), "\n")
  cat("I^2 (Heterog.):    ", round(max(summ$i2stat), 1), "%\n")
  cat("Cochran's Q:       ", round(summ$qstat$Q, 1), "\n")
}
cat("=================================================================\n")

#######################################################################

# =====================================================================
# 6. PUBLICATION PLOT (MAIN RESULTS - CAUSE-SPECIFIC, MODEL 1 ONLY)
# =====================================================================

df_plot_paper <- Df_pooled_ALL %>%
  
  # Nur das Hauptmodell für den Text filtern
  filter(Model == "model1") %>%
  
  # (Optional: Falls "AllAdmissions" im df ist, hier nur die 3 Diagnosen filtern)
  filter(Diagnosis %in% c("Resp_0", "Uri_0", "Infect_par")) %>%
  mutate(
    Diagnosis = case_when(
      Diagnosis == "Resp_0" ~ "Respiratory",
      Diagnosis == "Uri_0" ~ "Urogenital",
      Diagnosis == "Infect_par" ~ "Infectious and Parasitic"
    ),
    # Y-Achse sortieren (von oben nach unten im Plot)
    Diagnosis = factor(Diagnosis, levels = rev(c("Respiratory", "Urogenital", "Infectious and Parasitic"))),
    
    # Spalten tauschen
    Effect_Label = ifelse(Effect == "Alert", "Main Warning Effect (Alert)", "Added System Effect (DiD)"),
    Effect_Label = factor(Effect_Label, levels = c("Main Warning Effect (Alert)", "Added System Effect (DiD)"))
  )

pdf(
  "outputs/06_fdz_output_processing/06_05_cause_specific_second_stage/Plot_CauseSpecific_Model1_EN.pdf",
  width = 10,
  height = 4
)

ggplot(df_plot_paper, aes(x = Pooled_RR, y = Diagnosis, color = Effect_Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Pooled_CI_Low, xmax = Pooled_CI_Up), height = 0.2, linewidth = 1.2) +
  geom_point(size = 6, shape = 18) + 
  
  # Zwei saubere Spalten für die beiden Effekte
  facet_wrap(~ Effect_Label, scales = "free_x") +
  
  scale_color_manual(values = c("Main Warning Effect (Alert)" = "#2c7bb6", 
                                "Added System Effect (DiD)" = "#d7191c")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Cause-Specific Effects of Heat Alerts",
    subtitle = "Pooled Relative Risks for highly heat-sensitive diagnosis groups (Model 1)",
    x = "Relative Risk (RR) with 95% Confidence Interval",
    y = ""
  ) +
  theme(
    legend.position = "none", 
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
    axis.text.y = element_text(face = "bold", color = "black", size = 12)
  )

dev.off()
