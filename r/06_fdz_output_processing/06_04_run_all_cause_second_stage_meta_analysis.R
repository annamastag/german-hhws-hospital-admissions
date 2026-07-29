# ==============================================================================
# Script: 06_04_run_all_cause_second_stage_meta_analysis.R
#
# Purpose:
#   Conducts the national all-cause second-stage analysis by pooling
#   district-specific first-stage estimates using random-effects meta-analysis.
#   The script reconstructs a pointwise pooled temperature-response curve and
#   pools alert-day and Difference-in-Differences estimates across districts.
#
# Inputs:
#   - Released all-cause first-stage FDZ output table
#   - Environmental and heat-alert master dataset
#
# Outputs:
#   - CSV table of pooled alert-day and Difference-in-Differences estimates
#   - Pooled absolute temperature-response curve
#   - Temperature-percentile plot
#   - Pooled model-comparison plot
#   - Comparison plot of pre-implementation alert-day and DiD estimates
#   - Pooled numerical results printed to the console
#
# Primary analysis:
#   All-cause hospital admissions, RFM5, Model 1, May–September 2000–2009
#
# Required packages:
#   mixmeta, dplyr, ggplot2, readxl, dlnm, lubridate
#
# Dependencies:
#   Requires district-specific all-cause first-stage estimates and the
#   environmental master dataset.
#
# Repository note:
#   Only the file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

library(mixmeta)
library(dplyr)
library(ggplot2)
library(readxl)
library(dlnm)
library(lubridate)

rm(list = ls())

# =====================================================================
# 1. PARAMETER
# =====================================================================

varfun    <- "bs"
vardegree <- 2
varper    <- c(50, 90)

model_temperature <- "model1"
years_analysis    <- 2000:2009
months_analysis   <- 5:9
by_temperature    <- 0.1

# =====================================================================
# 2. DATEN EINLESEN UND PRÜFEN (Inkl. Jahres-Filter für Temp)
# =====================================================================

# A) First-Stage Exporte lesen & bereinigen
df <- read_excel(
  "data/processed/fdz_outputs/2026-05-13_Results_FirstStage_AllAdmissions_gh.xlsx"
) %>%
  filter(RFM == "heat_alerts_rfm5") %>%
  mutate(
    RR_Alert = as.numeric(gsub(",", ".", as.character(RR_Alert))),
    CI_Low_Alert = as.numeric(gsub(",", ".", as.character(CI_Low_Alert))),
    CI_Up_Alert = as.numeric(gsub(",", ".", as.character(CI_Up_Alert))),
    RR_DiD = as.numeric(gsub(",", ".", as.character(RR_DiD))),
    CI_Low_DiD = as.numeric(gsub(",", ".", as.character(CI_Low_DiD))),
    CI_Up_DiD = as.numeric(gsub(",", ".", as.character(CI_Up_DiD))),
    AGS_Zahlen = sprintf("%05d", as.numeric(gsub("AGS_|\\.csv", "", AGS_File)))
  )

if ("Diagnosis" %in% names(df)) {
  df <- df %>% filter(Diagnosis == "AllAdmissions")
}

firststage <- df %>% filter(Model == model_temperature)

# Für Plots und Sensitivitätsanalysen (alle Modelle)
df_plot <- df %>% filter(Model %in% c("model0", "model1", "model2")) %>%
  mutate(
    Model_Name = factor(case_when(
      Model == "model0" ~ "Model0 (no temperature dlnm)",
      Model == "model1" ~ "Model1 (time-invariant temperature dlnm)",
      Model == "model2" ~ "Model2 (time-variant temperature dlnm)"
    ), levels = c("Model0 (no temperature dlnm)", "Model1 (time-invariant temperature dlnm)", "Model2 (time-variant temperature dlnm)"))
  )

duplicate_check <- firststage %>% count(AGS_Zahlen, name = "N_Rows") %>% filter(N_Rows != 1L)
if (nrow(duplicate_check) > 0L) stop("Es liegt nicht genau eine Model-1-Zeile pro AGS vor.")

# B) Masterfile lesen 
temp_data <- read.csv(
  "data/external/fdz_input/Masterfile_Final_for_FDZ.csv",
  sep = ","
) %>%
  mutate(
    AGS   = sprintf("%05d", as.numeric(AGS)),
    Datum = as.Date(Datum, format = "%Y-%m-%d"),
    Jahr  = year(Datum),
    Monat = month(Datum)
  ) %>%
  filter(Jahr %in% years_analysis, Monat %in% months_analysis) %>%
  select(AGS, Datum, T_mean_popw) %>%
  distinct(AGS, Datum, .keep_all = TRUE) %>%
  filter(is.finite(T_mean_popw))

# =====================================================================
# 3. NATIONALE REFERENZEN & RASTER FÜR DIE ABSOLUTE KURVE
# =====================================================================

temp_all <- temp_data %>%
  filter(AGS %in% firststage$AGS_Zahlen) %>%
  pull(T_mean_popw)

reference_temperature <- as.numeric(quantile(temp_all, probs = 0.75, na.rm = TRUE, names = FALSE))
national_p99 <- as.numeric(quantile(temp_all, probs = 0.99, na.rm = TRUE, names = FALSE))
plot_limits <- as.numeric(quantile(temp_all, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE))

temperature_grid <- sort(unique(c(
  seq(plot_limits[1], plot_limits[2], by = by_temperature),
  reference_temperature, national_p99
)))

# =====================================================================
# 4. MATRIX-EXTRAKTION DER TEMPERATURKONTRASTE 
# =====================================================================

parse_pipe_numeric <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(as.character(x))) return(numeric(0))
  x <- gsub(",", ".", as.character(x), fixed = TRUE)
  suppressWarnings(as.numeric(strsplit(x, "\\|")[[1]]))
}

get_all_contrasts <- function(current_ags) {
  row_i <- firststage %>% filter(AGS_Zahlen == current_ags)
  
  beta_i <- parse_pipe_numeric(row_i$Temp_Coef_Vector[[1]])
  vcov_vec_i <- parse_pipe_numeric(row_i$Temp_Vcov_Matrix[[1]])
  
  if (length(beta_i) == 0L || anyNA(beta_i) || any(!is.finite(beta_i))) stop("Fehler Beta")
  if (length(vcov_vec_i) != length(beta_i)^2 || anyNA(vcov_vec_i) || any(!is.finite(vcov_vec_i))) stop("Fehler Vcov")
  
  V_i <- matrix(vcov_vec_i, nrow = length(beta_i), ncol = length(beta_i))
  temp_i <- temp_data %>% filter(AGS == current_ags) %>% arrange(Datum) %>% pull(T_mean_popw)
  
  knots_i <- as.numeric(quantile(temp_i, probs = varper / 100, na.rm = TRUE, names = FALSE))
  boundary_i <- range(temp_i, na.rm = TRUE)
  
  # A) LOKALER P99 vs P75 KONTRAST
  local_p75 <- as.numeric(quantile(temp_i, probs = 0.75, na.rm = TRUE, names = FALSE))
  local_p99 <- as.numeric(quantile(temp_i, probs = 0.99, na.rm = TRUE, names = FALSE))
  
  B_p99 <- onebasis(x = local_p99, fun = varfun, degree = vardegree, knots = knots_i, Boundary.knots = boundary_i, intercept = FALSE)
  B_p75 <- onebasis(x = local_p75, fun = varfun, degree = vardegree, knots = knots_i, Boundary.knots = boundary_i, intercept = FALSE)
  
  contrast_matrix_local <- B_p99 - B_p75
  df_local <- data.frame(
    AGS_Zahlen = current_ags, Local_P75 = local_p75, Local_P99 = local_p99,
    Log_RR = as.numeric(contrast_matrix_local %*% beta_i),
    Variance = as.numeric(contrast_matrix_local %*% V_i %*% t(contrast_matrix_local))
  )
  
  # B) ABSOLUTE KONTRASTE FÜR DIE KURVE
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
  return(list(local = df_local, absolute = df_abs))
}

print("Extrahiere lokale und absolute Kontraste für alle Landkreise...")
ags_list <- sort(unique(firststage$AGS_Zahlen))
list_local <- list(); list_abs <- list()

for (current_ags in ags_list) {
  res <- tryCatch(get_all_contrasts(current_ags), error = function(e) NULL)
  if (!is.null(res)) {
    list_local[[current_ags]] <- res$local
    if (!is.null(res$absolute)) list_abs[[current_ags]] <- res$absolute
  }
}

df_local_contrasts <- bind_rows(list_local) %>% filter(is.finite(Log_RR), is.finite(Variance), Variance > 0)
df_abs_contrasts <- bind_rows(list_abs) %>% filter(is.finite(Log_RR), is.finite(Variance))

# =====================================================================
# 5. POOLING DER TEMPERATUR (PUNKTWEISE) & TEXT-EXTRAKTION 
# =====================================================================

print("Poole absolute Kurve punktweise...")
pool_one_temp <- function(cur_temp) {
  dat_t <- df_abs_contrasts %>% filter(abs(Temperature - cur_temp) < 1e-8)
  if (abs(cur_temp - reference_temperature) < 1e-8) {
    return(data.frame(Temperature = cur_temp, Log_RR = 0, RR = 1, CI_Low = 1, CI_Up = 1, N_Districts = n_distinct(dat_t$AGS_Zahlen)))
  }
  dat_meta <- dat_t %>% filter(is.finite(Log_RR), is.finite(Variance), Variance > 0)
  if (nrow(dat_meta) < 2L) return(NULL)
  
  meta_t <- mixmeta(Log_RR ~ 1, S = Variance, data = dat_meta, method = "reml")
  pred_t <- predict(meta_t, ci = TRUE)
  data.frame(
    Temperature = cur_temp, Log_RR = as.numeric(pred_t[1, "fit"]), RR = exp(as.numeric(pred_t[1, "fit"])),
    CI_Low = exp(as.numeric(pred_t[1, "ci.lb"])), CI_Up = exp(as.numeric(pred_t[1, "ci.ub"])), N_Districts = n_distinct(dat_meta$AGS_Zahlen)
  )
}

pooled_curve <- bind_rows(lapply(temperature_grid, pool_one_temp)) %>% arrange(Temperature)
national_p99_result <- pooled_curve %>% filter(abs(Temperature - national_p99) < 1e-8)

meta_local <- mixmeta(Log_RR ~ 1, S = Variance, data = df_local_contrasts, method = "reml")
pred_local <- predict(meta_local, ci = TRUE)

cat("\n=== ZAHLEN FÜR DEN TEXT IN KAPITEL 3.2 (TEMPERATUR) ===\n")
cat("75. Perzentil (Baseline):", round(reference_temperature, 1), "°C\n")
cat("99. Perzentil (Extrem):  ", round(national_p99, 1), "°C\n")
cat("Gepooltes RR (Absolute Kurve P99 vs P75):", round(national_p99_result$RR, 3), "\n")
cat("95% CI Lower:            ", round(national_p99_result$CI_Low, 3), "\n")
cat("95% CI Upper:            ", round(national_p99_result$CI_Up, 3), "\n")
cat("--- ZUSATZINFO: LOKAL STANDARDISIERTER KONTRAST ---\n")
cat("Pooled RR:", round(exp(pred_local[1, "fit"]), 3), "[", round(exp(pred_local[1, "ci.lb"]), 3), "-", round(exp(pred_local[1, "ci.ub"]), 3), "]\n")
cat("=======================================================\n")

# =====================================================================
# 6. ALERT- UND DID-EFFEKTE (FOREST PLOTS & GEPOOLTE ZAHLEN)
# =====================================================================

print("Berechne DiD und Alert Effekte...")

get_pooled_numbers <- function(data, effect_type, mod) {
  dat_mod <- data %>% filter(Model == mod) %>% arrange(AGS_Zahlen)
  if(effect_type == "Alert") {
    rr_val <- dat_mod$RR_Alert; ci_l <- dat_mod$CI_Low_Alert; ci_u <- dat_mod$CI_Up_Alert
  } else {
    rr_val <- dat_mod$RR_DiD; ci_l <- dat_mod$CI_Low_DiD; ci_u <- dat_mod$CI_Up_DiD
  }
  
  yi <- log(rr_val); sei <- (log(ci_u) - log(ci_l)) / (2 * 1.96); vi <- sei^2
  valid_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0)
  
  meta_model <- mixmeta(yi[valid_idx] ~ 1, S = vi[valid_idx], method = "reml")
  RRpooled <- exp(predict(meta_model, ci=TRUE))
  i2_wert <- paste0(round(summary(meta_model)$i2stat[1], 1), " %")
  
  data.frame(Effect = effect_type, Model = mod, Pooled_RR = round(RRpooled[1, "fit"], 6), Pooled_CI_Low = round(RRpooled[1, "ci.lb"], 6), Pooled_CI_Up = round(RRpooled[1, "ci.ub"], 6), I2_Heterogenity = i2_wert)
}

Df_pooled_ALL <- bind_rows(lapply(c("model0", "model1", "model2"), function(m) {
  bind_rows(get_pooled_numbers(df_plot, "Alert", m), get_pooled_numbers(df_plot, "DiD", m))
}))
write.csv(
  Df_pooled_ALL,
  "outputs/06_fdz_output_processing/06_04_all_cause_second_stage/Gepoolte_Ergebnisse_GanzDeutschland.csv",
  row.names = FALSE
)

# --- ZAHLEN FÜR DEN FLIEßTEXT (DiD - HAUPTMODELL 1) ---

dat_did_m1 <- df %>% filter(Model == "model1") %>% arrange(AGS_Zahlen)
yi_did <- log(dat_did_m1$RR_DiD); sei_did <- (log(dat_did_m1$CI_Up_DiD) - log(dat_did_m1$CI_Low_DiD)) / (2 * 1.96); vi_did <- sei_did^2
valid_idx_did <- which(!is.na(yi_did) & !is.na(vi_did) & vi_did > 0)
meta_did_m1 <- mixmeta(yi_did[valid_idx_did] ~ 1, S = vi_did[valid_idx_did], method = "reml")
RRpooled_did <- exp(predict(meta_did_m1, ci=TRUE))
summ_did <- summary(meta_did_m1)

cat("\n=== ZAHLEN FÜR DEN TEXT IN KAPITEL 3.3 (HAUPTMODELL 1) ===\n")
cat("Gepooltes RR (DiD):", round(RRpooled_did[1, "fit"], 3), "\n")
cat("95% CI Lower:      ", round(RRpooled_did[1, "ci.lb"], 3), "\n")
cat("95% CI Upper:      ", round(RRpooled_did[1, "ci.ub"], 3), "\n")
cat("I^2 (Heterog.):    ", round(max(summ_did$i2stat), 1), "%\n")
cat("Cochran's Q:       ", round(summ_did$qstat$Q, 1), "\n")
cat("==========================================================\n")

# --- ZAHLEN FÜR SENSITIVITÄT (MODELS 0 & 2) ---

for (m in c("model0", "model2")) {
  dat_mod <- df %>% filter(Model == m) %>% arrange(AGS_Zahlen)
  yi <- log(dat_mod$RR_DiD); sei <- (log(dat_mod$CI_Up_DiD) - log(dat_mod$CI_Low_DiD)) / (2 * 1.96); vi <- sei^2
  v_idx <- which(!is.na(yi) & !is.na(vi) & vi > 0)
  meta_mod <- mixmeta(yi[v_idx] ~ 1, S = vi[v_idx], method = "reml")
  rr_p <- exp(predict(meta_mod, ci=TRUE)); summ <- summary(meta_mod)
  
  cat(sprintf("\n=== ZAHLEN FÜR %s (DiD EFFEKT) ===\n", toupper(m)))
  cat("Gepooltes RR:      ", round(rr_p[1, "fit"], 3), "\n")
  cat("95% CI Lower:      ", round(rr_p[1, "ci.lb"], 3), "\n")
  cat("95% CI Upper:      ", round(rr_p[1, "ci.ub"], 3), "\n")
  cat("I^2 (Heterog.):    ", round(max(summ$i2stat), 1), "%\n")
  cat("Cochran's Q:       ", round(summ$qstat$Q, 1), "\n")
  cat("=========================================\n")
}

# =====================================================================
# 7. ALLE PLOTS FÜR DAS PAPER
# =====================================================================

# A) TEMPERATUR-KURVE (Modell 1)

pdf(
  "outputs/06_fdz_output_processing/06_04_all_cause_second_stage/Temperature_Pooled_Absolute_Curve_Model1.pdf",
  width = 9,
  height = 5.5
)
print(ggplot(pooled_curve, aes(x = Temperature, y = RR)) +
        geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
        geom_vline(xintercept = reference_temperature, linetype = "dotted", linewidth = 0.8, color = "grey40") +
        geom_ribbon(aes(ymin = CI_Low, ymax = CI_Up), alpha = 0.20, fill = "#d73027") +
        geom_line(linewidth = 1.1, color = "#a50026") +
        coord_cartesian(ylim = c(NA, 1.04)) +
        theme_bw(base_size = 14) +
        theme(panel.grid.minor = element_blank()) +
        labs(x = "Daily Mean Air Temperature (°C)", y = paste0("Relative Risk (Ref: ", round(reference_temperature, 1), " °C)"), title = "Pooled Temperature–Admission Association"))
dev.off()

# B) PERZENTIL-PLOT FÜR TEMPERATUR (NEUE LOGIK, BASIEREND AUF POOLED_CURVE)
# Perzentil-Funktion basierend auf ALLEN nationalen Temperaturen (2000-2009)

percentile_funktion <- ecdf(temp_all)
pooled_curve$temp_percentile <- percentile_funktion(pooled_curve$Temperature) * 100

temp_breaks <- c(10, 20, 30)
perc_breaks <- percentile_funktion(temp_breaks) * 100

pdf(
  "outputs/06_fdz_output_processing/06_04_all_cause_second_stage/Gesamteffekt_ALLE_AGS_Percentiles.pdf",
  width = 9,
  height = 5.5
)
print(ggplot(pooled_curve, aes(x = temp_percentile, y = RR)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", linewidth=0.8) +
        geom_vline(xintercept = 75, linetype = "dotted", color = "grey50", linewidth=1) + 
        geom_ribbon(aes(ymin = CI_Low, ymax = CI_Up), fill = "#d7191c", alpha = 0.2) +
        geom_line(color = "#d7191c", linewidth = 1.2) +
        theme_minimal(base_size = 14) +
        labs(x = "Temperature Percentiles", y = "Relative Risk (RR)") +
        coord_cartesian(xlim = c(0, 100), ylim = c(0.96, max(pooled_curve$CI_Up, na.rm=TRUE) * 1.01)) +
        scale_x_continuous(breaks = seq(0, 100, by = 10), sec.axis = sec_axis(~ ., breaks = perc_breaks, labels = paste0(temp_breaks, "°C"), name = "Absolute Temperature")) +
        theme(plot.margin = margin(10, 25, 10, 10), panel.grid.minor = element_blank(), axis.title.x.top = element_text(margin = margin(b = 10), color = "grey40", face = "bold"), axis.text.x.top = element_text(color = "grey40", face = "bold"), axis.text.x.bottom = element_text(color = "black"), axis.title.x.bottom = element_text(margin = margin(t = 10))))
dev.off()

# C) ÜBERSICHTS-PLOT: POOLED EFFECTS (DiD vs. ALERT) FÜR ALLE MODELLE

df_plot_paper <- Df_pooled_ALL %>%
  mutate(
    Model_Label = factor(case_when(
      Model == "model0" ~ "Model 0\n(No temp. DLNM)", Model == "model1" ~ "Model 1\n(Time-invariant temp. DLNM)", Model == "model2" ~ "Model 2\n(Time-varying temp. DLNM)"
    ), levels = rev(c("Model 0\n(No temp. DLNM)", "Model 1\n(Time-invariant temp. DLNM)", "Model 2\n(Time-varying temp. DLNM)"))),
    Effect_Label = ifelse(Effect == "Alert", "Main Warning Effect (Alert)", "Added System Effect (DiD)")
  )

pdf(
  "outputs/06_fdz_output_processing/06_04_all_cause_second_stage/Plot_Pooled_Effects_Germany_EN.pdf",
  width = 9,
  height = 5
)
print(ggplot(df_plot_paper, aes(x = Pooled_RR, y = Model_Label, color = Effect_Label)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.8) +
        geom_errorbarh(aes(xmin = Pooled_CI_Low, xmax = Pooled_CI_Up), height = 0.2, linewidth = 1.2) +
        geom_point(size = 6, shape = 18) + 
        facet_wrap(~ Effect_Label, scales = "free_x") +
        scale_color_manual(values = c("Main Warning Effect (Alert)" = "#2c7bb6", "Added System Effect (DiD)" = "#d7191c")) +
        theme_bw(base_size = 14) +
        labs(x = "Relative Risk (RR) with 95% Confidence Interval", y = "") +
        theme(legend.position = "none", strip.background = element_rect(fill = "grey90", color = "white"), strip.text = element_text(face = "bold", size = 12), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold", color = "black")))
dev.off()

# D) PLACEBO VS. DiD VERGLEICHS-PLOT
# ALERT (Pre-Intervention / Placebo) extrahieren

dat_pre <- df %>% filter(Model == "model1") %>% arrange(AGS_Zahlen)
yi_pre <- log(dat_pre$RR_Alert); sei_pre <- (log(dat_pre$CI_Up_Alert) - log(dat_pre$CI_Low_Alert)) / (2 * 1.96); vi_pre <- sei_pre^2
v_pre <- which(!is.na(yi_pre) & !is.na(vi_pre) & vi_pre > 0)
m_pre <- mixmeta(yi_pre[v_pre] ~ 1, S = vi_pre[v_pre], method = "reml")
pred_pre <- predict(m_pre, ci=TRUE)

df_vergleich <- data.frame(
  Modell = factor(c("Simulated Alert (Placebo-Check, 2000-2004)", "Official Alert (DiD Estimator, 2005-2009)"), levels = c("Official Alert (DiD Estimator, 2005-2009)", "Simulated Alert (Placebo-Check, 2000-2004)")),
  RR = c(exp(pred_pre[1, "fit"]), RRpooled_did[1, "fit"]),
  CI_low = c(exp(pred_pre[1, "ci.lb"]), RRpooled_did[1, "ci.lb"]),
  CI_high = c(exp(pred_pre[1, "ci.ub"]), RRpooled_did[1, "ci.ub"])
)

pdf(
  "outputs/06_fdz_output_processing/06_04_all_cause_second_stage/Results_ForestPlot_Placebo_vs_DiD.pdf",
  width = 8,
  height = 3.5
)
print(ggplot(df_vergleich, aes(x = RR, y = Modell, color = Modell)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 1) +
        geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 1.2) +
        geom_point(size = 6, shape = 18) +
        scale_color_manual(values = c("Simulated Alert (Placebo-Check, 2000-2004)" = "#2c7bb6", "Official Alert (DiD Estimator, 2005-2009)" = "#d7191c")) +
        theme_minimal(base_size = 14) +
        labs(x = "Relative Risk (RR) with 95% Confidence Interval", y = "") +
        theme(legend.position = "none", axis.text.y = element_text(face = "bold", color = "black", size = 12), axis.text.x = element_text(size = 12), axis.title.x = element_text(margin = margin(t = 12)), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()))
dev.off()

print("--> All analyses and plots were saved to the configured output directory.")
