# ==============================================================================
# Script: 06_07_05_compare_period_specific_temperature_associations.R
#
# Purpose:
#   Compares pooled temperature-admission associations between the
#   pre-implementation and post-implementation periods. District-specific
#   Model 1 temperature coefficients are pooled separately for each period,
#   and the national relative risk at the 99th temperature percentile is
#   calculated relative to the 75th percentile.
#
# Inputs:
#   - Released period-specific first-stage FDZ output table
#   - Environmental and heat-alert master dataset
#
# Outputs:
#   - Pooled pre-implementation and post-implementation temperature estimates
#     printed to the console
#
# Outcome:
#   All-cause hospital admissions
#
# Analysis periods:
#   Pre-implementation:  2000–2004
#   Post-implementation: 2005–2009
#
# Model specification:
#   RFM5 and Model 1
#
# Required packages:
#   mixmeta, dplyr, readxl, dlnm, lubridate
#
# Repository note:
#   Personal absolute paths were replaced with neutral project-relative paths.
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

library(mixmeta)
library(dplyr)
library(readxl) 
library(dlnm)
library(lubridate) 

# 1. PARAMETER

varfun    <- "bs"              
vardegree <- 2                 
varper    <- c(50, 90)         

# 2. DATEN EINLESEN (Perioden-Split)

df_period <- read_excel(
  "data/processed/fdz_outputs/Results_FirstStage_GetrenntePerioden_gh.xlsx"
) %>%
  filter(RFM == "heat_alerts_rfm5" & Diagnosis == "AllAdmissions" & Model == "model1") %>%
  mutate(
    AGS_Zahlen = gsub("AGS_", "", AGS_File),
    AGS_Zahlen = gsub(".csv", "", AGS_Zahlen, fixed = TRUE)
  )

# 3. WETTERDATEN FÜR DIE REFERENZ-TEMPERATUREN 

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

temp_all <- df_temp$T_mean_popw
perc75_temp_all <- quantile(temp_all, 0.75, na.rm=TRUE)
perc99_temp_all <- quantile(temp_all, 0.99, na.rm=TRUE)

# Basis für die Temperaturkurve bauen
argvar_pool <- list(x = temp_all, fun = varfun, degree = vardegree, knots = quantile(temp_all, varper/100, na.rm=TRUE))
bvar_pool <- do.call(onebasis, argvar_pool)

# 4. HILFSFUNKTIONEN ZUM ENTPACKEN DER MATRIZEN

get_coef <- function(data) {
  val <- data$Temp_Coef_Vector
  if(length(val) == 0 || is.na(val)) return(NULL)
  as.numeric(unlist(strsplit(as.character(val), "\\|")))
}

get_vcov <- function(data) {
  val <- data$Temp_Vcov_Matrix
  if(length(val) == 0 || is.na(val)) return(NULL)
  vec <- as.numeric(unlist(strsplit(as.character(val), "\\|")))
  n <- sqrt(length(vec)) 
  matrix(vec, nrow=n, ncol=n)
}

# 5. META-ANALYSE FÜR PRE UND POST GETRENNT DURCHFÜHREN

perioden <- c("Pre", "Post")

cat("\n=== HITZE-EFFEKTE (99. PERZENTIL) NACH PERIODE ===\n")

for (p in perioden) {
  df_sub <- df_period %>% filter(Period == p)
  ags_list <- unique(df_sub$AGS_Zahlen)
  
  coef_list <- list()
  vcov_list <- list()
  
  for (ags in ags_list) {
    fdz_ags <- df_sub %>% filter(AGS_Zahlen == ags)
    if(nrow(fdz_ags) == 0) next
    
    c1 <- get_coef(fdz_ags)
    v1 <- get_vcov(fdz_ags)
    
    if(!is.null(c1) && !is.null(v1)) {
      coef_list[[ags]] <- c1
      vcov_list[[ags]] <- v1
    }
  }
  
  # Poolen
  coef_mat <- do.call(rbind, coef_list)
  meta_model <- mixmeta(coef_mat ~ 1, S = vcov_list, method = "reml")
  
  # Effekt für das 99. Perzentil berechnen
  cp <- crosspred(bvar_pool, 
                  coef = coef(meta_model), 
                  vcov = vcov(meta_model), 
                  model.link = "log", 
                  at = perc99_temp_all, 
                  cen = perc75_temp_all)
  
  rr_99 <- round(cp$allRRfit[1], 3)
  ci_low <- round(cp$allRRlow[1], 3)
  ci_up <- round(cp$allRRhigh[1], 3)
  
  cat(sprintf("Periode %s: RR = %.3f (95%% CI: %.3f - %.3f)\n", p, rr_99, ci_low, ci_up))
}
cat("====================================================\n")
