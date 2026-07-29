# ==============================================================================
# Script: 06_07_06_pool_post_period_reference_estimates.R
#
# Purpose:
#   Pools the results of the reverse-coded reference-period analysis. The alert
#   main effect represents the post-implementation alert-day association, while
#   the reverse-coded interaction represents the pre-implementation versus
#   post-implementation contrast using the post period as the reference.
#
# Inputs:
#   - Released reference-switch first-stage FDZ output table
#
# Outputs:
#   - Pooled post-implementation alert-day association printed to the console
#   - Pooled reverse-coded Difference-in-Differences estimate printed to the
#     console
#   - Heterogeneity statistics printed to the console
#
# Outcome:
#   All-cause hospital admissions
#
# Model specification:
#   RFM5 and Model 1
#
# Required packages:
#   mixmeta, dplyr, readxl
#
# Repository note:
#   The file name, path strings, and explanatory comments were standardized. 
# ==============================================================================

library(mixmeta)
library(dplyr)
library(readxl)

# =====================================================================
# 1. DATEN EINLESEN
# =====================================================================

df_switch <- read_excel(
  "data/processed/fdz_outputs/Results_FirstStage_ReferenceSwitch_gh.xlsx"
) %>%
  # Filter auf RFM5, All Admissions und Hauptmodell 1
  filter(RFM == "heat_alerts_rfm5" & Diagnosis == "AllAdmissions" & Model == "model1") %>%
  
  mutate(
    RR_Alert_Post2005 = as.numeric(gsub(",", ".", as.character(RR_Alert_Post2005))),
    CI_Low_Alert_Post2005 = as.numeric(gsub(",", ".", as.character(CI_Low_Alert_Post2005))),
    CI_Up_Alert_Post2005 = as.numeric(gsub(",", ".", as.character(CI_Up_Alert_Post2005))),
    
    RR_DiD_Rev = as.numeric(gsub(",", ".", as.character(RR_DiD_Rev))),
    CI_Low_DiD_Rev = as.numeric(gsub(",", ".", as.character(CI_Low_DiD_Rev))),
    CI_Up_DiD_Rev = as.numeric(gsub(",", ".", as.character(CI_Up_DiD_Rev)))
  ) %>%
  arrange(AGS_File)

# =====================================================================
# 2. META-ANALYSE: DER ABSOLUTE ALERT-EFFEKT IN DER POST-PERIODE
# =====================================================================

yi_alert <- log(df_switch$RR_Alert_Post2005)
sei_alert <- (log(df_switch$CI_Up_Alert_Post2005) - log(df_switch$CI_Low_Alert_Post2005)) / (2 * 1.96)
vi_alert <- sei_alert^2

valid_idx_alert <- which(!is.na(yi_alert) & !is.na(vi_alert) & vi_alert > 0)

meta_alert <- mixmeta(yi_alert[valid_idx_alert] ~ 1, S = vi_alert[valid_idx_alert], method = "reml")

RRpooled_alert <- exp(predict(meta_alert, ci=TRUE))
summ_alert <- summary(meta_alert)

alert_rr <- round(RRpooled_alert[1, "fit"], 3)
alert_ci_low <- round(RRpooled_alert[1, "ci.lb"], 3)
alert_ci_up <- round(RRpooled_alert[1, "ci.ub"], 3)
alert_i2 <- round(max(summ_alert$i2stat), 1)
alert_q <- round(summ_alert$qstat$Q, 1)

# =====================================================================
# 3. META-ANALYSE: DER UMGEKEHRTE DiD-EFFEKT
# =====================================================================

yi_did <- log(df_switch$RR_DiD_Rev)
sei_did <- (log(df_switch$CI_Up_DiD_Rev) - log(df_switch$CI_Low_DiD_Rev)) / (2 * 1.96)
vi_did <- sei_did^2

valid_idx_did <- which(!is.na(yi_did) & !is.na(vi_did) & vi_did > 0)

meta_did <- mixmeta(yi_did[valid_idx_did] ~ 1, S = vi_did[valid_idx_did], method = "reml")

RRpooled_did <- exp(predict(meta_did, ci=TRUE))
summ_did <- summary(meta_did)

did_rr <- round(RRpooled_did[1, "fit"], 3)
did_ci_low <- round(RRpooled_did[1, "ci.lb"], 3)
did_ci_up <- round(RRpooled_did[1, "ci.ub"], 3)
did_i2 <- round(max(summ_did$i2stat), 1)
did_q <- round(summ_did$qstat$Q, 1)

# =====================================================================
# 4. OUTPUT FÜR DEN TEXT IN WORD
# =====================================================================

cat("\n=== ZAHLEN FÜR DEN TEXT: REFERENCE SWITCH (ALL CAUSE) ===\n")
cat("\n--- 1. Absolutes RR an Warntagen (Post-2005) ---\n")
cat("Gepooltes RR:      ", alert_rr, "\n")
cat("95% CI Lower:      ", alert_ci_low, "\n")
cat("95% CI Upper:      ", alert_ci_up, "\n")
cat("I^2 (Heterog.):    ", alert_i2, "%\n")
cat("Cochran's Q:       ", alert_q, "\n")

cat("\n--- 2. Umgekehrter DiD-Effekt (Reference: Post) ---\n")
cat("Gepooltes RR:      ", did_rr, "\n")
cat("95% CI Lower:      ", did_ci_low, "\n")
cat("95% CI Upper:      ", did_ci_up, "\n")
cat("I^2 (Heterog.):    ", did_i2, "%\n")
cat("Cochran's Q:       ", did_q, "\n")
cat("==========================================================\n")
