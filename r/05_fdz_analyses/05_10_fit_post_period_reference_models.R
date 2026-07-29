# ==============================================================================
# Script: 05_10_fit_post_period_reference_models.R
#
# Purpose:
#   Fits district-specific quasi-Poisson time-series models using a reverse-coded
#   implementation-period indicator. The alert main effect directly represents
#   the post-implementation alert-day association, while the interaction
#   represents the pre-implementation to post-implementation contrast in the
#   reverse direction.
#
# Inputs:
#   - District-specific daily time-series CSV files containing hospital
#     admissions, meteorological exposures, calendar variables, reconstructed
#     pre-implementation alert days, and official post-implementation alerts
#
# Outputs:
#   - Results_FirstStage_ReferenceSwitch.csv
#   - Dropped_AGS_ReferenceSwitch.csv, where applicable
#   - Error_Log_ReferenceSwitch.csv, where applicable
#
# Analysis period:
#   May–September, 2000–2009
#
# Dependencies:
#   This script is sourced by
#   05_09_configure_post_period_reference_models.R. Required packages, paths,
#   outcomes, alert variables, and model parameters are defined in that
#   configuration script.
#
# Data access:
#   Hospital data and district-specific input files are restricted and are not
#   included in the repository.
#
# Repository note:
#   The file name and explanatory comments were standardized.
# ==============================================================================

# =====================================================================
# VORBEREITUNG
# =====================================================================

ags_files <- list.files(path = input_path, pattern = "^AGS_.*\\.csv$", full.names = TRUE)

if(length(ags_files) == 0) stop("FEHLER: Keine AGS-CSVs gefunden!")
print(paste(length(ags_files), "Landkreis-Dateien für die Analyse gefunden."))

error_log <- data.frame(Diagnosis = character(), AGS = character(), RFM = character(), Error_Message = character(), stringsAsFactors = FALSE)

# Days for Centering 
day1 <- as.numeric(as.Date("2000-07-15"))
day2 <- as.numeric(as.Date("2009-07-15"))

QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

results_list <- list()
dropped_ags_list <- list()

# =====================================================================
# HAUPT-SCHLEIFEN
# =====================================================================

for (current_diag in diagnoses_to_analyze) {
  for (current_rfm in rfm_list) {
    
    print(paste(">> Berechne Modelle (Reference Switch) für:", current_diag, "|", current_rfm))
    counter <- 0
    
    for (file_path in ags_files) {
      
      counter <- counter + 1
      if(counter %% 50 == 0) cat(paste0("...", counter, "\n")) 
      
      current_ags_name <- basename(file_path) 
      df <- fread(file_path)
      if ("Datum" %in% names(df)) setnames(df, "Datum", "date")
      
      # Dynamische AllAdmissions-Berechnung
      
      if (current_diag == "AllAdmissions") {
        valid_cols <- intersect(vars_to_aggregate, names(df))
        df[, AllAdmissions := rowSums(.SD, na.rm = TRUE), .SDcols = valid_cols]
      }
      
      if(!current_rfm %in% names(df) || !current_diag %in% names(df)) next
      
      df[, date := as.Date(date)]
      
      # Aggregation auf Diagnose
      
      agsdat <- df[, .(Outcome_Count = sum(get(current_diag), na.rm = TRUE)), 
                   by = .(date, T_mean_popw, holiday, get(current_rfm), heat_alerts_dwd)]
      
      setnames(agsdat, "get", current_rfm) 
      
      agsdat[, year  := as.numeric(format(date, "%Y"))]
      agsdat[, month := as.numeric(format(date, "%m"))]
      agsdat[, dow   := as.factor(format(date, "%u"))]
      agsdat[, yday.sm := as.numeric(format(date, "%j"))]
      
      # combined_alerts bauen (Pre = RFM, Post = DWD)
      agsdat[, combined_alerts := ifelse(year < 2005, get(current_rfm), heat_alerts_dwd)]
      
      # =====================================================================
      # REFERENCE SWITCH
      # =====================================================================
      
      # Bisher: ifelse(year >= 2005, 1, 0)
      # NEU: Vor 2005 = 1, Ab 2005 = 0. -> Post-2005 wird zur Referenzkategorie
      agsdat[, HW_Status_rev := ifelse(year < 2005, 1, 0)]
      
      # Zeitraum filtern 
      agsdat <- agsdat[year %in% years_analysis & month %in% months_analysis]
      setorder(agsdat, date)
      
      total_cases <- sum(agsdat$Outcome_Count, na.rm=TRUE)
      n_days <- nrow(agsdat)
      
      if(n_days < 100 || total_cases < 10) {
        dropped_ags_list[[length(dropped_ags_list) + 1]] <- data.frame(
          AGS_File = current_ags_name, RFM = current_rfm, Diagnosis = current_diag,
          Total_Cases = total_cases, Total_Days = n_days
        )
        next 
      }
      
      if(sum(is.na(agsdat$T_mean_popw)) > 0) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = "NAs in Temperatur"))
        next
      }
      
      if(length(unique(agsdat$year)) < 2) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Kontrast-Fehler"))
        next
      }
      
      warnings_pre_2005 <- sum(agsdat$combined_alerts[agsdat$year < 2005], na.rm = TRUE)
      warnings_post_2005 <- sum(agsdat$combined_alerts[agsdat$year >= 2005], na.rm = TRUE)
      
      if(warnings_pre_2005 == 0 || warnings_post_2005 == 0) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Warnungen fehlen in einer der Perioden"))
        next
      }
      
      outcome_data  <- agsdat$Outcome_Count
      temp_data     <- agsdat$T_mean_popw
      eligible_data <- agsdat$combined_alerts
      holiday_data  <- agsdat$holiday
      impl_data_rev <- agsdat$HW_Status_rev 
      
      tryCatch({
        
        cbtmean <- crossbasis(temp_data, lag=lag, 
                              argvar=list(fun=varfun, degree=vardegree, knots=quantile(temp_data, varper/100, na.rm=TRUE)), 
                              arglag=list(fun="integer"), group=agsdat$year)
        
        datenum <- as.numeric(agsdat$date)
        intfirst <- ((datenum-day1)/(day2-day1))*cbtmean
        intlast  <- ((datenum-day2)/(day1-day2))*cbtmean 
        
        cbalert <- crossbasis(eligible_data, lag=lag,
                              argvar=list(fun="lin"),
                              arglag=list(fun="integer"), group=agsdat$year)
        
        # DiD-Term (Referenz ist jetzt Post-2005)
        
        alertint_rev <- impl_data_rev * cbalert
        
        m0 <- glm(outcome_data ~ cbalert + alertint_rev + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        m1 <- glm(outcome_data ~ cbtmean + cbalert + alertint_rev + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        m2 <- glm(outcome_data ~ cbtmean + intfirst + cbalert + alertint_rev + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        m3 <- glm(outcome_data ~ cbtmean + intlast + cbalert + alertint_rev + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        
        mean_temp <- mean(temp_data, na.rm=TRUE) 
        
        red_temp_1 <- crossreduce(cbtmean, m1, cen=mean_temp)
        red_temp_2 <- crossreduce(cbtmean, m2, cen=mean_temp)
        red_temp_3 <- crossreduce(cbtmean, m3, cen=mean_temp)
        
        # The alert main effect represents the post-implementation alert association
        
        p_alert_0 <- crossreduce(cbalert, m0, at=1)
        p_alert_1 <- crossreduce(cbalert, m1, at=1)
        p_alert_2 <- crossreduce(cbalert, m2, at=1)
        p_alert_3 <- crossreduce(cbalert, m3, at=1)
        
        # p_int zeigt jetzt den Unterschied "Pre minus Post" 
        
        p_int_0 <- crossreduce(alertint_rev, m0, at=1)
        p_int_1 <- crossreduce(alertint_rev, m1, at=1)
        p_int_2 <- crossreduce(alertint_rev, m2, at=1)
        p_int_3 <- crossreduce(alertint_rev, m3, at=1)
        
        res_tmp <- data.frame(
          AGS_File = current_ags_name,
          RFM = current_rfm,
          Diagnosis = current_diag, 
          Model = c("model0", "model1", "model2", "model3"),
          Total_Cases = total_cases,      
          QAIC = c(QAIC(m0), QAIC(m1), QAIC(m2), QAIC(m3)), 
          
          Temp_Mean_Cen = mean_temp, 
          Temp_Coef_Vector = c(NA, paste(coef(red_temp_1), collapse="|"), paste(coef(red_temp_2), collapse="|"), paste(coef(red_temp_3), collapse="|")),
          Temp_Vcov_Matrix = c(NA, paste(as.vector(vcov(red_temp_1)), collapse="|"), paste(as.vector(vcov(red_temp_2)), collapse="|"), paste(as.vector(vcov(red_temp_3)), collapse="|")),
          
          # Store the direct post-implementation alert association
          Alert_Coef = c(coef(p_alert_0), coef(p_alert_1), coef(p_alert_2), coef(p_alert_3)),
          Alert_Variance = c(vcov(p_alert_0), vcov(p_alert_1), vcov(p_alert_2), vcov(p_alert_3)),
          RR_Alert_Post2005 = c(p_alert_0$RRfit, p_alert_1$RRfit, p_alert_2$RRfit, p_alert_3$RRfit),                
          CI_Low_Alert_Post2005 = c(p_alert_0$RRlow, p_alert_1$RRlow, p_alert_2$RRlow, p_alert_3$RRlow),   
          CI_Up_Alert_Post2005 = c(p_alert_0$RRhigh, p_alert_1$RRhigh, p_alert_2$RRhigh, p_alert_3$RRhigh),
          
          # umgedrehter DiD-Term
          DiD_Rev_Coef = c(coef(p_int_0), coef(p_int_1), coef(p_int_2), coef(p_int_3)),
          DiD_Rev_Variance = c(vcov(p_int_0), vcov(p_int_1), vcov(p_int_2), vcov(p_int_3)),
          RR_DiD_Rev = c(p_int_0$RRfit, p_int_1$RRfit, p_int_2$RRfit, p_int_3$RRfit),                
          CI_Low_DiD_Rev = c(p_int_0$RRlow, p_int_1$RRlow, p_int_2$RRlow, p_int_3$RRlow),   
          CI_Up_DiD_Rev = c(p_int_0$RRhigh, p_int_1$RRhigh, p_int_2$RRhigh, p_int_3$RRhigh)
        )
        
        results_list[[length(results_list) + 1]] <- res_tmp
        
      }, error = function(e) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = e$message))
      }) 
    } 
  } 
} 

# =====================================================================
# EXPORT
# =====================================================================

if (length(results_list) > 0) {
  final_df <- rbindlist(results_list)
  out_name <- "Results_FirstStage_ReferenceSwitch.csv"
  write.csv(final_df, file.path(output_path, out_name), row.names = FALSE)
  print(paste("Ergebnisse gespeichert unter:", out_name))
}

if (length(dropped_ags_list) > 0) {
  dropped_df <- unique(rbindlist(dropped_ags_list), by = c("Diagnosis", "AGS_File", "RFM")) 
  write.csv(dropped_df, file.path(output_path, "Dropped_AGS_ReferenceSwitch.csv"), row.names = FALSE)
}

if (nrow(error_log) > 0) {
  write.csv(error_log, file.path(output_path, "Error_Log_ReferenceSwitch.csv"), row.names = FALSE)
}
