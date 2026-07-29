# ==============================================================================
# Script: 05_08_fit_period_specific_first_stage_models.R
#
# Purpose:
#   Fits district-specific quasi-Poisson time-series models separately for the
#   pre-implementation period and the post-implementation period. Temperature
#   and alert-day associations are estimated within each period without a
#   Difference-in-Differences interaction.
#
# Inputs:
#   - District-specific daily time-series CSV files containing hospital
#     admissions, meteorological exposures, calendar variables, reconstructed
#     pre-implementation alert days, and official post-implementation alerts
#
# Outputs:
#   - Results_FirstStage_GetrenntePerioden.csv
#   - Dropped_AGS_GetrenntePerioden.csv, where applicable
#   - Error_Log_GetrenntePerioden.csv, where applicable
#
# Analysis periods:
#   Pre-implementation:  May–September, 2000–2004
#   Post-implementation: May–September, 2005–2009
#
# Dependencies:
#   This script is sourced by
#   05_07_configure_period_specific_first_stage_models.R. Required packages,
#   paths, outcomes, alert variables, and model parameters are defined in that
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

error_log <- data.frame(Diagnosis = character(), Period = character(), AGS = character(), RFM = character(), Error_Message = character(), stringsAsFactors = FALSE)

QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

results_list <- list()
dropped_ags_list <- list()

periods_to_analyze <- list(
  "Pre"  = 2000:2004,
  "Post" = 2005:2009
)

# =====================================================================
# HAUPT-SCHLEIFEN
# =====================================================================

for (current_diag in diagnoses_to_analyze) {
  for (current_rfm in rfm_list) {
    for (period_name in names(periods_to_analyze)) {
      
      current_years <- periods_to_analyze[[period_name]]
      print(paste(">> Berechne:", current_diag, "|", current_rfm, "| Periode:", period_name))
      counter <- 0
      
      if(period_name == "Pre") {
        day1 <- as.numeric(as.Date("2000-07-15"))
        day2 <- as.numeric(as.Date("2004-07-15"))
      } else {
        day1 <- as.numeric(as.Date("2005-07-15"))
        day2 <- as.numeric(as.Date("2009-07-15"))
      }
      
      for (file_path in ags_files) {
        
        counter <- counter + 1
        if(counter %% 50 == 0) cat(paste0("...", counter, "\n")) 
        
        current_ags_name <- basename(file_path) 
        df <- fread(file_path)
        if ("Datum" %in% names(df)) setnames(df, "Datum", "date")
        
        # Dynamische Spaltenerstellung für AllAdmissions:
        if (current_diag == "AllAdmissions") {
          valid_cols <- intersect(vars_to_aggregate, names(df))
          df[, AllAdmissions := rowSums(.SD, na.rm = TRUE), .SDcols = valid_cols]
        }
        
        # Sicherheitscheck, ob die Spalte (jetzt) existiert
        if(!current_rfm %in% names(df) || !current_diag %in% names(df)) next
        
        df[, date := as.Date(date)]
        
        # Aggregation auf die spezifische Diagnose
        agsdat <- df[, .(Outcome_Count = sum(get(current_diag), na.rm = TRUE)), 
                     by = .(date, T_mean_popw, holiday, get(current_rfm), heat_alerts_dwd)]
        
        setnames(agsdat, "get", current_rfm) 
        
        agsdat[, year  := as.numeric(format(date, "%Y"))]
        agsdat[, month := as.numeric(format(date, "%m"))]
        agsdat[, dow   := as.factor(format(date, "%u"))]
        agsdat[, yday.sm := as.numeric(format(date, "%j"))]
        
        agsdat[, combined_alerts := ifelse(year < 2005, get(current_rfm), heat_alerts_dwd)]
        agsdat <- agsdat[year %in% current_years & month %in% months_analysis]
        setorder(agsdat, date)
        
        total_cases <- sum(agsdat$Outcome_Count, na.rm=TRUE)
        n_days <- nrow(agsdat)
        
        if(n_days < 100 || total_cases < 10) {
          dropped_ags_list[[length(dropped_ags_list) + 1]] <- data.frame(
            AGS_File = current_ags_name, RFM = current_rfm, Diagnosis = current_diag, Period = period_name,
            Total_Cases = total_cases, Total_Days = n_days
          )
          next 
        }
        
        if(sum(is.na(agsdat$T_mean_popw)) > 0) {
          error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, Period = period_name, AGS = current_ags_name, RFM = current_rfm, Error_Message = "NAs in Temperatur"))
          next
        }
        
        if(length(unique(agsdat$year)) < 2) {
          error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, Period = period_name, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Nur 1 Jahr Daten (Kontrast-Fehler)"))
          next
        }
        
        if(sum(agsdat$combined_alerts, na.rm = TRUE) == 0) {
          error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, Period = period_name, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Keine Warnungen in dieser Periode"))
          next
        }
        
        outcome_data  <- agsdat$Outcome_Count
        temp_data     <- agsdat$T_mean_popw
        eligible_data <- agsdat$combined_alerts
        holiday_data  <- agsdat$holiday
        
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
          
          m0 <- glm(outcome_data ~ cbalert + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
          m1 <- glm(outcome_data ~ cbtmean + cbalert + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
          m2 <- glm(outcome_data ~ cbtmean + intfirst + cbalert + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
          m3 <- glm(outcome_data ~ cbtmean + intlast + cbalert + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
          
          mean_temp <- mean(temp_data, na.rm=TRUE) 
          
          red_temp_1 <- crossreduce(cbtmean, m1, cen=mean_temp)
          red_temp_2 <- crossreduce(cbtmean, m2, cen=mean_temp)
          red_temp_3 <- crossreduce(cbtmean, m3, cen=mean_temp)
          
          p_alert_0 <- crossreduce(cbalert, m0, at=1)
          p_alert_1 <- crossreduce(cbalert, m1, at=1)
          p_alert_2 <- crossreduce(cbalert, m2, at=1)
          p_alert_3 <- crossreduce(cbalert, m3, at=1)
          
          res_tmp <- data.frame(
            AGS_File = current_ags_name,
            RFM = current_rfm,
            Diagnosis = current_diag, 
            Period = period_name,           
            Model = c("model0", "model1", "model2", "model3"),
            Total_Cases = total_cases,      
            QAIC = c(QAIC(m0), QAIC(m1), QAIC(m2), QAIC(m3)), 
            
            Temp_Mean_Cen = mean_temp, 
            Temp_Coef_Vector = c(NA, paste(coef(red_temp_1), collapse="|"), paste(coef(red_temp_2), collapse="|"), paste(coef(red_temp_3), collapse="|")),
            Temp_Vcov_Matrix = c(NA, paste(as.vector(vcov(red_temp_1)), collapse="|"), paste(as.vector(vcov(red_temp_2)), collapse="|"), paste(as.vector(vcov(red_temp_3)), collapse="|")),
            
            Alert_Coef = c(coef(p_alert_0), coef(p_alert_1), coef(p_alert_2), coef(p_alert_3)),
            Alert_Variance = c(vcov(p_alert_0), vcov(p_alert_1), vcov(p_alert_2), vcov(p_alert_3)),
            
            RR_Alert = c(p_alert_0$RRfit, p_alert_1$RRfit, p_alert_2$RRfit, p_alert_3$RRfit),                
            CI_Low_Alert = c(p_alert_0$RRlow, p_alert_1$RRlow, p_alert_2$RRlow, p_alert_3$RRlow),   
            CI_Up_Alert = c(p_alert_0$RRhigh, p_alert_1$RRhigh, p_alert_2$RRhigh, p_alert_3$RRhigh)
          )
          
          results_list[[length(results_list) + 1]] <- res_tmp
          
        }, error = function(e) {
          error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, Period = period_name, AGS = current_ags_name, RFM = current_rfm, Error_Message = e$message))
        }) 
      } 
    } 
  } 
} 

# =====================================================================
# EXPORT
# =====================================================================

if (length(results_list) > 0) {
  final_df <- rbindlist(results_list)
  out_name <- "Results_FirstStage_GetrenntePerioden.csv"
  write.csv(final_df, file.path(output_path, out_name), row.names = FALSE)
  print(paste("Ergebnisse gespeichert unter:", out_name))
}

if (length(dropped_ags_list) > 0) {
  dropped_df <- unique(rbindlist(dropped_ags_list), by = c("Diagnosis", "Period", "AGS_File", "RFM")) 
  write.csv(dropped_df, file.path(output_path, "Dropped_AGS_GetrenntePerioden.csv"), row.names = FALSE)
}

if (nrow(error_log) > 0) {
  write.csv(error_log, file.path(output_path, "Error_Log_GetrenntePerioden.csv"), row.names = FALSE)
}
