# ==============================================================================
# Script: 05_12_fit_combined_first_stage_and_sensitivity_models.R
#
# Purpose:
#   Fits district-specific quasi-Poisson time-series models for the outcomes
#   specified in the paired configuration script. The models estimate
#   cumulative temperature associations, pre-implementation alert-day
#   associations, Difference-in-Differences interaction estimates, and
#   lag-specific alert-day associations. Additional models assess sensitivity
#   to alternative degrees of freedom for the temperature function.
#
# Inputs:
#   - District-specific daily time-series CSV files containing hospital
#     admissions, meteorological exposures, calendar variables, reconstructed
#     pre-implementation alert days, and official post-implementation warnings
#
# Outputs:
#   - District-specific primary first-stage estimates
#   - Degrees-of-freedom sensitivity estimates
#   - District-exclusion file
#   - Model error log
#
# Configured outcomes:
#   Outcomes listed in diagnoses_to_analyze in
#   05_11_configure_combined_first_stage_and_sensitivity_models.R
#
# Primary analysis period:
#   May–September, 2000–2009
#
# Dependencies:
#   This script is intended to be sourced by
#   05_11_configure_combined_first_stage_and_sensitivity_models.R. Required
#   packages, paths, outcomes, alert variables, analysis periods, and
#   sensitivity parameters are defined in that configuration script.
#
# Data access:
#   Hospital data and district-specific input files are restricted and are not
#   included in the repository.
#
# Repository note:
#   The file name and explanatory comments were standardized.
# ==============================================================================

# =====================================================================
# 1. SETUP
# =====================================================================

# CSV-Dateien im Input-Ordner finden
ags_files <- list.files(path = input_path, pattern = "^AGS_.*\\.csv$", full.names = TRUE)
if(length(ags_files) == 0) {
  stop("FEHLER: Keine AGS-CSVs im input_path gefunden!")
}
print(paste(length(ags_files), "Landkreis-Dateien für die Analyse gefunden."))

# Fehlerprotokoll initialisieren
error_log <- data.frame(Diagnosis = character(), AGS = character(), RFM = character(), Error_Message = character(), stringsAsFactors = FALSE)

# Ankerpunkte für Model 2 und Model 3 (Gasparrini)
day1 <- as.numeric(as.Date("2002-07-15"))
day2 <- as.numeric(as.Date("2007-07-15")) 

QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

results_list <- list()
sens_results_list <- list() # Neue Liste für den Sensitivitätstest
dropped_ags_list <- list()

# Hilfsfunktion, um Lags sauber auszulesen
get_lag_val <- function(cp, lag_idx, stat = "RRfit") {
  if (stat == "RRfit") return(cp$matRRfit[1, lag_idx])
  if (stat == "RRlow") return(cp$matRRlow[1, lag_idx])
  if (stat == "RRhigh") return(cp$matRRhigh[1, lag_idx])
}

# =====================================================================
# 2. SCHLEIFE (DIAGNOSEN -> RFM -> AGS)
# =====================================================================

for (current_diag in diagnoses_to_analyze) {
  
  for (current_rfm in rfm_list) {
    
    print(paste(">> Berechne Modelle für:", current_rfm, "| Diagnose:", current_diag))
    counter <- 0
    
    for (file_path in ags_files) {
      
      counter <- counter + 1
      if(counter %% 50 == 0) cat(paste0("...", counter, "\n")) 
      
      current_ags_name <- basename(file_path) 
      
      # Daten laden 
      df <- fread(file_path)
      
      if ("Datum" %in% names(df)) setnames(df, "Datum", "date")
      
      # Check
      if(!current_rfm %in% names(df) || !current_diag %in% names(df)) next
      
      # Datums-Korrektur vor Aggregation
      df[, date := as.Date(date)]
      
      # Aggregation auf die SPEZIFISCHE Diagnose pro Tag
      agsdat <- df[, .(Outcome_Count = sum(get(current_diag), na.rm = TRUE)), 
                   by = .(date, T_mean_popw, holiday, get(current_rfm), heat_alerts_dwd)]
      
      setnames(agsdat, "get", current_rfm) 
      
      # Variablen für Modell aufbauen
      agsdat[, year    := as.numeric(format(date, "%Y"))]
      agsdat[, month   := as.numeric(format(date, "%m"))]
      agsdat[, dow     := as.factor(format(date, "%u"))]
      agsdat[, yday.sm := as.numeric(format(date, "%j"))]
      
      # HW_Status (Hitzewarnsystem-Implementierung) 
      agsdat[, HW_Status := ifelse(year >= 2005, 1, 0)]
      
      # Zeitraum filtern (2000 - 2009)
      agsdat <- agsdat[year %in% years_analysis & month %in% months_analysis]
      setorder(agsdat, date)
      
      # combined_alerts bauen (vor 2005 = RFM, ab 2005 = echter DWD Alert)
      agsdat[, combined_alerts := ifelse(year < 2005, get(current_rfm), heat_alerts_dwd)]
      
      # SICHERHEITSCHECKS
      total_cases      <- sum(agsdat$Outcome_Count, na.rm=TRUE)
      total_cases_pre  <- sum(agsdat$Outcome_Count[agsdat$year < 2005], na.rm=TRUE)
      total_cases_post <- sum(agsdat$Outcome_Count[agsdat$year >= 2005], na.rm=TRUE)
      n_days           <- nrow(agsdat)
      
      if(n_days < 100 || total_cases < 10) {
        dropped_ags_list[[length(dropped_ags_list) + 1]] <- data.frame(
          AGS_File = current_ags_name, RFM = current_rfm, Diagnosis = current_diag, 
          Total_Cases = total_cases, Total_Days = n_days
        )
        next 
      }
      
      if(sum(is.na(agsdat$T_mean_popw)) > 0) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Abbruch: NAs in Temperatur gefunden"))
        next
      }
      
      if(length(unique(agsdat$year)) < 2) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Abbruch: Kontrast-Fehler"))
        next
      }
      
      warnings_pre_2005  <- sum(agsdat$combined_alerts[agsdat$year < 2005], na.rm = TRUE)
      warnings_post_2005 <- sum(agsdat$combined_alerts[agsdat$year >= 2005], na.rm = TRUE)
      
      if(warnings_pre_2005 == 0 || warnings_post_2005 == 0) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = "Abbruch: Warnungen fehlen"))
        next
      }
      
      outcome_data  <- agsdat$Outcome_Count
      temp_data     <- agsdat$T_mean_popw
      eligible_data <- agsdat$combined_alerts
      holiday_data  <- agsdat$holiday
      impl_data     <- agsdat$HW_Status
      
      # =====================================================================
      # 3. MODELLIERUNG (DLNM & GLM)
      # =====================================================================
      
      tryCatch({
        
        cbtmean <- crossbasis(temp_data, lag=lag, 
                              argvar=list(fun=varfun, degree=vardegree, knots=quantile(temp_data, varper/100, na.rm=TRUE)), 
                              arglag=list(fun="integer"), group=agsdat$year)
        
        datenum  <- as.numeric(agsdat$date)
        intfirst <- ((datenum-day1)/(day2-day1))*cbtmean
        intlast  <- ((datenum-day2)/(day1-day2))*cbtmean 
        
        cbalert <- crossbasis(eligible_data, lag=lag,
                              argvar=list(fun="lin"),
                              arglag=list(fun="integer"), group=agsdat$year)
        
        alertint <- impl_data * cbalert
        
        trend_df <- round(length(unique(agsdat$year))/dftrend/10)
        if(trend_df < 1) trend_df <- 1 
        
        m0 <- glm(outcome_data ~ cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=trend_df) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        m1 <- glm(outcome_data ~ cbtmean + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=trend_df) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        m2 <- glm(outcome_data ~ cbtmean + intfirst + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=trend_df) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        m3 <- glm(outcome_data ~ cbtmean + intlast + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=trend_df) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
        
        # EFFEKTE EXTRAHIEREN 
        mean_temp <- mean(temp_data, na.rm=TRUE) 
        
        red_temp_1 <- crossreduce(cbtmean, m1, cen=mean_temp)
        red_temp_2 <- crossreduce(cbtmean, m2, cen=mean_temp)
        red_temp_3 <- crossreduce(cbtmean, m3, cen=mean_temp)
        
        coef_temp_1 <- paste(coef(red_temp_1), collapse="|")
        coef_temp_2 <- paste(coef(red_temp_2), collapse="|")
        coef_temp_3 <- paste(coef(red_temp_3), collapse="|")
        
        vcov_temp_1 <- paste(as.vector(vcov(red_temp_1)), collapse="|")
        vcov_temp_2 <- paste(as.vector(vcov(red_temp_2)), collapse="|")
        vcov_temp_3 <- paste(as.vector(vcov(red_temp_3)), collapse="|")
        
        # Für den Gesamt-Effekt
        p_alert_0 <- crossreduce(cbalert, m0, at=1)
        p_alert_1 <- crossreduce(cbalert, m1, at=1)
        p_alert_2 <- crossreduce(cbalert, m2, at=1)
        p_alert_3 <- crossreduce(cbalert, m3, at=1)
        
        p_int_0 <- crossreduce(alertint, m0, at=1)
        p_int_1 <- crossreduce(alertint, m1, at=1)
        p_int_2 <- crossreduce(alertint, m2, at=1)
        p_int_3 <- crossreduce(alertint, m3, at=1)
        
        # NEU: Lag-spezifische RRs für Alert extrahieren 
        cp_alert_0 <- crosspred(cbalert, m0, at=1)
        cp_alert_1 <- crosspred(cbalert, m1, at=1)
        cp_alert_2 <- crosspred(cbalert, m2, at=1)
        cp_alert_3 <- crosspred(cbalert, m3, at=1)
        
        # ZAHLEN SPEICHERN
        res_tmp <- data.frame(
          AGS_File = current_ags_name,
          RFM = current_rfm,
          Diagnosis = current_diag, 
          Model = c("model0", "model1", "model2", "model3"),
          Total_Cases = total_cases,                                       
          Total_Cases_Pre_2005 = total_cases_pre,  
          Total_Cases_Post_2005 = total_cases_post,
          QAIC = c(QAIC(m0), QAIC(m1), QAIC(m2), QAIC(m3)), 
          
          Temp_Mean_Cen = mean_temp, 
          Temp_Coef_Vector = c(NA, coef_temp_1, coef_temp_2, coef_temp_3),
          Temp_Vcov_Matrix = c(NA, vcov_temp_1, vcov_temp_2, vcov_temp_3),
          
          RR_Alert = c(p_alert_0$RRfit, p_alert_1$RRfit, p_alert_2$RRfit, p_alert_3$RRfit),                
          CI_Low_Alert = c(p_alert_0$RRlow, p_alert_1$RRlow, p_alert_2$RRlow, p_alert_3$RRlow),   
          CI_Up_Alert = c(p_alert_0$RRhigh, p_alert_1$RRhigh, p_alert_2$RRhigh, p_alert_3$RRhigh), 
          
          # NEU: Lags separat speichern
          RR_Alert_Lag0 = c(get_lag_val(cp_alert_0, 1), get_lag_val(cp_alert_1, 1), get_lag_val(cp_alert_2, 1), get_lag_val(cp_alert_3, 1)),
          CI_Low_Alert_Lag0 = c(get_lag_val(cp_alert_0, 1, "RRlow"), get_lag_val(cp_alert_1, 1, "RRlow"), get_lag_val(cp_alert_2, 1, "RRlow"), get_lag_val(cp_alert_3, 1, "RRlow")),
          CI_Up_Alert_Lag0 = c(get_lag_val(cp_alert_0, 1, "RRhigh"), get_lag_val(cp_alert_1, 1, "RRhigh"), get_lag_val(cp_alert_2, 1, "RRhigh"), get_lag_val(cp_alert_3, 1, "RRhigh")),
          
          RR_Alert_Lag1 = c(get_lag_val(cp_alert_0, 2), get_lag_val(cp_alert_1, 2), get_lag_val(cp_alert_2, 2), get_lag_val(cp_alert_3, 2)),
          RR_Alert_Lag2 = c(get_lag_val(cp_alert_0, 3), get_lag_val(cp_alert_1, 3), get_lag_val(cp_alert_2, 3), get_lag_val(cp_alert_3, 3)),
          RR_Alert_Lag3 = c(get_lag_val(cp_alert_0, 4), get_lag_val(cp_alert_1, 4), get_lag_val(cp_alert_2, 4), get_lag_val(cp_alert_3, 4)),
          
          DiD_Coef = c(coef(p_int_0), coef(p_int_1), coef(p_int_2), coef(p_int_3)),                            
          DiD_Variance = c(vcov(p_int_0), vcov(p_int_1), vcov(p_int_2), vcov(p_int_3)),                        
          RR_DiD = c(p_int_0$RRfit, p_int_1$RRfit, p_int_2$RRfit, p_int_3$RRfit),                
          CI_Low_DiD = c(p_int_0$RRlow, p_int_1$RRlow, p_int_2$RRlow, p_int_3$RRlow),   
          CI_Up_DiD = c(p_int_0$RRhigh, p_int_1$RRhigh, p_int_2$RRhigh, p_int_3$RRhigh)
        )
        
        results_list[[length(results_list) + 1]] <- res_tmp
        
        # =====================================================================
        # SENSITIVITÄTSTEST FÜR OVER-ADJUSTMENT (Nur für Model 1)
        # =====================================================================
        
        for(sens_df in df_sens_list) {
          
          cbtmean_sens <- crossbasis(temp_data, lag=lag, 
                                     argvar=list(fun=varfun, degree=vardegree, df=sens_df), 
                                     arglag=list(fun="integer"), group=agsdat$year)
          
          m1_sens <- glm(outcome_data ~ cbtmean_sens + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=trend_df) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
          
          p_alert_sens <- crossreduce(cbalert, m1_sens, at=1)
          p_int_sens   <- crossreduce(alertint, m1_sens, at=1)
          
          # Auch hier die Lags per crosspred ziehen, um absolut sicher zu sein
          cp_alert_sens <- crosspred(cbalert, m1_sens, at=1)
          
          sens_tmp <- data.frame(
            AGS_File = current_ags_name,
            RFM = current_rfm,
            Diagnosis = current_diag,
            Tested_DF = sens_df,
            
            RR_Alert = p_alert_sens$RRfit,
            CI_Low_Alert = p_alert_sens$RRlow,
            CI_Up_Alert = p_alert_sens$RRhigh,
            
            # NEU: Lag-spezifische RRs für das Sensitivitätsmodell
            RR_Alert_Lag0 = get_lag_val(cp_alert_sens, 1),
            RR_Alert_Lag1 = get_lag_val(cp_alert_sens, 2),
            RR_Alert_Lag2 = get_lag_val(cp_alert_sens, 3),
            RR_Alert_Lag3 = get_lag_val(cp_alert_sens, 4),
            
            RR_DiD = p_int_sens$RRfit,
            CI_Low_DiD = p_int_sens$RRlow,
            CI_Up_DiD = p_int_sens$RRhigh
          )
          
          sens_results_list[[length(sens_results_list) + 1]] <- sens_tmp
        }
        
      }, error = function(e) {
        error_log <<- rbind(error_log, data.frame(Diagnosis = current_diag, AGS = current_ags_name, RFM = current_rfm, Error_Message = e$message))
      }) # Ende tryCatch   
    } # Ende AGS Loop
  } # Ende RFM Loop
} # Ende Diagnose Loop

# =====================================================================
# 4. EXPORT DER ERGEBNISSE
# =====================================================================

prefix <- gsub(".R", "", syntax_name)

if (length(results_list) > 0) {
  final_df <- rbindlist(results_list)
  out_name <- paste0(prefix, "_Results.csv")
  write.csv(final_df, file.path(output_path, out_name), row.names = FALSE)
  print(paste("Haupt-Ergebnisse gespeichert unter:", out_name))
}

# Export Sensitivitätstest
if (length(sens_results_list) > 0) {
  sens_df <- rbindlist(sens_results_list)
  sens_out_name <- paste0(prefix, "_Sensitivity_DF.csv")
  write.csv(sens_df, file.path(output_path, sens_out_name), row.names = FALSE)
  print(paste("Sensitivitäts-Ergebnisse gespeichert unter:", sens_out_name))
}

if (length(dropped_ags_list) > 0) {
  dropped_df <- rbindlist(dropped_ags_list)
  dropped_df <- unique(dropped_df, by = c("Diagnosis", "AGS_File", "RFM")) 
  dropped_out_name <- paste0(prefix, "_Dropped_AGS.csv")
  write.csv(dropped_df, file.path(output_path, dropped_out_name), row.names = FALSE)
  print(paste("Liste der ausgeschlossenen Kreise liegt hier:", dropped_out_name))
} else {
  print("Kein Landkreis musste ausgeschlossen werden.")
}

if (nrow(error_log) > 0) {
  error_out_name <- paste0(prefix, "_Error_Log.csv")
  write.csv(error_log, file.path(output_path, error_out_name), row.names = FALSE)
  print(paste("Abstürze vorhanden, Fehlerprotokoll liegt hier:", error_out_name))
} else {
  print("Keine Modell-Abstürze")
}
