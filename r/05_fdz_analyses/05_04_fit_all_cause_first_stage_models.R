# ==============================================================================
# Script: 05_04_fit_all_cause_first_stage_models.R
#
# Purpose:
#   Fits district-specific quasi-Poisson time-series models for all-cause
#   hospital admissions. The analysis combines reconstructed Random Forest
#   alert days before 2005 with official DWD heat-alert days from 2005 onward
#   and estimates temperature, alert-day, and Difference-in-Differences
#   associations.
#
# Inputs:
#   - District-specific daily time-series CSV files containing hospital
#     admissions, meteorological exposures, calendar variables, reconstructed
#     heat-alert classifications, and official DWD heat alerts
#
# Outputs:
#   - Results_FirstStage_AllAdmissions.csv
#   - Dropped_AGS_AllAdmissions.csv, where applicable
#   - Error_Log_FirstStage_all_admissions.csv, where applicable
#
# Outcome:
#   All-cause hospital admissions derived by summing the configured diagnosis
#   groups
#
# Analysis period:
#   May–September, 2000–2009
#
# Dependencies:
#   This script is sourced by
#   05_03_configure_all_cause_first_stage_models.R. Required packages, paths,
#   analysis periods, outcome variables, alert variables, and model parameters
#   are defined in that configuration script.
#
# Data access:
#   Hospital data and district-specific input files are restricted and are not
#   included in the repository.
#
# Repository note:
#   The file name and explanatory comments were standardized.
# ==============================================================================

# CSV-Dateien im Input-Ordner finden
ags_files <- list.files(path = input_path, pattern = "^AGS_.*\\.csv$", full.names = TRUE)
if(length(ags_files) == 0) {
  stop("FEHLER: Keine AGS-CSVs im input_path gefunden!")
}
print(paste(length(ags_files), "Landkreis-Dateien für die Analyse gefunden."))

# Fehlerprotokoll einfügen
error_log <- data.frame(AGS = character(), Error_Message = character(), stringsAsFactors = FALSE)

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

# -------------------------------------------------------------------------
# DLNM LOOP 

for (current_rfm in rfm_list) {
  
  print(paste(">> Berechne Modelle für:", current_rfm))
  
  counter <- 0
  
  # INNER LOOP: PRO AGS
  for (file_path in ags_files) {
    
    counter <- counter + 1
    if(counter %% 50 == 0) cat(paste0("...", counter, "\n")) 
    
    current_ags_name <- basename(file_path) 
    
    # Daten laden
    df <- fread(file_path)
    
    if ("Datum" %in% names(df)) setnames(df, "Datum", "date")
    
    # Check
    if(!current_rfm %in% names(df)) next
    
    # DATUMS-KORREKTUR VOR DER AGGREGATION
    df[, date := as.Date(date)]
    
    # Aggregation auf "All Admissions" pro Tag
    # Summieren aller Diagnosen zur Gesamtzahl
    df[, All_Cases := rowSums(.SD, na.rm = TRUE), .SDcols = vars_to_aggregate]
    
    # Tages-Summen erstellen 
    agsdat <- df[, .(Outcome_Count = sum(All_Cases, na.rm = TRUE)), 
                 by = .(date, T_mean_popw, holiday, get(current_rfm), heat_alerts_dwd)]
    
    setnames(agsdat, "get", current_rfm)
    
    # Variablen für Modell aufbauen
    agsdat[, year  := as.numeric(format(date, "%Y"))]
    agsdat[, month := as.numeric(format(date, "%m"))]
    agsdat[, dow   := as.factor(format(date, "%u"))]
    agsdat[, yday.sm := as.numeric(format(date, "%j"))]
    
    # HW_Status (Hitzewarnsystem-Implementierung) 
    # Ab 2005 = 1 (Aktiv), davor = 0 (Inaktiv)
    agsdat[, HW_Status := ifelse(year >= 2005, 1, 0)]
    
    # Zeitraum filtern
    agsdat <- agsdat[year %in% years_analysis & month %in% months_analysis]
    setorder(agsdat, date)
    
    # SICHERHEITSCHECKS VOR DER MODELLIERUNG
    # Mindestfallzahl-Check (Gegen Modell-Abstürze)
    total_cases <- sum(agsdat$Outcome_Count, na.rm=TRUE)
    n_days <- nrow(agsdat)
    
    if(n_days < 100 || total_cases < 10) {
      dropped_ags_list[[length(dropped_ags_list) + 1]] <- data.frame(
        AGS_File = current_ags_name,
        RFM = current_rfm,
        Diagnosis = "AllAdmissions", 
        Total_Cases = total_cases,
        Total_Days = n_days
      )
      next 
    }
    
    # NAs in der Temperatur Check:
    if(sum(is.na(agsdat$T_mean_popw)) > 0) {
      error_log <<- rbind(error_log, data.frame(AGS = current_ags_name, RFM = current_rfm, Error_Message = "Abbruch: NAs in Temperatur gefunden"))
      next
    }
    
    # Kontrast-Check
    if(length(unique(agsdat$year)) < 2) {
      error_log <<- rbind(error_log, data.frame(AGS = current_ags_name, RFM = current_rfm, Error_Message = "Abbruch: Daten umfassen nur 1 Jahr (Kontrast-Fehler)"))
      next
    }
    
    # NEU Default: 2000-2004 RFM Warnungen 2005-2009 DWD Warnungen 
    
    agsdat[, combined_alerts := ifelse(year < 2005, get(current_rfm), heat_alerts_dwd)]
    
    # DiD/Kollinearitäts-Check anpassen auf die neue kombinierte Variable
    warnings_pre_2005 <- sum(agsdat$combined_alerts[agsdat$year < 2005], na.rm = TRUE)
    warnings_post_2005 <- sum(agsdat$combined_alerts[agsdat$year >= 2005], na.rm = TRUE)
    
    if(warnings_pre_2005 == 0 || warnings_post_2005 == 0) {
      error_log <<- rbind(error_log, data.frame(AGS = current_ags_name, RFM = current_rfm, Error_Message = paste("Abbruch: Warnungen fehlen (Vor 2005:", warnings_pre_2005, "| Nach 2005:", warnings_post_2005, ")")))
      next
    }
    
    # Extrahieren für DLNM
    outcome_data  <- agsdat$Outcome_Count
    temp_data     <- agsdat$T_mean_popw
    eligible_data <- agsdat$combined_alerts
    holiday_data  <- agsdat$holiday
    impl_data     <- agsdat$HW_Status
    
    ## MODELLIERUNG ##
    
    tryCatch({
      
      # Basen erstellen 
      cbtmean <- crossbasis(temp_data, lag=lag, 
                            argvar=list(fun=varfun, degree=vardegree, knots=quantile(temp_data, varper/100, na.rm=TRUE)), 
                            arglag=list(fun="integer"), group=agsdat$year)
      
      datenum <- as.numeric(agsdat$date)
      intfirst <- ((datenum-day1)/(day2-day1))*cbtmean
      intlast  <- ((datenum-day2)/(day1-day2))*cbtmean # intlast für model3
      
      cbalert <- crossbasis(eligible_data, lag=lag,
                            argvar=list(fun="lin"),
                            arglag=list(fun="integer"), group=agsdat$year)
      
      alertint <- impl_data * cbalert
      
      ## QUASI-POISSON MODELS ##
      
      # WITHOUT DLNM FOR TEMPERATURE
      m0 <- glm(outcome_data ~ cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
      
      # WITHOUT TIME-VARYING DLNM FOR TEMPERATURE
      m1 <- glm(outcome_data ~ cbtmean + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
      
      # WITH TIME-VARYING DLNM FOR TEMPERATURE (ACCORDING TO GASPARRINI ET AL 2015 EHP)
      m2 <- glm(outcome_data ~ cbtmean + intfirst + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
      
      
      m3 <- glm(outcome_data ~ cbtmean + intlast + cbalert + alertint + ns(yday.sm, df=dfseas):factor(year) + ns(date, df=round(length(unique(year))/dftrend/10)) + as.factor(dow) + as.factor(holiday_data), family=quasipoisson, data=agsdat, na.action="na.exclude")
      
      # EFFEKTE EXTRAHIEREN 
      mean_temp <- mean(temp_data, na.rm=TRUE) 
      
      
      # A) Temperatur DLNM (für Plots)
      # Model 0 hat keine Temperatur -> extrahieren 1, 2 und 3
      red_temp_1 <- crossreduce(cbtmean, m1, cen=mean_temp)
      red_temp_2 <- crossreduce(cbtmean, m2, cen=mean_temp)
      red_temp_3 <- crossreduce(cbtmean, m3, cen=mean_temp)
      
      # Vektoren in Textstrings verwandeln, mit Pipe
      coef_temp_1 <- paste(coef(red_temp_1), collapse="|")
      coef_temp_2 <- paste(coef(red_temp_2), collapse="|")
      coef_temp_3 <- paste(coef(red_temp_3), collapse="|")
      
      vcov_temp_1 <- paste(as.vector(vcov(red_temp_1)), collapse="|")
      vcov_temp_2 <- paste(as.vector(vcov(red_temp_2)), collapse="|")
      vcov_temp_3 <- paste(as.vector(vcov(red_temp_3)), collapse="|")
      
      # B) Alert-Effekt
      p_alert_0 <- crossreduce(cbalert, m0, at=1)
      p_alert_1 <- crossreduce(cbalert, m1, at=1)
      p_alert_2 <- crossreduce(cbalert, m2, at=1)
      p_alert_3 <- crossreduce(cbalert, m3, at=1)
      
      # C) DiD-Interaktion 
      p_int_0 <- crossreduce(alertint, m0, at=1)
      p_int_1 <- crossreduce(alertint, m1, at=1)
      p_int_2 <- crossreduce(alertint, m2, at=1)
      p_int_3 <- crossreduce(alertint, m3, at=1)
      
      # ZAHLEN IN DIE CSV SCHREIBEN
      res_tmp <- data.frame(
        AGS_File = current_ags_name,
        RFM = current_rfm,
        Diagnosis = "AllAdmissions",
        Model = c("model0", "model1", "model2", "model3"),
        Total_Cases = total_cases,                                     
        QAIC = c(QAIC(m0), QAIC(m1), QAIC(m2), QAIC(m3)), 
        
        # lokale R-Plots bauen
        Temp_Mean_Cen = mean_temp, 
        Temp_Coef_Vector = c(NA, coef_temp_1, coef_temp_2, coef_temp_3), # Model 0 ist NA
        Temp_Vcov_Matrix = c(NA, vcov_temp_1, vcov_temp_2, vcov_temp_3),
        
        # Alert-Effekt
        RR_Alert = c(p_alert_0$RRfit, p_alert_1$RRfit, p_alert_2$RRfit, p_alert_3$RRfit),                
        CI_Low_Alert = c(p_alert_0$RRlow, p_alert_1$RRlow, p_alert_2$RRlow, p_alert_3$RRlow),   
        CI_Up_Alert = c(p_alert_0$RRhigh, p_alert_1$RRhigh, p_alert_2$RRhigh, p_alert_3$RRhigh), 
        
        # DiD Interaktion
        DiD_Coef = c(coef(p_int_0), coef(p_int_1), coef(p_int_2), coef(p_int_3)),                          
        DiD_Variance = c(vcov(p_int_0), vcov(p_int_1), vcov(p_int_2), vcov(p_int_3)),                      
        RR_DiD = c(p_int_0$RRfit, p_int_1$RRfit, p_int_2$RRfit, p_int_3$RRfit),                
        CI_Low_DiD = c(p_int_0$RRlow, p_int_1$RRlow, p_int_2$RRlow, p_int_3$RRlow),   
        CI_Up_DiD = c(p_int_0$RRhigh, p_int_1$RRhigh, p_int_2$RRhigh, p_int_3$RRhigh)
      )
      
      results_list[[length(results_list) + 1]] <- res_tmp
      
    }, error = function(e) {
      
      # Error File
      error_log <<- rbind(error_log, data.frame(AGS = current_ags_name, RFM = current_rfm, Error_Message = e$message))
      
    }) # Ende tryCatch   
  } # Ende AGS Loop
} # Ende RFM Loop

# -------------------------------------------------------------------------
# EXPORT DER ERGEBNISSE

if (length(results_list) > 0) {
  final_df <- rbindlist(results_list)
  out_name <- paste0("Results_FirstStage_AllAdmissions.csv")
  write.csv(final_df, file.path(output_path, out_name), row.names = FALSE)
  print(paste("Ergebnisse gespeichert unter:", out_name))
}

# EXPORT DER AUSGESCHLOSSENEN LANDKREISE
if (length(dropped_ags_list) > 0) {
  dropped_df <- rbindlist(dropped_ags_list)
  dropped_df <- unique(dropped_df, by = c("AGS_File", "Diagnosis")) 
  
  dropped_out_name <- "Dropped_AGS_AllAdmissions.csv"
  write.csv(dropped_df, file.path(output_path, dropped_out_name), row.names = FALSE)
  print(paste("Liste der ausgeschlossenen Kreise liegt hier:", dropped_out_name))
} else {
  print("Kein Landkreis musste wegen zu geringer Fallzahlen ausgeschlossen werden.")
}

# Export Fehlerprotokoll
if (nrow(error_log) > 0) {
  error_out_name <- "Error_Log_FirstStage_all_admissions.csv"
  write.csv(error_log, file.path(output_path, error_out_name), row.names = FALSE)
  print(paste("Abstürze vorhanden, Fehlerprotokoll liegt hier:", error_out_name))
} else {
  print("Keine Modell-Abstürze")
}
