# ----------------------------------------
# YEP Time-to-Anoxia Analysis with Mass/Volume Normalization
# ----------------------------------------
rm(list=ls(all=T))
library(deSolve)   
library(FME)       
library(ggplot2)   
library(dplyr)
library(tidyr)
library(gridExtra)

# -------------------------------
# Create output directory
# -------------------------------
out_dir <- "modeling_outputs"
if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -------------------------------
# Load and process data
# -------------------------------
file_path <- "YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv"
mapping_file_path <- "C:/Users/gara009/OneDrive - PNNL/RC-SFA - Documents/Study_YEP/DO/05_PublishReadyData/Merged_Firesting_Mapping_2025-10-24.csv"
mass_volume_file_path <- "YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Water_Mass_Volume.csv"

# Load mass/volume data
mass_volume_df <- read.csv(mass_volume_file_path, stringsAsFactors = FALSE, skip = 2)%>%
  filter(grepl('YEP', Sample_Name))

# Load DO and mapping data
mapping_df <- read.csv(mapping_file_path, stringsAsFactors = FALSE)
df <- read.csv(file_path, skip = 14, header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("Field_Name","Sample_Name","IGSN","Material","DateTime",
                  "Elapsed_Seconds","Temperature_degreesC","DO_mg_per_L",
                  "Firesting_Serial_Number","Methods_Deviation")

df <- df %>%
  inner_join(mapping_df %>% select(Sample_Name, Start_Time_Elapsed_Seconds = Elapsed_Seconds_Start, 
                                   End_Time_Elapsed_Seconds = Elapsed_Seconds_End),
             by = "Sample_Name") %>%
  filter(Elapsed_Seconds >= Start_Time_Elapsed_Seconds,
         Elapsed_Seconds <= End_Time_Elapsed_Seconds) %>%
  select(-Start_Time_Elapsed_Seconds, -End_Time_Elapsed_Seconds)

df <- df %>%
  mutate(Elapsed_Seconds = as.numeric(Elapsed_Seconds),
         DO_mg_per_L = as.numeric(DO_mg_per_L)) %>%
  filter(!is.na(Elapsed_Seconds), !is.na(DO_mg_per_L)) %>%
  filter(DO_mg_per_L > 0)

# -------------------------------
# Constants (from paper context)
# -------------------------------
Km <- 10      # Half-saturation constant (mg/L)
CB <- 1.0     # Initial microbial biomass (dimensionless)

# -------------------------------
# YEP Treatment classification
# -------------------------------
classify_yep_treatment <- function(sample_names) {
  treatment_types <- case_when(
    grepl("YEP1.*S", sample_names, ignore.case = TRUE) ~ "Dry_Control",
    grepl("YEP1.*U", sample_names, ignore.case = TRUE) ~ "Dry_Unburned_DOC",
    grepl("YEP1.*H", sample_names, ignore.case = TRUE) ~ "Dry_HighBurn_DOC",
    grepl("YEP2.*S", sample_names, ignore.case = TRUE) ~ "Wet_Control",
    grepl("YEP2.*U", sample_names, ignore.case = TRUE) ~ "Wet_Unburned_DOC",
    grepl("YEP2.*H", sample_names, ignore.case = TRUE) ~ "Wet_HighBurn_DOC",
    TRUE ~ "Unknown"
  )
  return(treatment_types)
}

# -------------------------------
# Calculate sample-specific conversion factors
# -------------------------------
calculate_conversion_factors <- function(mass_volume_df) {
  # Filter for incubation samples only
  incubation_samples <- mass_volume_df %>%
    filter(grepl("_INC-", Sample_Name)) %>%
    mutate(
      Dry_Sediment_Mass_g = as.numeric(Dry_Sediment_Mass_g),
      Water_Mass_g = as.numeric(Water_Mass_g),
      # Calculate conversion factor: mg/L to mg/g dry sediment
      # Conversion = Water_Volume_L / (Dry_Sediment_Mass_g)
      Water_Volume_L = Water_Mass_g / 1000,  # Convert g to L (assuming water density = 1)
      Conversion_Factor = Water_Volume_L / Dry_Sediment_Mass_g
    ) %>%
    select(Sample_Name, Dry_Sediment_Mass_g, Water_Volume_L, Conversion_Factor)
  
  return(incubation_samples)
}

conversion_factors <- calculate_conversion_factors(mass_volume_df)

# -------------------------------
# Data trimming function - Remove first and last 2 minutes
# -------------------------------
trim_time_series <- function(data_subset) {
  trim_hours <- 2/60  # 2 minutes = 0.033 hours
  
  data_trimmed <- data_subset %>%
    filter(time_hr >= (min(time_hr) + trim_hours) & 
             time_hr <= (max(time_hr) - trim_hours))
  
  return(data_trimmed)
}

# -------------------------------
# Model solver functions - Correct units from paper
# -------------------------------

# 1. Linear abiotic: dDO/dt = kL × [DO] (kL in h⁻¹)
solveLinear <- function(pars, times, DO0) {
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      dDO <- -kL * DO
      list(c(dDO))
    })
  }
  state <- c(DO = DO0)
  tryCatch({
    out <- ode(y = state, times = times, func = derivs, parms = pars, method = "lsoda")
    result <- as.data.frame(out[, c("time","DO")])
    if(any(is.na(result$DO)) || any(result$DO < 0)) {
      stop("Invalid solution")
    }
    return(result)
  }, error = function(e) {
    result <- data.frame(
      time = times,
      DO = pmax(0, DO0 * exp(-pars$kL * times))
    )
    return(result)
  })
}

# 2. Biotic: dDO/dt = Vmax × [DO]/([DO] + Km) × CB (Vmax in h⁻¹)
solveBiotic <- function(pars, times, DO0, CB, Km) {
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      r1 <- Vmax * DO/(DO + Km + 1e-10)
      dDO <- -r1 * CB1
      dCB1 <- r1 * CB1
      list(c(dDO, dCB1))
    })
  }
  state <- c(DO = DO0, CB1 = CB)
  tryCatch({
    out <- ode(y = state, times = times, func = derivs, parms = pars, method = "lsoda")
    result <- as.data.frame(out[, c("time","DO")])
    if(any(is.na(result$DO)) || any(result$DO < 0)) {
      stop("Invalid solution")
    }
    return(result)
  }, error = function(e) {
    result <- data.frame(
      time = times,
      DO = pmax(0, DO0 * exp(-pars$Vmax * (DO0/(DO0 + Km)) * CB * times))
    )
    return(result)
  })
}

# 3. Combined: dDO/dt = (Vmax × [DO]/([DO] + Km) × CB) + kL × [DO]
solveCombined <- function(pars, times, DO0, CB, Km) {
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      r1 <- Vmax * DO/(DO + Km + 1e-10)
      dDO <- -r1 * CB1 - kL * DO
      dCB1 <- r1 * CB1
      list(c(dDO, dCB1))
    })
  }
  state <- c(DO = DO0, CB1 = CB)
  tryCatch({
    out <- ode(y = state, times = times, func = derivs, parms = pars, method = "lsoda")
    result <- as.data.frame(out[, c("time","DO")])
    if(any(is.na(result$DO)) || any(result$DO < 0)) {
      stop("Invalid solution")
    }
    return(result)
  }, error = function(e) {
    total_rate <- pars$kL + pars$Vmax * (DO0/(DO0 + Km)) * CB
    result <- data.frame(
      time = times,
      DO = pmax(0, DO0 * exp(-total_rate * times))
    )
    return(result)
  })
}

# -------------------------------
# Helper functions
# -------------------------------
calculate_time_to_anoxia <- function(time_vec, DO_vec, threshold = 0.1) {
  anoxic_idx <- which(DO_vec <= threshold)[1]
  if(is.na(anoxic_idx)) {
    return(max(time_vec))
  }
  return(time_vec[anoxic_idx])
}

compute_info <- function(SSR, k, n) {
  AIC <- n * log(SSR / n) + 2 * k
  BIC <- n * log(SSR / n) + k * log(n)
  list(AIC = AIC, BIC = BIC)
}

# -------------------------------
# Complete fitting function with mass/volume corrections
# -------------------------------
fit_all_yep_models <- function(dat, sample_name, treatment_type = "Unknown", conversion_info = NULL) {
  
  # Trim first and last 2 minutes, then clean data
  dat_trimmed <- trim_time_series(dat)
  
  Data <- dat_trimmed %>%
    arrange(time_hr) %>%
    select(time = time_hr, DO = DO_mg_per_L) %>%
    group_by(time) %>%
    summarise(DO = mean(DO, na.rm = TRUE), .groups = 'drop') %>%
    filter(DO > 0, !is.na(DO), !is.na(time)) %>%
    arrange(time) %>%
    as.data.frame()
  
  n <- nrow(Data)
  if (n < 5) {
    message("Skipping ", sample_name, ": insufficient data points after trimming (n=", n, ")")
    return(NULL)
  }
  
  # Rebase time to start from 0
  Data$time <- Data$time - min(Data$time)
  
  DO0 <- Data$DO[1]
  time_to_anoxia_obs <- calculate_time_to_anoxia(Data$time, Data$DO)
  
  # Get conversion factor for this sample
  if(!is.null(conversion_info)) {
    dry_mass_g <- conversion_info$Dry_Sediment_Mass_g
    water_vol_L <- conversion_info$Water_Volume_L
    conversion_factor <- conversion_info$Conversion_Factor
  } else {
    # Use default values if no specific mass/volume data
    dry_mass_g <- 10
    water_vol_L <- 0.030
    conversion_factor <- 0.033
  }
  
  message("Fitting YEP models for ", sample_name, " (n=", n, ", treatment=", treatment_type, ")")
  message("  Initial DO: ", round(DO0, 2), " mg/L, Time to anoxia: ", round(time_to_anoxia_obs, 1), " h")
  message("  Dry mass: ", round(dry_mass_g, 2), " g, Water vol: ", round(water_vol_L*1000, 1), " mL")
  message("  Conversion factor: ", round(conversion_factor, 4), " (vs. standard 0.050)")
  
  # Estimate starting parameters
  linear_rate_est <- abs(coef(lm(DO ~ time, data = Data))[["time"]])
  
  # Define objective functions
  
  # 1. Linear abiotic objective
  Objective_LIN <- function(x) {
    if(x[["kL"]] <= 0 || x[["kL"]] > 50) return(1e10)
    pars <- list(kL = x[["kL"]])
    out <- solveLinear(pars, times = Data$time, DO0 = DO0)
    if(any(is.na(out$DO))) return(1e10)
    modCost(model = out, obs = Data)
  }
  
  # 2. Biotic objective
  Objective_BIO <- function(x) {
    if(x[["Vmax"]] <= 0 || x[["Vmax"]] > 100) return(1e10)
    pars <- list(Vmax = x[["Vmax"]])
    out <- solveBiotic(pars, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    if(any(is.na(out$DO))) return(1e10)
    modCost(model = out, obs = Data)
  }
  
  # 3. Combined objective
  Objective_COM <- function(x) {
    if(x[["Vmax"]] <= 0 || x[["Vmax"]] > 100 || x[["kL"]] <= 0 || x[["kL"]] > 50) return(1e10)
    pars <- list(Vmax = x[["Vmax"]], kL = x[["kL"]])
    out <- solveCombined(pars, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    if(any(is.na(out$DO))) return(1e10)
    modCost(model = out, obs = Data)
  }
  
  # Fit all models
  fit_LIN <- fit_BIO <- fit_COM <- NULL
  
  # Linear fit
  try({
    fit_LIN <- modFit(p = c(kL = linear_rate_est/DO0), f = Objective_LIN,
                      lower = c(kL = 0.001), upper = c(kL = 10),
                      control = list(maxiter = 200))
  }, silent = TRUE)
  
  # Biotic fit
  try({
    fit_BIO <- modFit(p = c(Vmax = linear_rate_est), f = Objective_BIO,
                      lower = c(Vmax = 0.01), upper = c(Vmax = 50),
                      control = list(maxiter = 200))
  }, silent = TRUE)
  
  # Combined fit
  try({
    fit_COM <- modFit(p = c(Vmax = linear_rate_est/2, kL = linear_rate_est/DO0/2), 
                      f = Objective_COM,
                      lower = c(Vmax = 0.01, kL = 0.001), 
                      upper = c(Vmax = 50, kL = 10),
                      control = list(maxiter = 200))
  }, silent = TRUE)
  
  # Generate predictions
  model_LIN <- model_BIO <- model_COM <- NULL
  
  if(!is.null(fit_LIN)) {
    pars_LIN <- list(kL = fit_LIN$par[["kL"]])
    model_LIN <- solveLinear(pars_LIN, times = Data$time, DO0 = DO0)
    message("  Linear fit: kL = ", round(fit_LIN$par[["kL"]], 4), " h⁻¹, SSR = ", round(fit_LIN$ssr, 2))
  }
  
  if(!is.null(fit_BIO)) {
    pars_BIO <- list(Vmax = fit_BIO$par[["Vmax"]])
    model_BIO <- solveBiotic(pars_BIO, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    message("  Biotic fit: Vmax = ", round(fit_BIO$par[["Vmax"]], 4), " h⁻¹, SSR = ", round(fit_BIO$ssr, 2))
  }
  
  if(!is.null(fit_COM)) {
    pars_COM <- list(Vmax = fit_COM$par[["Vmax"]], kL = fit_COM$par[["kL"]])
    model_COM <- solveCombined(pars_COM, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    message("  Combined fit: Vmax = ", round(fit_COM$par[["Vmax"]], 4), " h⁻¹",
            ", kL = ", round(fit_COM$par[["kL"]], 4), " h⁻¹, SSR = ", round(fit_COM$ssr, 2))
  }
  
  # Create plots with sample name only
  tryCatch({
    pdf_file <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", sample_name), ".pdf"))
    pdf(pdf_file, width = 14, height = 10)
    
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    
    # All models plot
    plot(Data$time, Data$DO, pch = 16, col = "orange", cex = 1.5,
         xlab = "Time (h)", ylab = "DO (mg/L)",
         main = paste(sample_name, "\nTreatment:", treatment_type),
         ylim = c(0, max(Data$DO) * 1.1))
    
    if(!is.null(model_LIN)) lines(model_LIN$time, model_LIN$DO, col = "blue", lwd = 3, lty = 1)
    if(!is.null(model_BIO)) lines(model_BIO$time, model_BIO$DO, col = "red", lwd = 3, lty = 2)
    if(!is.null(model_COM)) lines(model_COM$time, model_COM$DO, col = "green", lwd = 3, lty = 3)
    
    # SSR text box
    ssr_text <- "SSR Values:\n"
    if(!is.null(fit_LIN)) ssr_text <- paste0(ssr_text, sprintf("Linear: %.1f\n", fit_LIN$ssr))
    if(!is.null(fit_BIO)) ssr_text <- paste0(ssr_text, sprintf("Biotic: %.1f\n", fit_BIO$ssr))
    if(!is.null(fit_COM)) ssr_text <- paste0(ssr_text, sprintf("Combined: %.1f", fit_COM$ssr))
    
    # Mass normalization info
    mass_text <- paste0("Normalization:\n", 
                        sprintf("Dry mass: %.1f g\n", dry_mass_g),
                        sprintf("Water: %.1f mL\n", water_vol_L*1000),
                        sprintf("Conv. factor: %.3f", conversion_factor))
    
    text(x = max(Data$time) * 0.6, y = max(Data$DO) * 0.8,
         labels = ssr_text, adj = c(0, 1), bg = "white", cex = 0.9, family = "mono")
    
    text(x = max(Data$time) * 0.05, y = max(Data$DO) * 0.8,
         labels = mass_text, adj = c(0, 1), bg = "white", cex = 0.8, family = "mono")
    
    legend("topright", c("Observed", "Linear Abiotic", "Nonlinear Biotic", "Combined"),
           pch = c(16, NA, NA, NA), lty = c(NA, 1, 2, 3),
           col = c("orange", "blue", "red", "green"), lwd = c(NA, 3, 3, 3), bty = "n")
    
    # Individual model plots
    models_list <- list(
      list(fit = fit_LIN, model = model_LIN, title = "Linear Abiotic", color = "blue"),
      list(fit = fit_BIO, model = model_BIO, title = "Nonlinear Biotic", color = "red"),
      list(fit = fit_COM, model = model_COM, title = "Combined", color = "green")
    )
    
    for(i in 1:3) {
      if(!is.null(models_list[[i]]$fit)) {
        plot(Data$time, Data$DO, pch = 16, col = "orange", cex = 1.2,
             xlab = "Time (h)", ylab = "DO (mg/L)",
             main = paste(models_list[[i]]$title, "\nSSR =", round(models_list[[i]]$fit$ssr, 1)))
        lines(models_list[[i]]$model$time, models_list[[i]]$model$DO, 
              col = models_list[[i]]$color, lwd = 3)
        
        # Add parameter values
        if(i == 1 && !is.null(fit_LIN)) {
          text(x = max(Data$time) * 0.05, y = max(Data$DO) * 0.8,
               labels = paste("kL =", round(fit_LIN$par[["kL"]], 3), "h⁻¹"),
               adj = c(0, 1), bg = "white", cex = 0.9)
        } else if(i == 2 && !is.null(fit_BIO)) {
          text(x = max(Data$time) * 0.05, y = max(Data$DO) * 0.8,
               labels = paste("Vmax =", round(fit_BIO$par[["Vmax"]], 3), "h⁻¹"),
               adj = c(0, 1), bg = "white", cex = 0.9)
        } else if(i == 3 && !is.null(fit_COM)) {
          param_text <- paste("Vmax =", round(fit_COM$par[["Vmax"]], 3), "h⁻¹\n",
                              "kL =", round(fit_COM$par[["kL"]], 3), "h⁻¹")
          text(x = max(Data$time) * 0.05, y = max(Data$DO) * 0.8,
               labels = param_text, adj = c(0, 1), bg = "white", cex = 0.9)
        }
        
        legend("topright", c("Observed", "Model"), 
               pch = c(16, NA), lty = c(NA, 1),
               col = c("orange", models_list[[i]]$color), lwd = c(NA, 3), bty = "n")
      } else {
        plot(1, 1, type = "n", xlab = "", ylab = "", main = paste(models_list[[i]]$title, "- Fit Failed"))
        text(1, 1, "Model fit failed", cex = 1.5, col = "red")
      }
    }
    
    dev.off()
    message("  Wrote: ", pdf_file)
    
  }, error = function(e) {
    message("  Error creating plots: ", e$message)
  })
  
  # Calculate normalized rates (per gram dry sediment)
  DO_rate_linear_per_g <- if(!is.null(fit_LIN)) fit_LIN$par[["kL"]] * DO0 * conversion_factor else NA
  DO_rate_biotic_per_g <- if(!is.null(fit_BIO)) fit_BIO$par[["Vmax"]] * (DO0/(DO0 + Km)) * CB * conversion_factor else NA
  DO_rate_combined_per_g <- if(!is.null(fit_COM)) {
    bio_rate <- fit_COM$par[["Vmax"]] * (DO0/(DO0 + Km)) * CB
    abi_rate <- fit_COM$par[["kL"]] * DO0
    (bio_rate + abi_rate) * conversion_factor
  } else NA
  
  # Create comprehensive summary
  summary_row <- data.frame(
    Sample_Name = sample_name,
    Treatment_Type = treatment_type,
    n_points = n,
    Initial_DO_mg_per_L = round(DO0, 2),
    Time_to_Anoxia_h = round(time_to_anoxia_obs, 2),
    Dry_Sediment_Mass_g = round(dry_mass_g, 2),
    Water_Volume_mL = round(water_vol_L * 1000, 1),
    Conversion_Factor = round(conversion_factor, 4),
    # Linear model (kL in h⁻¹)
    Linear_kL_per_h = ifelse(!is.null(fit_LIN), round(fit_LIN$par[["kL"]], 4), NA),
    Linear_Rate_mg_per_L_per_h = ifelse(!is.null(fit_LIN), round(fit_LIN$par[["kL"]] * DO0, 4), NA),
    Linear_Rate_mg_per_g_per_h = round(DO_rate_linear_per_g, 4),
    Linear_SSR = ifelse(!is.null(fit_LIN), round(fit_LIN$ssr, 2), NA),
    Linear_AIC = ifelse(!is.null(fit_LIN), round(compute_info(fit_LIN$ssr, 1, n)$AIC, 2), NA),
    # Biotic model (Vmax in h⁻¹)
    Biotic_Vmax_per_h = ifelse(!is.null(fit_BIO), round(fit_BIO$par[["Vmax"]], 4), NA),
    Biotic_Rate_mg_per_L_per_h = round(DO_rate_biotic_per_g / conversion_factor, 4),
    Biotic_Rate_mg_per_g_per_h = round(DO_rate_biotic_per_g, 4),
    Biotic_SSR = ifelse(!is.null(fit_BIO), round(fit_BIO$ssr, 2), NA),
    Biotic_AIC = ifelse(!is.null(fit_BIO), round(compute_info(fit_BIO$ssr, 1, n)$AIC, 2), NA),
    # Combined model
    Combined_Vmax_per_h = ifelse(!is.null(fit_COM), round(fit_COM$par[["Vmax"]], 4), NA),
    Combined_kL_per_h = ifelse(!is.null(fit_COM), round(fit_COM$par[["kL"]], 4), NA),
    Combined_Rate_mg_per_L_per_h = round(DO_rate_combined_per_g / conversion_factor, 4),
    Combined_Rate_mg_per_g_per_h = round(DO_rate_combined_per_g, 4),
    Combined_SSR = ifelse(!is.null(fit_COM), round(fit_COM$ssr, 2), NA),
    Combined_AIC = ifelse(!is.null(fit_COM), round(compute_info(fit_COM$ssr, 2, n)$AIC, 2), NA),
    stringsAsFactors = FALSE
  )
  
  # Best model by AIC
  aic_values <- c(summary_row$Linear_AIC, summary_row$Biotic_AIC, summary_row$Combined_AIC)
  aic_names <- c("Linear", "Biotic", "Combined")
  valid_aic <- !is.na(aic_values)
  
  if(any(valid_aic)) {
    best_idx <- which.min(aic_values[valid_aic])
    summary_row$Best_Model <- aic_names[valid_aic][best_idx]
    summary_row$Best_AIC <- min(aic_values[valid_aic])
  } else {
    summary_row$Best_Model <- "None"
    summary_row$Best_AIC <- NA
  }
  
  return(list(summary = summary_row, 
              fits = list(linear = fit_LIN, biotic = fit_BIO, combined = fit_COM),
              data = Data))
}

# -------------------------------
# Create comprehensive plots
# -------------------------------
create_yep_summary_plots <- function(summary_df) {
  
  # Filter valid data
  plot_data <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    mutate(
      Moisture = ifelse(grepl("Dry", Treatment_Type), "Dry", "Wet"),
      DOC_Treatment = case_when(
        grepl("Control", Treatment_Type) ~ "Control",
        grepl("Unburned", Treatment_Type) ~ "Unburned_DOC",
        grepl("HighBurn", Treatment_Type) ~ "HighBurn_DOC",
        TRUE ~ "Other"
      )
    )
  
  # 1. Mass-normalized rate comparison
  pdf(file.path(out_dir, "Rate_Comparison_Mass_Normalized.pdf"), width = 12, height = 8)
  
  rate_summary <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Linear_Rate_Mean = mean(Linear_Rate_mg_per_g_per_h, na.rm = TRUE),
      Linear_Rate_SE = sd(Linear_Rate_mg_per_g_per_h, na.rm = TRUE) / sqrt(n()),
      Biotic_Rate_Mean = mean(Biotic_Rate_mg_per_g_per_h, na.rm = TRUE),
      Biotic_Rate_SE = sd(Biotic_Rate_mg_per_g_per_h, na.rm = TRUE) / sqrt(n()),
      Combined_Rate_Mean = mean(Combined_Rate_mg_per_g_per_h, na.rm = TRUE),
      Combined_Rate_SE = sd(Combined_Rate_mg_per_g_per_h, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  par(mfrow = c(1, 2))
  
  # Bar plot 1: Linear vs Combined abiotic rates
  rates_1 <- as.matrix(rate_summary[, c("Linear_Rate_Mean", "Combined_Rate_Mean")])
  se_1 <- as.matrix(rate_summary[, c("Linear_Rate_SE", "Combined_Rate_SE")])
  rownames(rates_1) <- rate_summary$Treatment_Type
  
  bp1 <- barplot(t(rates_1), beside = TRUE, 
                 col = c("lightblue", "blue"),
                 main = "Abiotic Rates (Mass Normalized)",
                 ylab = "Rate (mg DO/g dry sediment/h)",
                 legend.text = c("Linear", "Combined"),
                 args.legend = list(x = "topright", bty = "n"),
                 las = 2)
  
  # Bar plot 2: Biotic vs Combined biotic rates  
  rates_2 <- as.matrix(rate_summary[, c("Biotic_Rate_Mean", "Combined_Rate_Mean")])
  se_2 <- as.matrix(rate_summary[, c("Biotic_Rate_SE", "Combined_Rate_SE")])
  
  bp2 <- barplot(t(rates_2), beside = TRUE,
                 col = c("lightcoral", "red"),
                 main = "Biotic Rates (Mass Normalized)",
                 ylab = "Rate (mg DO/g dry sediment/h)",
                 legend.text = c("Biotic", "Combined"),
                 args.legend = list(x = "topright", bty = "n"),
                 las = 2)
  
  dev.off()
  
  # 2. Time to anoxia comparison with mass info
  pdf(file.path(out_dir, "Time_to_Anoxia_with_Mass_Info.pdf"), width = 12, height = 8)
  
  anoxia_data <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Mean_Time = mean(Time_to_Anoxia_h, na.rm = TRUE),
      SE_Time = sd(Time_to_Anoxia_h, na.rm = TRUE) / sqrt(n()),
      Mean_Dry_Mass = mean(Dry_Sediment_Mass_g, na.rm = TRUE),
      Mean_Water_Vol = mean(Water_Volume_mL, na.rm = TRUE),
      .groups = 'drop'
    )
  
  par(mfrow = c(1, 2))
  
  # Time to anoxia
  bp3 <- barplot(anoxia_data$Mean_Time, 
                 names.arg = gsub("_", "\n", anoxia_data$Treatment_Type),
                 col = rainbow(nrow(anoxia_data)),
                 main = "Time to Anoxia by Treatment",
                 ylab = "Time to Anoxia (h)",
                 ylim = c(0, max(anoxia_data$Mean_Time + anoxia_data$SE_Time) * 1.1))
  
  arrows(bp3, anoxia_data$Mean_Time - anoxia_data$SE_Time,
         bp3, anoxia_data$Mean_Time + anoxia_data$SE_Time,
         length = 0.05, angle = 90, code = 3)
  
  # Sample mass comparison
  mass_matrix <- as.matrix(anoxia_data[, c("Mean_Dry_Mass", "Mean_Water_Vol")])
  rownames(mass_matrix) <- anoxia_data$Treatment_Type
  
  barplot(t(mass_matrix), beside = TRUE,
          col = c("brown", "lightblue"),
          main = "Sample Mass/Volume by Treatment", 
          ylab = "Mass/Volume",
          legend.text = c("Dry Sediment (g)", "Water (mL)"),
          args.legend = list(x = "topright", bty = "n"),
          las = 2)
  
  dev.off()
  
  # 3. Conversion factor comparison
  pdf(file.path(out_dir, "Conversion_Factor_Comparison.pdf"), width = 10, height = 6)
  
  conv_data <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      Mean_Conversion = mean(Conversion_Factor, na.rm = TRUE),
      SE_Conversion = sd(Conversion_Factor, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  bp4 <- barplot(conv_data$Mean_Conversion,
                 names.arg = gsub("_", "\n", conv_data$Treatment_Type),
                 col = "lightgreen",
                 main = "Sample-Specific Conversion Factors",
                 ylab = "Conversion Factor (vs. standard 0.050)",
                 ylim = c(0, max(conv_data$Mean_Conversion + conv_data$SE_Conversion) * 1.1))
  
  arrows(bp4, conv_data$Mean_Conversion - conv_data$SE_Conversion,
         bp4, conv_data$Mean_Conversion + conv_data$SE_Conversion,
         length = 0.05, angle = 90, code = 3)
  
  # Add reference line for standard conversion
  abline(h = 0.050, col = "red", lwd = 2, lty = 2)
  text(max(bp4) * 0.8, 0.055, "Standard conversion\n(0.050)", col = "red", cex = 0.8)
  
  dev.off()
  
  message("Created mass-normalized summary plots")
}

# -------------------------------
# Main analysis execution
# -------------------------------

# Process data with mass/volume corrections
df2 <- df %>%
  group_by(Sample_Name) %>%
  arrange(Elapsed_Seconds, .by_group = TRUE) %>%
  mutate(time_hr = (Elapsed_Seconds - min(Elapsed_Seconds)) / 3600) %>%
  ungroup() %>%
  mutate(Treatment_Type = classify_yep_treatment(Sample_Name))

# Check treatment classification
treatment_summary <- df2 %>% count(Treatment_Type, sort = TRUE)
print("YEP Treatment Classification:")
print(treatment_summary)

# Check conversion factors available
print("Available Conversion Factors:")
print(conversion_factors)

samples <- unique(df2$Sample_Name)

message("\n=== YEP Analysis with Mass/Volume Corrections ===")
message("Samples: ", length(samples))
message("Output directory: ", out_dir)
message("Trimming: First and last 2 minutes removed")
message("Mass/Volume normalization: Sample-specific")

# Run analysis for all samples
all_results <- lapply(samples, function(sname) {
  dat_s <- df2 %>% filter(Sample_Name == sname)
  treatment_type <- unique(dat_s$Treatment_Type)[1]
  
  # Get conversion info for this sample
  conversion_info <- conversion_factors %>% filter(Sample_Name == sname)
  if(nrow(conversion_info) == 0) {
    conversion_info <- NULL
    message("Warning: No mass/volume data found for ", sname, ", using defaults")
  } else {
    conversion_info <- conversion_info[1, ]  # Take first match
  }
  
  fit_all_yep_models(dat_s, sname, treatment_type, conversion_info)
})

# Combine results
summary_df <- bind_rows(lapply(all_results, function(x) if(!is.null(x)) x$summary else NULL))

# Save results
if(nrow(summary_df) > 0) {
  summary_file <- file.path(out_dir, "YEP_Complete_Analysis_Mass_Normalized.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  message("\nResults saved: ", summary_file)
  
  # Create summary plots
  create_yep_summary_plots(summary_df)
  
  # Print enhanced summary
  message("\n=== YEP MASS-NORMALIZED RESULTS SUMMARY ===")
  
  # Time to anoxia with mass info
  anoxia_summary <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Mean_Time_h = round(mean(Time_to_Anoxia_h, na.rm = TRUE), 1),
      SD_Time_h = round(sd(Time_to_Anoxia_h, na.rm = TRUE), 1),
      Mean_Dry_Mass_g = round(mean(Dry_Sediment_Mass_g, na.rm = TRUE), 1),
      Mean_Conversion = round(mean(Conversion_Factor, na.rm = TRUE), 3),
      .groups = 'drop'
    )
  
  print("Time to Anoxia with Sample Info:")
  print(anoxia_summary)
  
  # Mass-normalized rates
  rate_summary <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      Linear_Rate_mg_g_h = paste0(round(mean(Linear_Rate_mg_per_g_per_h, na.rm = TRUE), 3), " ± ", 
                                  round(sd(Linear_Rate_mg_per_g_per_h, na.rm = TRUE), 3)),
      Biotic_Rate_mg_g_h = paste0(round(mean(Biotic_Rate_mg_per_g_per_h, na.rm = TRUE), 3), " ± ", 
                                  round(sd(Biotic_Rate_mg_per_g_per_h, na.rm = TRUE), 3)),
      Combined_Rate_mg_g_h = paste0(round(mean(Combined_Rate_mg_per_g_per_h, na.rm = TRUE), 3), " ± ",
                                    round(sd(Combined_Rate_mg_per_g_per_h, na.rm = TRUE), 3)),
      .groups = 'drop'
    )
  
  print("\nMass-Normalized Rates (mg DO/g dry sediment/h):")
  print(rate_summary)
  
} else {
  message("No successful analyses completed.")
}

message("\n=== Analysis Complete! ===")
message("All outputs saved in: ", out_dir, "/")
message("Files generated:")
message("- Individual sample PDFs: [SampleName].pdf")
message("- Complete analysis CSV: YEP_Complete_Analysis_Mass_Normalized.csv") 
message("- Mass-normalized rate plots: Rate_Comparison_Mass_Normalized.pdf")
message("- Time to anoxia plots: Time_to_Anoxia_with_Mass_Info.pdf")
message("- Conversion factor plot: Conversion_Factor_Comparison.pdf")


####################################################################
# ----------------------------------------
# Generate Slides 6, 7, and 8 Figures
# ----------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Assume you have summary_df from the main analysis

# -------------------------------
# Slide 6: Rate Comparison Figures (h⁻¹ units)
# -------------------------------
create_slide6_figures <- function(summary_df, out_dir = "modeling_outputs") {
  
  # Filter and prepare data
  plot_data <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      # Abiotic rates (h⁻¹)
      Linear_kL_Mean = mean(Linear_kL_per_h, na.rm = TRUE),
      Linear_kL_SE = sd(Linear_kL_per_h, na.rm = TRUE) / sqrt(n()),
      Combined_kL_Mean = mean(Combined_kL_per_h, na.rm = TRUE),
      Combined_kL_SE = sd(Combined_kL_per_h, na.rm = TRUE) / sqrt(n()),
      # Biotic rates (h⁻¹)
      Biotic_Vmax_Mean = mean(Biotic_Vmax_per_h, na.rm = TRUE),
      Biotic_Vmax_SE = sd(Biotic_Vmax_per_h, na.rm = TRUE) / sqrt(n()),
      Combined_Vmax_Mean = mean(Combined_Vmax_per_h, na.rm = TRUE),
      Combined_Vmax_SE = sd(Combined_Vmax_per_h, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Figure A: Abiotic Rates (kL, h⁻¹)
  p1 <- ggplot(plot_data) +
    geom_col(aes(x = Treatment_Type, y = Linear_kL_Mean), 
             fill = "lightblue", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(x = Treatment_Type, 
                      ymin = Linear_kL_Mean - Linear_kL_SE,
                      ymax = Linear_kL_Mean + Linear_kL_SE),
                  width = 0.2) +
    labs(title = "Abiotic Rates (kL)",
         x = "YEP Treatment",
         y = "Fitted Rate (h⁻¹)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Figure B: Biotic Rates (Vmax, h⁻¹)
  p2 <- ggplot(plot_data) +
    geom_col(aes(x = Treatment_Type, y = Biotic_Vmax_Mean), 
             fill = "lightcoral", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(x = Treatment_Type, 
                      ymin = Biotic_Vmax_Mean - Biotic_Vmax_SE,
                      ymax = Biotic_Vmax_Mean + Biotic_Vmax_SE),
                  width = 0.2) +
    labs(title = "Biotic Rates (Vmax)",
         x = "YEP Treatment", 
         y = "Fitted Rate (h⁻¹)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save combined plot
  pdf(file.path(out_dir, "Slide6_Rate_Comparison_h_units.pdf"), width = 12, height = 6)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  
  return(list(abiotic_plot = p1, biotic_plot = p2, data = plot_data))
}

# -------------------------------
# Slide 7: Factorial Effects (Interaction Plots)
# -------------------------------
create_slide7_factorial <- function(summary_df, out_dir = "modeling_outputs") {
  
  # Prepare factorial data
  factorial_data <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    mutate(
      Moisture = ifelse(grepl("Dry", Treatment_Type), "Dry", "Wet"),
      DOC_Treatment = case_when(
        grepl("Control", Treatment_Type) ~ "Control",
        grepl("Unburned", Treatment_Type) ~ "Unburned",
        grepl("HighBurn", Treatment_Type) ~ "HighBurn",
        TRUE ~ "Other"
      )
    ) %>%
    group_by(Moisture, DOC_Treatment) %>%
    summarise(
      n = n(),
      Mean_Time_Anoxia = mean(Time_to_Anoxia_h, na.rm = TRUE),
      SE_Time_Anoxia = sd(Time_to_Anoxia_h, na.rm = TRUE) / sqrt(n()),
      Mean_Linear_Rate = mean(Linear_kL_per_h, na.rm = TRUE),
      Mean_Biotic_Rate = mean(Biotic_Vmax_per_h, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Interaction plot for Time to Anoxia
  p1 <- ggplot(factorial_data, aes(x = DOC_Treatment, y = Mean_Time_Anoxia, 
                                   color = Moisture, group = Moisture)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = Mean_Time_Anoxia - SE_Time_Anoxia,
                      ymax = Mean_Time_Anoxia + SE_Time_Anoxia),
                  width = 0.1) +
    labs(title = "Moisture × DOC Amendment Effects",
         subtitle = "Time to Anoxia",
         x = "DOC Treatment",
         y = "Time to Anoxia (h)") +
    theme_minimal() +
    scale_color_manual(values = c("Dry" = "brown", "Wet" = "blue"))
  
  # Main effects comparison
  moisture_effect <- factorial_data %>%
    group_by(Moisture) %>%
    summarise(Mean_Time = mean(Mean_Time_Anoxia), .groups = 'drop')
  
  doc_effect <- factorial_data %>%
    group_by(DOC_Treatment) %>%
    summarise(Mean_Time = mean(Mean_Time_Anoxia), .groups = 'drop')
  
  p2 <- ggplot(moisture_effect, aes(x = Moisture, y = Mean_Time)) +
    geom_col(fill = "lightgreen", alpha = 0.7) +
    labs(title = "Main Effect: Moisture", 
         y = "Mean Time to Anoxia (h)") +
    theme_minimal()
  
  p3 <- ggplot(doc_effect, aes(x = DOC_Treatment, y = Mean_Time)) +
    geom_col(fill = "orange", alpha = 0.7) +
    labs(title = "Main Effect: DOC Treatment",
         y = "Mean Time to Anoxia (h)") +
    theme_minimal()
  
  # Save plots
  pdf(file.path(out_dir, "Slide7_Factorial_Effects.pdf"), width = 15, height = 5)
  grid.arrange(p1, p2, p3, ncol = 3)
  dev.off()
  
  # Calculate effect magnitudes for conclusions
  moisture_range <- diff(range(moisture_effect$Mean_Time))
  doc_range <- diff(range(doc_effect$Mean_Time))
  
  stronger_effect <- ifelse(moisture_range > doc_range, "Moisture", "DOC Treatment")
  
  return(list(
    interaction_plot = p1,
    moisture_effect = moisture_effect,
    doc_effect = doc_effect, 
    stronger_effect = stronger_effect,
    factorial_data = factorial_data
  ))
}

# -------------------------------
# Slide 8: Key Findings Summary
# -------------------------------
generate_slide8_summary <- function(summary_df, factorial_results) {
  
  # Find fastest and slowest treatments
  treatment_summary <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      Mean_Time_Anoxia = mean(Time_to_Anoxia_h, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(Mean_Time_Anoxia)
  
  fastest_treatment <- treatment_summary$Treatment_Type[1]
  slowest_treatment <- treatment_summary$Treatment_Type[nrow(treatment_summary)]
  
  # Model performance
  best_model_summary <- summary_df %>%
    count(Best_Model) %>%
    arrange(desc(n)) %>%
    mutate(Percentage = round(n/sum(n) * 100, 1))
  
  best_model <- best_model_summary$Best_Model[1]
  best_model_pct <- best_model_summary$Percentage[1]
  
  # Process dominance analysis
  process_dominance <- summary_df %>%
    filter(!is.na(Combined_kL_per_h), !is.na(Combined_Vmax_per_h)) %>%
    mutate(
      Rate_Ratio = Combined_kL_per_h / Combined_Vmax_per_h,
      Dominant_Process = ifelse(Rate_Ratio > 1, "Abiotic", "Biotic")
    ) %>%
    group_by(Treatment_Type, Dominant_Process) %>%
    count() %>%
    group_by(Treatment_Type) %>%
    slice_max(n) %>%
    select(Treatment_Type, Dominant_Process)
  
  # Create summary text
  findings_summary <- list(
    fastest_treatment = fastest_treatment,
    slowest_treatment = slowest_treatment,
    best_model = paste0(best_model, " model fits ", best_model_pct, "% of samples"),
    stronger_effect = factorial_results$stronger_effect,
    process_dominance = process_dominance
  )
  
  return(findings_summary)
}

# -------------------------------
# Run all analyses
# -------------------------------
# Assuming you have summary_df from your main analysis

# Generate Slide 6
slide6_results <- create_slide6_figures(summary_df)
print("Generated Slide 6: Rate Comparison Figures")

# Generate Slide 7  
slide7_results <- create_slide7_factorial(summary_df)
print("Generated Slide 7: Factorial Effects")
print(paste("Stronger effect:", slide7_results$stronger_effect))

# Generate Slide 8 summary
slide8_summary <- generate_slide8_summary(summary_df, slide7_results)
print("Generated Slide 8: Key Findings")
print("Key Findings:")
print(paste("- Fastest oxygen consumption:", slide8_summary$fastest_treatment))
print(paste("- Slowest oxygen consumption:", slide8_summary$slowest_treatment))
print(paste("- Best model:", slide8_summary$best_model))
print(paste("- Stronger factorial effect:", slide8_summary$stronger_effect))

# Print process dominance
print("Process Dominance by Treatment:")
print(slide8_summary$process_dominance)