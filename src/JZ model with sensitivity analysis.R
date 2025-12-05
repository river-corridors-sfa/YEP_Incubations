# ----------------------------------------
# Complete YEP Time-to-Anoxia Analysis Following Original Paper Methodology
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
  message("Created output directory: ", out_dir)
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
# Constants from original paper
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
  incubation_samples <- mass_volume_df %>%
    filter(grepl("_INC-", Sample_Name)) %>%
    mutate(
      Dry_Sediment_Mass_g = as.numeric(Dry_Sediment_Mass_g),
      Water_Mass_g = as.numeric(Water_Mass_g),
      Water_Volume_L = Water_Mass_g / 1000,  # Convert g to L
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
# Model solver functions following original paper equations
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
  if(SSR <= 0 || n <= k) {
    return(list(AIC = Inf, BIC = Inf))
  }
  AIC <- n * log(SSR / n) + 2 * k
  BIC <- n * log(SSR / n) + k * log(n)
  list(AIC = AIC, BIC = BIC)
}

# -------------------------------
# Complete fitting function - Original paper style
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
    # Use YEP defaults based on treatment
    if(grepl("YEP1", sample_name)) {
      dry_mass_g <- 10.2
      water_vol_L <- 0.0304
      conversion_factor <- 0.0298
    } else {
      dry_mass_g <- 14.1
      water_vol_L <- 0.0366
      conversion_factor <- 0.0259
    }
  }
  
  message("Fitting YEP models for ", sample_name, " (n=", n, ", treatment=", treatment_type, ")")
  message("  Initial DO: ", round(DO0, 2), " mg/L, Time to anoxia: ", round(time_to_anoxia_obs, 1), " h")
  message("  Dry mass: ", round(dry_mass_g, 2), " g, Water vol: ", round(water_vol_L*1000, 1), " mL")
  
  # Estimate starting parameters
  linear_rate_est <- abs(coef(lm(DO ~ time, data = Data))[["time"]])
  
  # Define objective functions
  
  # 1. Linear abiotic objective
  Objective_LIN <- function(x) {
    if(x[["kL"]] <= 0 || x[["kL"]] > 50) return(1e10)
    pars <- list(kL = x[["kL"]])
    out <- solveLinear(pars, times = Data$time, DO0 = DO0)
    if(any(is.na(out$DO))) return(1e10)
    tryCatch(modCost(model = out, obs = Data), error = function(e) 1e10)
  }
  
  # 2. Biotic objective
  Objective_BIO <- function(x) {
    if(x[["Vmax"]] <= 0 || x[["Vmax"]] > 100) return(1e10)
    pars <- list(Vmax = x[["Vmax"]])
    out <- solveBiotic(pars, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    if(any(is.na(out$DO))) return(1e10)
    tryCatch(modCost(model = out, obs = Data), error = function(e) 1e10)
  }
  
  # 3. Combined objective
  Objective_COM <- function(x) {
    if(x[["Vmax"]] <= 0 || x[["Vmax"]] > 100 || x[["kL"]] <= 0 || x[["kL"]] > 50) return(1e10)
    pars <- list(Vmax = x[["Vmax"]], kL = x[["kL"]])
    out <- solveCombined(pars, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    if(any(is.na(out$DO))) return(1e10)
    tryCatch(modCost(model = out, obs = Data), error = function(e) 1e10)
  }
  
  # Fit all models
  fit_LIN <- fit_BIO <- fit_COM <- NULL
  
  # Linear fit
  try({
    fit_LIN <- modFit(p = c(kL = linear_rate_est/DO0), f = Objective_LIN,
                      lower = c(kL = 0.001), upper = c(kL = 10),
                      control = list(maxiter = 500))
  }, silent = TRUE)
  
  # Biotic fit
  try({
    fit_BIO <- modFit(p = c(Vmax = linear_rate_est), f = Objective_BIO,
                      lower = c(Vmax = 0.01), upper = c(Vmax = 50),
                      control = list(maxiter = 500))
  }, silent = TRUE)
  
  # Combined fit
  try({
    fit_COM <- modFit(p = c(Vmax = linear_rate_est/2, kL = linear_rate_est/DO0/2), 
                      f = Objective_COM,
                      lower = c(Vmax = 0.01, kL = 0.001), 
                      upper = c(Vmax = 50, kL = 10),
                      control = list(maxiter = 500))
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
  
  # Create individual sample plots (original paper style)
  tryCatch({
    pdf_file <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", sample_name), ".pdf"))
    pdf(pdf_file, width = 14, height = 10)
    
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    
    # Plot 1: All models together (PowerPoint style)
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
    
    text(x = max(Data$time) * 0.6, y = max(Data$DO) * 0.8,
         labels = ssr_text, adj = c(0, 1), bg = "white", cex = 0.9, family = "mono")
    
    legend("topright", c("Observed", "Linear Abiotic", "Nonlinear Biotic", "Combined"),
           pch = c(16, NA, NA, NA), lty = c(NA, 1, 2, 3),
           col = c("orange", "blue", "red", "green"), lwd = c(NA, 3, 3, 3), bty = "n")
    
    # Individual model plots (original paper style)
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
        plot(1, 1, type = "n", xlab = "", ylab = "", 
             main = paste(models_list[[i]]$title, "- Fit Failed"))
        text(1, 1, "Model fit failed", cex = 1.5, col = "red")
      }
    }
    
    dev.off()
    message("  Wrote: ", pdf_file)
    
  }, error = function(e) {
    message("  Error creating plots: ", e$message)
  })
  
  # Calculate DO consumption rates (mg/L/h)
  DO_rate_linear <- if(!is.null(fit_LIN)) fit_LIN$par[["kL"]] * DO0 else NA
  DO_rate_biotic <- if(!is.null(fit_BIO)) fit_BIO$par[["Vmax"]] * (DO0/(DO0 + Km)) * CB else NA
  DO_rate_combined <- if(!is.null(fit_COM)) {
    bio_rate <- fit_COM$par[["Vmax"]] * (DO0/(DO0 + Km)) * CB
    abi_rate <- fit_COM$par[["kL"]] * DO0
    bio_rate + abi_rate
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
    Linear_Rate_mg_per_L_per_h = round(DO_rate_linear, 4),
    Linear_SSR = ifelse(!is.null(fit_LIN), round(fit_LIN$ssr, 2), NA),
    Linear_AIC = ifelse(!is.null(fit_LIN), round(compute_info(fit_LIN$ssr, 1, n)$AIC, 2), NA),
    # Biotic model (Vmax in h⁻¹)
    Biotic_Vmax_per_h = ifelse(!is.null(fit_BIO), round(fit_BIO$par[["Vmax"]], 4), NA),
    Biotic_Rate_mg_per_L_per_h = round(DO_rate_biotic, 4),
    Biotic_SSR = ifelse(!is.null(fit_BIO), round(fit_BIO$ssr, 2), NA),
    Biotic_AIC = ifelse(!is.null(fit_BIO), round(compute_info(fit_BIO$ssr, 1, n)$AIC, 2), NA),
    # Combined model
    Combined_Vmax_per_h = ifelse(!is.null(fit_COM), round(fit_COM$par[["Vmax"]], 4), NA),
    Combined_kL_per_h = ifelse(!is.null(fit_COM), round(fit_COM$par[["kL"]], 4), NA),
    Combined_Rate_mg_per_L_per_h = round(DO_rate_combined, 4),
    Combined_SSR = ifelse(!is.null(fit_COM), round(fit_COM$ssr, 2), NA),
    Combined_AIC = ifelse(!is.null(fit_COM), round(compute_info(fit_COM$ssr, 2, n)$AIC, 2), NA),
    stringsAsFactors = FALSE
  )
  
  # Best model by AIC
  aic_values <- c(summary_row$Linear_AIC, summary_row$Biotic_AIC, summary_row$Combined_AIC)
  aic_names <- c("Linear", "Biotic", "Combined")
  valid_aic <- !is.na(aic_values) & is.finite(aic_values)
  
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
# Create comprehensive evaluation plots
# -------------------------------
create_evaluation_plots <- function(summary_df) {
  
  # Filter valid data
  plot_data <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    mutate(
      Moisture = ifelse(grepl("Dry", Treatment_Type), "Dry", "Wet"),
      DOC_Treatment = case_when(
        grepl("Control", Treatment_Type) ~ "Control",
        grepl("Unburned", Treatment_Type) ~ "Unburned",
        grepl("HighBurn", Treatment_Type) ~ "HighBurn",
        TRUE ~ "Other"
      )
    )
  
  # 1. Model Performance Summary (SSR comparison like original paper)
  pdf(file.path(out_dir, "Model_Performance_Summary.pdf"), width = 12, height = 8)
  
  ssr_summary <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Linear_SSR_Mean = mean(Linear_SSR, na.rm = TRUE),
      Biotic_SSR_Mean = mean(Biotic_SSR, na.rm = TRUE),
      Combined_SSR_Mean = mean(Combined_SSR, na.rm = TRUE),
      .groups = 'drop'
    )
  
  ssr_matrix <- as.matrix(ssr_summary[, c("Linear_SSR_Mean", "Biotic_SSR_Mean", "Combined_SSR_Mean")])
  rownames(ssr_matrix) <- ssr_summary$Treatment_Type
  
  barplot(t(ssr_matrix), beside = TRUE,
          col = c("lightblue", "lightcoral", "lightgreen"),
          main = "Model Performance (SSR) by YEP Treatment\n(Lower = Better Fit)",
          ylab = "Sum of Squares Regression (SSR)",
          legend.text = c("Linear Abiotic", "Nonlinear Biotic", "Combined"),
          args.legend = list(x = "topright", bty = "n"),
          las = 2)
  
  dev.off()
  
  # 2. Rate Comparison (h⁻¹ units like original)
  pdf(file.path(out_dir, "Rate_Comparison_h_units.pdf"), width = 12, height = 6)
  
  rate_data <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Linear_kL_Mean = mean(Linear_kL_per_h, na.rm = TRUE),
      Linear_kL_SE = sd(Linear_kL_per_h, na.rm = TRUE) / sqrt(n()),
      Biotic_Vmax_Mean = mean(Biotic_Vmax_per_h, na.rm = TRUE),
      Biotic_Vmax_SE = sd(Biotic_Vmax_per_h, na.rm = TRUE) / sqrt(n()),
      Combined_kL_Mean = mean(Combined_kL_per_h, na.rm = TRUE),
      Combined_Vmax_Mean = mean(Combined_Vmax_per_h, na.rm = TRUE),
      .groups = 'drop'
    )
  
  par(mfrow = c(1, 2))
  
  # Abiotic rates (kL)
  bp1 <- barplot(rate_data$Linear_kL_Mean, 
                 names.arg = gsub("_", "\n", rate_data$Treatment_Type),
                 col = "lightblue",
                 main = "Abiotic Rates (kL)", 
                 ylab = "Fitted Rate (h⁻¹)",
                 ylim = c(0, max(rate_data$Linear_kL_Mean + rate_data$Linear_kL_SE, na.rm = TRUE) * 1.1))
  
  arrows(bp1, rate_data$Linear_kL_Mean - rate_data$Linear_kL_SE,
         bp1, rate_data$Linear_kL_Mean + rate_data$Linear_kL_SE,
         length = 0.05, angle = 90, code = 3)
  
  # Biotic rates (Vmax)
  bp2 <- barplot(rate_data$Biotic_Vmax_Mean,
                 names.arg = gsub("_", "\n", rate_data$Treatment_Type), 
                 col = "lightcoral",
                 main = "Biotic Rates (Vmax)",
                 ylab = "Fitted Rate (h⁻¹)",
                 ylim = c(0, max(rate_data$Biotic_Vmax_Mean + rate_data$Biotic_Vmax_SE, na.rm = TRUE) * 1.1))
  
  arrows(bp2, rate_data$Biotic_Vmax_Mean - rate_data$Biotic_Vmax_SE,
         bp2, rate_data$Biotic_Vmax_Mean + rate_data$Biotic_Vmax_SE,
         length = 0.05, angle = 90, code = 3)
  
  dev.off()
  
  # 3. Time to Anoxia Summary (like original Table 2)
  pdf(file.path(out_dir, "Time_to_Anoxia_Summary.pdf"), width = 10, height = 6)
  
  anoxia_summary <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Mean_Time_h = mean(Time_to_Anoxia_h, na.rm = TRUE),
      SE_Time_h = sd(Time_to_Anoxia_h, na.rm = TRUE) / sqrt(n()),
      Mean_Rate = mean(Linear_Rate_mg_per_L_per_h, na.rm = TRUE),
      SE_Rate = sd(Linear_Rate_mg_per_L_per_h, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  bp3 <- barplot(anoxia_summary$Mean_Time_h,
                 names.arg = gsub("_", "\n", anoxia_summary$Treatment_Type),
                 col = rainbow(nrow(anoxia_summary)),
                 main = "Time to Anoxia by YEP Treatment",
                 ylab = "Time to Anoxia (h)",
                 ylim = c(0, max(anoxia_summary$Mean_Time_h + anoxia_summary$SE_Time_h, na.rm = TRUE) * 1.1))
  
  arrows(bp3, anoxia_summary$Mean_Time_h - anoxia_summary$SE_Time_h,
         bp3, anoxia_summary$Mean_Time_h + anoxia_summary$SE_Time_h,
         length = 0.05, angle = 90, code = 3)
  
  # Add sample sizes
  text(bp3, anoxia_summary$Mean_Time_h + anoxia_summary$SE_Time_h + max(anoxia_summary$Mean_Time_h) * 0.05,
       paste("n =", anoxia_summary$n), cex = 0.8)
  
  dev.off()
  
  # 4. Factorial Effects (Interaction Plots)
  pdf(file.path(out_dir, "Factorial_Effects_YEP.pdf"), width = 12, height = 8)
  
  factorial_data <- plot_data %>%
    group_by(Moisture, DOC_Treatment) %>%
    summarise(
      n = n(),
      Mean_Time_Anoxia = mean(Time_to_Anoxia_h, na.rm = TRUE),
      SE_Time_Anoxia = sd(Time_to_Anoxia_h, na.rm = TRUE) / sqrt(n()),
      Mean_Linear_Rate = mean(Linear_kL_per_h, na.rm = TRUE),
      Mean_Biotic_Rate = mean(Biotic_Vmax_per_h, na.rm = TRUE),
      .groups = 'drop'
    )
  
  par(mfrow = c(2, 2))
  
  # Interaction plot for time to anoxia
  interaction.plot(x.factor = factorial_data$DOC_Treatment,
                   trace.factor = factorial_data$Moisture,
                   response = factorial_data$Mean_Time_Anoxia,
                   type = "b", pch = c(16, 17), col = c("blue", "red"),
                   lwd = 2, main = "Moisture × DOC Treatment Effects\n(Time to Anoxia)",
                   xlab = "DOC Treatment", ylab = "Time to Anoxia (h)")
  
  # Interaction plot for linear rates  
  interaction.plot(x.factor = factorial_data$DOC_Treatment,
                   trace.factor = factorial_data$Moisture,
                   response = factorial_data$Mean_Linear_Rate,
                   type = "b", pch = c(16, 17), col = c("darkblue", "darkred"),
                   lwd = 2, main = "Moisture × DOC Treatment Effects\n(Linear Rates)",
                   xlab = "DOC Treatment", ylab = "kL (h⁻¹)")
  
  # Main effect: Moisture
  moisture_effect <- factorial_data %>%
    group_by(Moisture) %>%
    summarise(Mean_Time = mean(Mean_Time_Anoxia), .groups = 'drop')
  
  barplot(moisture_effect$Mean_Time, names.arg = moisture_effect$Moisture,
          col = c("brown", "lightblue"), main = "Main Effect: Moisture",
          ylab = "Mean Time to Anoxia (h)")
  
  # Main effect: DOC Treatment
  doc_effect <- factorial_data %>%
    group_by(DOC_Treatment) %>%
    summarise(Mean_Time = mean(Mean_Time_Anoxia), .groups = 'drop')
  
  barplot(doc_effect$Mean_Time, names.arg = doc_effect$DOC_Treatment,
          col = c("orange", "yellow", "pink"), main = "Main Effect: DOC Treatment",
          ylab = "Mean Time to Anoxia (h)")
  
  dev.off()
  
  # 5. Rate Ratio Sensitivity (like original Figure 6A)
  pdf(file.path(out_dir, "Rate_Ratio_Sensitivity_YEP.pdf"), width = 10, height = 6)
  
  # Calculate actual rate ratios from YEP data
  actual_ratios <- plot_data %>%
    filter(!is.na(Combined_kL_per_h) & !is.na(Combined_Vmax_per_h) & Combined_Vmax_per_h > 0) %>%
    mutate(Rate_Ratio = Combined_kL_per_h / Combined_Vmax_per_h) %>%
    group_by(Treatment_Type) %>%
    summarise(Mean_Rate_Ratio = mean(Rate_Ratio, na.rm = TRUE), .groups = 'drop')
  
  # Simulate DO curves for different rate ratios
  rate_ratios <- c(0, 0.1, 0.5, 1, 2, 5)
  time_seq <- seq(0, 30, by = 0.5)
  colors <- rainbow(length(rate_ratios))
  
  plot(0, 0, type = "n", xlim = c(0, 30), ylim = c(0, 10),
       xlab = "Time (h)", ylab = "DO (mg/L)", 
       main = "Rate Ratio Sensitivity Analysis\n(kL/Vmax = Abiotic/Biotic)")
  
  for(i in 1:length(rate_ratios)) {
    kL_val <- rate_ratios[i] * 0.1
    vmax_val <- 0.5
    initial_DO <- 8.5
    
    bio_rate <- vmax_val * (initial_DO/(initial_DO + Km)) * CB
    abi_rate <- kL_val * initial_DO
    total_rate <- (bio_rate + abi_rate) / initial_DO
    
    do_curve <- initial_DO * exp(-total_rate * time_seq)
    do_curve <- pmax(0, do_curve)
    
    lines(time_seq, do_curve, col = colors[i], lwd = 3)
  }
  
  legend("topright", paste("kL/Vmax =", rate_ratios), 
         col = colors, lwd = 3, bty = "n", cex = 0.9)
  
  # Add YEP actual rate ratios text
  if(nrow(actual_ratios) > 0) {
    actual_text <- "YEP Rate Ratios:\n"
    for(i in 1:nrow(actual_ratios)) {
      actual_text <- paste0(actual_text, 
                            gsub("_", " ", actual_ratios$Treatment_Type[i]), 
                            ": ", round(actual_ratios$Mean_Rate_Ratio[i], 2), "\n")
    }
    
    text(x = 25, y = 8, labels = actual_text, adj = c(1, 1), 
         bg = "white", cex = 0.8, family = "mono")
  }
  
  dev.off()
  
  message("Created all evaluation plots")
}

# -------------------------------
# Create summary tables for PowerPoint
# -------------------------------
create_summary_tables <- function(summary_df) {
  
  # Table 1: Time to Anoxia Results (like original Table 2)
  table1 <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Time_to_Anoxia_h = paste0(round(mean(Time_to_Anoxia_h, na.rm = TRUE), 1), " ± ", 
                                round(sd(Time_to_Anoxia_h, na.rm = TRUE), 1)),
      DO_Consumption_Rate = paste0(round(mean(Linear_Rate_mg_per_L_per_h, na.rm = TRUE), 2), " ± ",
                                   round(sd(Linear_Rate_mg_per_L_per_h, na.rm = TRUE), 2)),
      Linear_kL = paste0(round(mean(Linear_kL_per_h, na.rm = TRUE), 3), " ± ",
                         round(sd(Linear_kL_per_h, na.rm = TRUE), 3)),
      Biotic_Vmax = paste0(round(mean(Biotic_Vmax_per_h, na.rm = TRUE), 3), " ± ",
                           round(sd(Biotic_Vmax_per_h, na.rm = TRUE), 3)),
      .groups = 'drop'
    )
  
  write.csv(table1, file.path(out_dir, "YEP_TimeToAnoxia_Summary_Table.csv"), row.names = FALSE)
  
  # Table 2: Model Performance (SSR comparison like original)
  table2 <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      Linear_SSR = round(mean(Linear_SSR, na.rm = TRUE), 1),
      Biotic_SSR = round(mean(Biotic_SSR, na.rm = TRUE), 1),
      Combined_SSR = round(mean(Combined_SSR, na.rm = TRUE), 1),
      Best_Model = names(sort(table(Best_Model)))[length(names(sort(table(Best_Model))))],
      .groups = 'drop'
    )
  
  write.csv(table2, file.path(out_dir, "YEP_Model_Performance_Table.csv"), row.names = FALSE)
  
  message("Created summary tables")
  return(list(time_to_anoxia = table1, model_performance = table2))
}

# -------------------------------
# Main analysis execution
# -------------------------------

# Process data with trimming and treatment classification
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
print(head(conversion_factors))

samples <- unique(df2$Sample_Name)

message("\n=== YEP Time-to-Anoxia Analysis (Following Original Paper) ===")
message("Samples: ", length(samples))
message("Output directory: ", out_dir)
message("Trimming: First and last 2 minutes removed")
message("Mass/Volume normalization: Sample-specific from CSV")

# Run analysis for all samples
all_results <- lapply(samples, function(sname) {
  dat_s <- df2 %>% filter(Sample_Name == sname)
  treatment_type <- unique(dat_s$Treatment_Type)[1]
  
  # Get conversion info for this sample
  conversion_info <- conversion_factors %>% filter(Sample_Name == sname)
  if(nrow(conversion_info) == 0) {
    conversion_info <- NULL
    message("Warning: No mass/volume data found for ", sname, ", using treatment defaults")
  } else {
    conversion_info <- conversion_info[1, ]
  }
  
  fit_all_yep_models(dat_s, sname, treatment_type, conversion_info)
})

# Combine results
summary_df <- bind_rows(lapply(all_results, function(x) if(!is.null(x)) x$summary else NULL))

# Save and analyze results
if(nrow(summary_df) > 0) {
  summary_file <- file.path(out_dir, "YEP_Complete_Analysis_Following_Paper.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  message("\nResults saved: ", summary_file)
  
  # Create all evaluation plots
  create_evaluation_plots(summary_df)
  
  # Create summary tables for PowerPoint
  tables <- create_summary_tables(summary_df)
  
  # Print comprehensive analysis summary
  message("\n=== YEP ANALYSIS SUMMARY (Following Original Paper) ===")
  
  # Treatment summary
  treatment_results <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      Mean_Time_h = round(mean(Time_to_Anoxia_h, na.rm = TRUE), 1),
      SD_Time_h = round(sd(Time_to_Anoxia_h, na.rm = TRUE), 1),
      Mean_Linear_SSR = round(mean(Linear_SSR, na.rm = TRUE), 1),
      Mean_Biotic_SSR = round(mean(Biotic_SSR, na.rm = TRUE), 1),
      Mean_Combined_SSR = round(mean(Combined_SSR, na.rm = TRUE), 1),
      .groups = 'drop'
    )
  
  print("YEP Treatment Results:")
  print(treatment_results)
  
  # Model preference
  model_preference <- summary_df %>%
    count(Best_Model, sort = TRUE) %>%
    mutate(Percentage = round(n / sum(n) * 100, 1))
  
  print("\nModel Selection (by AIC):")
  print(model_preference)
  
  # Key findings
  fastest_treatment <- treatment_results$Treatment_Type[which.min(treatment_results$Mean_Time_h)]
  slowest_treatment <- treatment_results$Treatment_Type[which.max(treatment_results$Mean_Time_h)]
  best_model_overall <- model_preference$Best_Model[1]
  
  message("\n=== KEY FINDINGS ===")
  message("Fastest oxygen consumption: ", fastest_treatment)
  message("Slowest oxygen consumption: ", slowest_treatment) 
  message("Best model overall: ", best_model_overall, " (", model_preference$Percentage[1], "% of samples)")
  
} else {
  message("No successful analyses completed. Check data and sample naming.")
}

message("\n=== Analysis Complete! ===")
message("All outputs saved in: ", out_dir, "/")
message("\nFiles for evaluation:")
message("- Individual sample plots: [SampleName].pdf")
message("- Complete results: YEP_Complete_Analysis_Following_Paper.csv") 
message("- Model performance: Model_Performance_Summary.pdf")
message("- Rate comparison: Rate_Comparison_h_units.pdf")
message("- Time to anoxia: Time_to_Anoxia_Summary.pdf")
message("- Factorial effects: Factorial_Effects_YEP.pdf")
message("- Sensitivity analysis: Rate_Ratio_Sensitivity_YEP.pdf")
message("- Summary tables: YEP_TimeToAnoxia_Summary_Table.csv, YEP_Model_Performance_Table.csv")