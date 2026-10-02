# ----------------------------------------
# Complete YEP Analysis with Rates in mg/kg/h
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

save_sample_fit_plots <- FALSE
save_development_summary_plots <- FALSE

# -------------------------------
# Load and process data
# -------------------------------
file_path <- file.path(
  "v2_data", "v2_YEP_Sample_Data",
  "YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv"
)
mapping_file_path <- file.path("v2_data", "Merged_Firesting_Mapping_2025-10-24.csv")
mass_volume_file_path <- file.path(
  "v2_data", "v2_YEP_Sample_Data",
  "v2_YEP_Sediment_Water_Mass_Volume.csv"
)

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
# Constants / model assumptions
# -------------------------------
# Primary formulation follows Patel et al. (2024):
# Km is a fixed apparent O2 half-saturation parameter and
# CB starts at 1 as a normalized microbial biomass state.
Km_primary <- 10       # mg O2/L
Km <- Km_primary
CB <- 1.0              # initial normalized microbial biomass (dimensionless)

# Km values used to address uncertainty in transferring the
# soil-derived Patel value to freshwater sediment systems.
Km_sensitivity_values <- c(0.2, 0.5, 1.0, 5.28, 10)

# Optimization bounds
Vmax_lower_bound <- 0.01
Vmax_upper_bound <- 50
kL_lower_bound <- 0.001
kL_upper_bound <- 10
bound_tolerance <- 1e-5

# Practical guardrail for accepting the extra parameter in the
# combined model. This is an effect-size criterion, not a p-value.
minimum_rmse_improvement_pct <- 5

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

# 1. Linear abiotic: dDO/dt = kL * [DO] (kL in h^-1)
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

# 2. Biotic: dDO/dt = Vmax * [DO]/([DO] + Km) * CB (Vmax in h^-1)
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

# 3. Combined: dDO/dt = (Vmax * [DO]/([DO] + Km) * CB) + kL * [DO]
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
# Model diagnostic helper functions
# -------------------------------
calc_rmse <- function(obs, pred) {
  if(is.null(pred) || length(obs) != length(pred)) return(NA_real_)
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

calc_lag1_residual_ac <- function(obs, pred) {
  if(is.null(pred) || length(obs) != length(pred)) return(NA_real_)
  residuals <- obs - pred
  residuals <- residuals[is.finite(residuals)]
  if(length(residuals) < 3 || sd(residuals) == 0) return(NA_real_)
  cor(residuals[-length(residuals)], residuals[-1], use = "complete.obs")
}

at_lower_bound <- function(value, lower_bound, tolerance = bound_tolerance) {
  if(length(value) == 0 || is.na(value) || !is.finite(value)) return(NA)
  value <= (lower_bound + tolerance)
}

at_upper_bound <- function(value, upper_bound, tolerance = bound_tolerance) {
  if(length(value) == 0 || is.na(value) || !is.finite(value)) return(NA)
  value >= (upper_bound - tolerance)
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
  
  # Get conversion factor for this sample.
  # Samples without sediment mass/volume information (e.g., water-only controls)
  # retain concentration-space fits, but sediment-normalized quantities are NA.
  if(!is.null(conversion_info)) {
    dry_mass_g <- conversion_info$Dry_Sediment_Mass_g
    water_vol_L <- conversion_info$Water_Volume_L
    conversion_factor <- conversion_info$Conversion_Factor
  } else {
    dry_mass_g <- NA_real_
    water_vol_L <- NA_real_
    conversion_factor <- NA_real_
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
                      lower = c(kL = kL_lower_bound), upper = c(kL = kL_upper_bound),
                      control = list(maxiter = 500))
  }, silent = TRUE)
  
  # Biotic fit
  try({
    fit_BIO <- modFit(p = c(Vmax = linear_rate_est), f = Objective_BIO,
                      lower = c(Vmax = Vmax_lower_bound), upper = c(Vmax = Vmax_upper_bound),
                      control = list(maxiter = 500))
  }, silent = TRUE)
  
  # Combined fit
  try({
    fit_COM <- modFit(p = c(Vmax = linear_rate_est/2, kL = linear_rate_est/DO0/2),
                      f = Objective_COM,
                      lower = c(Vmax = Vmax_lower_bound, kL = kL_lower_bound),
                      upper = c(Vmax = Vmax_upper_bound, kL = kL_upper_bound),
                      control = list(maxiter = 500))
  }, silent = TRUE)
  
  # Generate predictions
  model_LIN <- model_BIO <- model_COM <- NULL
  
  if(!is.null(fit_LIN)) {
    pars_LIN <- list(kL = fit_LIN$par[["kL"]])
    model_LIN <- solveLinear(pars_LIN, times = Data$time, DO0 = DO0)
    message("  Linear fit: kL = ", round(fit_LIN$par[["kL"]], 4), " h^-1, SSR = ", round(fit_LIN$ssr, 2))
  }
  
  if(!is.null(fit_BIO)) {
    pars_BIO <- list(Vmax = fit_BIO$par[["Vmax"]])
    model_BIO <- solveBiotic(pars_BIO, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    message("  Biotic fit: Vmax = ", round(fit_BIO$par[["Vmax"]], 4), " h^-1, SSR = ", round(fit_BIO$ssr, 2))
  }
  
  if(!is.null(fit_COM)) {
    pars_COM <- list(Vmax = fit_COM$par[["Vmax"]], kL = fit_COM$par[["kL"]])
    model_COM <- solveCombined(pars_COM, times = Data$time, DO0 = DO0, CB = CB, Km = Km)
    message("  Combined fit: Vmax = ", round(fit_COM$par[["Vmax"]], 4), " h^-1",
            ", kL = ", round(fit_COM$par[["kL"]], 4), " h^-1, SSR = ", round(fit_COM$ssr, 2))
  }
  
  # -------------------------------
  # Model diagnostics
  # -------------------------------
  # AIC is retained as a descriptive diagnostic, but is not used alone for
  # model selection because observations occur every 5 s and residuals are
  # temporally autocorrelated.
  Linear_RMSE <- if(!is.null(model_LIN)) calc_rmse(Data$DO, model_LIN$DO) else NA_real_
  Biotic_RMSE <- if(!is.null(model_BIO)) calc_rmse(Data$DO, model_BIO$DO) else NA_real_
  Combined_RMSE <- if(!is.null(model_COM)) calc_rmse(Data$DO, model_COM$DO) else NA_real_
  
  Linear_Residual_Lag1_AC <- if(!is.null(model_LIN)) calc_lag1_residual_ac(Data$DO, model_LIN$DO) else NA_real_
  Biotic_Residual_Lag1_AC <- if(!is.null(model_BIO)) calc_lag1_residual_ac(Data$DO, model_BIO$DO) else NA_real_
  Combined_Residual_Lag1_AC <- if(!is.null(model_COM)) calc_lag1_residual_ac(Data$DO, model_COM$DO) else NA_real_
  
  Combined_kL_at_lower_bound <- if(!is.null(fit_COM)) {
    at_lower_bound(fit_COM$par[["kL"]], kL_lower_bound)
  } else NA
  
  Combined_Vmax_at_lower_bound <- if(!is.null(fit_COM)) {
    at_lower_bound(fit_COM$par[["Vmax"]], Vmax_lower_bound)
  } else NA
  
  Combined_kL_at_upper_bound <- if(!is.null(fit_COM)) {
    at_upper_bound(fit_COM$par[["kL"]], kL_upper_bound)
  } else NA
  
  Combined_Vmax_at_upper_bound <- if(!is.null(fit_COM)) {
    at_upper_bound(fit_COM$par[["Vmax"]], Vmax_upper_bound)
  } else NA
  
  Best_Single_RMSE <- suppressWarnings(min(c(Linear_RMSE, Biotic_RMSE), na.rm = TRUE))
  if(!is.finite(Best_Single_RMSE)) Best_Single_RMSE <- NA_real_
  
  Combined_RMSE_Improvement_vs_BestSingle_pct <-
    if(!is.na(Combined_RMSE) && !is.na(Best_Single_RMSE) && Best_Single_RMSE > 0) {
      100 * (Best_Single_RMSE - Combined_RMSE) / Best_Single_RMSE
    } else {
      NA_real_
    }
  
  # Create individual sample plots only when troubleshooting model fits.
  if (save_sample_fit_plots) {
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
                 labels = paste("kL =", round(fit_LIN$par[["kL"]], 3), "h^-1"),
                 adj = c(0, 1), bg = "white", cex = 0.9)
          } else if(i == 2 && !is.null(fit_BIO)) {
            text(x = max(Data$time) * 0.05, y = max(Data$DO) * 0.8,
                 labels = paste("Vmax =", round(fit_BIO$par[["Vmax"]], 3), "h^-1"),
                 adj = c(0, 1), bg = "white", cex = 0.9)
          } else if(i == 3 && !is.null(fit_COM)) {
            param_text <- paste("Vmax =", round(fit_COM$par[["Vmax"]], 3), "h^-1\n",
                                "kL =", round(fit_COM$par[["kL"]], 3), "h^-1")
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
  }
  
  # Calculate DO consumption rates (mg/L/h) - VOLUMETRIC RATES
  DO_rate_linear <- if(!is.null(fit_LIN)) fit_LIN$par[["kL"]] * DO0 else NA
  DO_rate_biotic <- if(!is.null(fit_BIO)) fit_BIO$par[["Vmax"]] * (DO0/(DO0 + Km)) * CB else NA
  DO_rate_combined <- if(!is.null(fit_COM)) {
    bio_rate <- fit_COM$par[["Vmax"]] * (DO0/(DO0 + Km)) * CB
    abi_rate <- fit_COM$par[["kL"]] * DO0
    bio_rate + abi_rate
  } else NA
  
  # Calculate DO consumption rates per dry sediment mass.
  DO_rate_linear_per_g <- if(!is.null(fit_LIN)) DO_rate_linear * conversion_factor else NA
  DO_rate_biotic_per_g <- if(!is.null(fit_BIO)) DO_rate_biotic * conversion_factor else NA
  DO_rate_combined_per_g <- if(!is.null(fit_COM)) DO_rate_combined * conversion_factor else NA
  
  DO_rate_linear_per_kg <- DO_rate_linear_per_g * 1000
  DO_rate_biotic_per_kg <- DO_rate_biotic_per_g * 1000
  DO_rate_combined_per_kg <- DO_rate_combined_per_g * 1000
  
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
    Km_mg_per_L = Km,
    CB_initial_dimensionless = CB,
    # Linear model (kL in h^-1)
    Linear_kL_per_h = ifelse(!is.null(fit_LIN), round(fit_LIN$par[["kL"]], 4), NA),
    Linear_Rate_mg_per_L_per_h = round(DO_rate_linear, 4),
    Linear_Rate_mg_per_g_per_h = round(DO_rate_linear_per_g, 4),
    Linear_Rate_mg_per_kg_per_h = round(DO_rate_linear_per_kg, 4),
    Linear_SSR = ifelse(!is.null(fit_LIN), round(fit_LIN$ssr, 2), NA),
    Linear_AIC = ifelse(!is.null(fit_LIN), round(compute_info(fit_LIN$ssr, 1, n)$AIC, 2), NA),
    Linear_RMSE_mg_per_L = round(Linear_RMSE, 6),
    Linear_Residual_Lag1_AC = round(Linear_Residual_Lag1_AC, 6),
    # Biotic model (Vmax in h^-1)
    Biotic_Vmax_per_h = ifelse(!is.null(fit_BIO), round(fit_BIO$par[["Vmax"]], 4), NA),
    Biotic_Rate_mg_per_L_per_h = round(DO_rate_biotic, 4),
    Biotic_Rate_mg_per_g_per_h = round(DO_rate_biotic_per_g, 4),
    Biotic_Rate_mg_per_kg_per_h = round(DO_rate_biotic_per_kg, 4),
    Biotic_SSR = ifelse(!is.null(fit_BIO), round(fit_BIO$ssr, 2), NA),
    Biotic_AIC = ifelse(!is.null(fit_BIO), round(compute_info(fit_BIO$ssr, 1, n)$AIC, 2), NA),
    Biotic_RMSE_mg_per_L = round(Biotic_RMSE, 6),
    Biotic_Residual_Lag1_AC = round(Biotic_Residual_Lag1_AC, 6),
    # Combined model
    Combined_Vmax_per_h = ifelse(!is.null(fit_COM), round(fit_COM$par[["Vmax"]], 4), NA),
    Combined_kL_per_h = ifelse(!is.null(fit_COM), round(fit_COM$par[["kL"]], 4), NA),
    Combined_Rate_mg_per_L_per_h = round(DO_rate_combined, 4),
    Combined_Rate_mg_per_g_per_h = round(DO_rate_combined_per_g, 4),
    Combined_Rate_mg_per_kg_per_h = round(DO_rate_combined_per_kg, 4),
    Combined_SSR = ifelse(!is.null(fit_COM), round(fit_COM$ssr, 2), NA),
    Combined_AIC = ifelse(!is.null(fit_COM), round(compute_info(fit_COM$ssr, 2, n)$AIC, 2), NA),
    Combined_RMSE_mg_per_L = round(Combined_RMSE, 6),
    Combined_Residual_Lag1_AC = round(Combined_Residual_Lag1_AC, 6),
    Combined_RMSE_Improvement_vs_BestSingle_pct = round(Combined_RMSE_Improvement_vs_BestSingle_pct, 4),
    Combined_kL_at_lower_bound = Combined_kL_at_lower_bound,
    Combined_Vmax_at_lower_bound = Combined_Vmax_at_lower_bound,
    Combined_kL_at_upper_bound = Combined_kL_at_upper_bound,
    Combined_Vmax_at_upper_bound = Combined_Vmax_at_upper_bound,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------
  # Model interpretation
  # -------------------------------
  # AIC is retained only as a diagnostic. Because the 5-second DO observations
  # are serially autocorrelated, treating every point as an independent
  # observation can make AIC differences look much more decisive than they are.
  aic_values <- c(summary_row$Linear_AIC, summary_row$Biotic_AIC, summary_row$Combined_AIC)
  aic_names <- c("FirstOrder", "Biotic", "Combined")
  valid_aic <- !is.na(aic_values) & is.finite(aic_values)
  
  if(any(valid_aic)) {
    lowest_aic_idx <- which.min(aic_values[valid_aic])
    summary_row$Lowest_AIC_Model <- aic_names[valid_aic][lowest_aic_idx]
    summary_row$Lowest_AIC <- min(aic_values[valid_aic])
  } else {
    summary_row$Lowest_AIC_Model <- "None"
    summary_row$Lowest_AIC <- NA_real_
  }
  
  # Best single-component model by RMSE
  if(!is.na(Linear_RMSE) && !is.na(Biotic_RMSE)) {
    Best_Single_Model_RMSE <- ifelse(Linear_RMSE <= Biotic_RMSE, "FirstOrder", "Biotic")
  } else if(!is.na(Linear_RMSE)) {
    Best_Single_Model_RMSE <- "FirstOrder"
  } else if(!is.na(Biotic_RMSE)) {
    Best_Single_Model_RMSE <- "Biotic"
  } else {
    Best_Single_Model_RMSE <- "None"
  }
  
  summary_row$Best_Single_Model_RMSE <- Best_Single_Model_RMSE
  
  combined_lower_rmse_than_both <-
    !is.na(Combined_RMSE) && !is.na(Linear_RMSE) && !is.na(Biotic_RMSE) &&
    Combined_RMSE < Linear_RMSE && Combined_RMSE < Biotic_RMSE
  
  combined_lower_abs_ac_than_both <-
    !is.na(Combined_Residual_Lag1_AC) &&
    !is.na(Linear_Residual_Lag1_AC) &&
    !is.na(Biotic_Residual_Lag1_AC) &&
    abs(Combined_Residual_Lag1_AC) < abs(Linear_Residual_Lag1_AC) &&
    abs(Combined_Residual_Lag1_AC) < abs(Biotic_Residual_Lag1_AC)
  
  summary_row$Combined_Lower_RMSE_Than_Both_Singles <- combined_lower_rmse_than_both
  summary_row$Combined_Lower_Abs_Lag1_AC_Than_Both_Singles <- combined_lower_abs_ac_than_both
  
  # Interpretation hierarchy:
  # 1) Boundary collapse -> simpler model.
  # 2) Upper-bound estimate -> manual review because the parameter may be poorly constrained.
  # 3) Otherwise call the combined model only if the extra parameter produces a
  #    practically meaningful (>=5%) RMSE improvement AND improves residual
  #    temporal structure relative to both single-component models.
  # 4) If the combined model does not meet that guardrail, retain the lower-RMSE
  #    single-component model.
  if(is.null(fit_COM)) {
    summary_row$Model_Interpretation <- Best_Single_Model_RMSE
    summary_row$Model_Interpretation_Basis <-
      "Combined fit failed; retained lower-RMSE single-component model"
    
  } else if(isTRUE(Combined_kL_at_lower_bound) && isTRUE(Combined_Vmax_at_lower_bound)) {
    summary_row$Model_Interpretation <- "ManualReview"
    summary_row$Model_Interpretation_Basis <-
      "Both combined-model parameters are at their lower bounds"
    
  } else if(isTRUE(Combined_kL_at_upper_bound) || isTRUE(Combined_Vmax_at_upper_bound)) {
    summary_row$Model_Interpretation <- "ManualReview"
    summary_row$Model_Interpretation_Basis <-
      "At least one combined-model parameter is at its upper optimization bound"
    
  } else if(isTRUE(Combined_kL_at_lower_bound)) {
    summary_row$Model_Interpretation <- "Biotic"
    summary_row$Model_Interpretation_Basis <-
      "Combined kL is at its lower bound; the combined model collapses to the biotic model"
    
  } else if(isTRUE(Combined_Vmax_at_lower_bound)) {
    summary_row$Model_Interpretation <- "FirstOrder"
    summary_row$Model_Interpretation_Basis <-
      "Combined Vmax is at its lower bound; the combined model collapses to the first-order model"
    
  } else if(
    isTRUE(combined_lower_rmse_than_both) &&
    isTRUE(combined_lower_abs_ac_than_both) &&
    !is.na(Combined_RMSE_Improvement_vs_BestSingle_pct) &&
    Combined_RMSE_Improvement_vs_BestSingle_pct >= minimum_rmse_improvement_pct
  ) {
    summary_row$Model_Interpretation <- "Combined"
    summary_row$Model_Interpretation_Basis <- paste0(
      "Both parameters are interior; combined model improves RMSE by at least ",
      minimum_rmse_improvement_pct,
      "% and has lower absolute lag-1 residual autocorrelation than both single models"
    )
    
  } else {
    summary_row$Model_Interpretation <- Best_Single_Model_RMSE
    summary_row$Model_Interpretation_Basis <- paste0(
      "Combined model did not meet the ",
      minimum_rmse_improvement_pct,
      "% RMSE-improvement plus residual-structure criterion; retained lower-RMSE single-component model"
    )
  }
  
  return(list(summary = summary_row,
              fits = list(linear = fit_LIN, biotic = fit_BIO, combined = fit_COM),
              data = Data))
}

# -------------------------------
# Create comprehensive evaluation plots - Enhanced for kg units
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
  
  # 1. Model Performance Summary
  # RMSE is used for across-sample comparison because it is in DO units and
  # does not increase simply because a time series contains more observations.
  pdf(file.path(out_dir, "Model_Performance_Summary.pdf"), width = 12, height = 8)
  
  rmse_summary <- plot_data %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      FirstOrder_RMSE_Mean = mean(Linear_RMSE_mg_per_L, na.rm = TRUE),
      Biotic_RMSE_Mean = mean(Biotic_RMSE_mg_per_L, na.rm = TRUE),
      Combined_RMSE_Mean = mean(Combined_RMSE_mg_per_L, na.rm = TRUE),
      .groups = 'drop'
    )
  
  rmse_matrix <- as.matrix(rmse_summary[, c("FirstOrder_RMSE_Mean", "Biotic_RMSE_Mean", "Combined_RMSE_Mean")])
  rownames(rmse_matrix) <- rmse_summary$Treatment_Type
  
  barplot(t(rmse_matrix), beside = TRUE,
          col = c("lightblue", "lightcoral", "lightgreen"),
          main = "Model Performance by YEP Treatment\n(Lower RMSE = Better Fit)",
          ylab = "RMSE (mg O2/L)",
          legend.text = c("First-order", "Biotic", "Combined"),
          args.legend = list(x = "topright", bty = "n"),
          las = 2)
  
  dev.off()
  
  if (save_development_summary_plots) {
    # 2. Rate Comparison (h^-1 units like original)
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
                   ylab = "Fitted Rate (h^-1)",
                   ylim = c(0, max(rate_data$Linear_kL_Mean + rate_data$Linear_kL_SE, na.rm = TRUE) * 1.1))
    
    arrows(bp1, rate_data$Linear_kL_Mean - rate_data$Linear_kL_SE,
           bp1, rate_data$Linear_kL_Mean + rate_data$Linear_kL_SE,
           length = 0.05, angle = 90, code = 3)
    
    # Biotic rates (Vmax)
    bp2 <- barplot(rate_data$Biotic_Vmax_Mean,
                   names.arg = gsub("_", "\n", rate_data$Treatment_Type),
                   col = "lightcoral",
                   main = "Biotic Rates (Vmax)",
                   ylab = "Fitted Rate (h^-1)",
                   ylim = c(0, max(rate_data$Biotic_Vmax_Mean + rate_data$Biotic_Vmax_SE, na.rm = TRUE) * 1.1))
    
    arrows(bp2, rate_data$Biotic_Vmax_Mean - rate_data$Biotic_Vmax_SE,
           bp2, rate_data$Biotic_Vmax_Mean + rate_data$Biotic_Vmax_SE,
           length = 0.05, angle = 90, code = 3)
    
    dev.off()
    
    # 3. NEW: Mass-Normalized Rate Comparison (mg/kg/h units)
    pdf(file.path(out_dir, "Rate_Comparison_Mass_Normalized_kg.pdf"), width = 12, height = 6)
    
    mass_rate_data <- plot_data %>%
      group_by(Treatment_Type) %>%
      summarise(
        n = n(),
        Linear_Rate_Mean = mean(Linear_Rate_mg_per_kg_per_h, na.rm = TRUE),
        Linear_Rate_SE = sd(Linear_Rate_mg_per_kg_per_h, na.rm = TRUE) / sqrt(n()),
        Biotic_Rate_Mean = mean(Biotic_Rate_mg_per_kg_per_h, na.rm = TRUE),
        Biotic_Rate_SE = sd(Biotic_Rate_mg_per_kg_per_h, na.rm = TRUE) / sqrt(n()),
        Combined_Rate_Mean = mean(Combined_Rate_mg_per_kg_per_h, na.rm = TRUE),
        Combined_Rate_SE = sd(Combined_Rate_mg_per_kg_per_h, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    par(mfrow = c(1, 3))
    
    # Linear rates per kg
    bp4 <- barplot(mass_rate_data$Linear_Rate_Mean,
                   names.arg = gsub("_", "\n", mass_rate_data$Treatment_Type),
                   col = "lightblue",
                   main = "Linear Abiotic Rates\n(Mass Normalized)",
                   ylab = "Rate (mg O2/kg dry sediment/h)",
                   ylim = c(0, max(mass_rate_data$Linear_Rate_Mean + mass_rate_data$Linear_Rate_SE, na.rm = TRUE) * 1.1))
    
    arrows(bp4, mass_rate_data$Linear_Rate_Mean - mass_rate_data$Linear_Rate_SE,
           bp4, mass_rate_data$Linear_Rate_Mean + mass_rate_data$Linear_Rate_SE,
           length = 0.05, angle = 90, code = 3)
    
    # Biotic rates per kg
    bp5 <- barplot(mass_rate_data$Biotic_Rate_Mean,
                   names.arg = gsub("_", "\n", mass_rate_data$Treatment_Type),
                   col = "lightcoral",
                   main = "Biotic Rates\n(Mass Normalized)",
                   ylab = "Rate (mg O2/kg dry sediment/h)",
                   ylim = c(0, max(mass_rate_data$Biotic_Rate_Mean + mass_rate_data$Biotic_Rate_SE, na.rm = TRUE) * 1.1))
    
    arrows(bp5, mass_rate_data$Biotic_Rate_Mean - mass_rate_data$Biotic_Rate_SE,
           bp5, mass_rate_data$Biotic_Rate_Mean + mass_rate_data$Biotic_Rate_SE,
           length = 0.05, angle = 90, code = 3)
    
    # Combined rates per kg
    bp6 <- barplot(mass_rate_data$Combined_Rate_Mean,
                   names.arg = gsub("_", "\n", mass_rate_data$Treatment_Type),
                   col = "lightgreen",
                   main = "Combined Rates\n(Mass Normalized)",
                   ylab = "Rate (mg O2/kg dry sediment/h)",
                   ylim = c(0, max(mass_rate_data$Combined_Rate_Mean + mass_rate_data$Combined_Rate_SE, na.rm = TRUE) * 1.1))
    
    arrows(bp6, mass_rate_data$Combined_Rate_Mean - mass_rate_data$Combined_Rate_SE,
           bp6, mass_rate_data$Combined_Rate_Mean + mass_rate_data$Combined_Rate_SE,
           length = 0.05, angle = 90, code = 3)
    
    dev.off()
    
    # 4. Time to Anoxia Summary (like original Table 2)
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
    
  }
  
  message("Created model performance summary plot")
}

# -------------------------------
# Create model-comparison summary table
# -------------------------------
create_summary_tables <- function(summary_df) {
  table2 <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    group_by(Treatment_Type) %>%
    summarise(
      n = n(),
      FirstOrder_RMSE = round(mean(Linear_RMSE_mg_per_L, na.rm = TRUE), 4),
      Biotic_RMSE = round(mean(Biotic_RMSE_mg_per_L, na.rm = TRUE), 4),
      Combined_RMSE = round(mean(Combined_RMSE_mg_per_L, na.rm = TRUE), 4),
      FirstOrder_Mean_Abs_Lag1_AC = round(mean(abs(Linear_Residual_Lag1_AC), na.rm = TRUE), 4),
      Biotic_Mean_Abs_Lag1_AC = round(mean(abs(Biotic_Residual_Lag1_AC), na.rm = TRUE), 4),
      Combined_Mean_Abs_Lag1_AC = round(mean(abs(Combined_Residual_Lag1_AC), na.rm = TRUE), 4),
      Combined_kL_at_lower_bound_n = sum(Combined_kL_at_lower_bound %in% TRUE),
      Combined_Vmax_at_lower_bound_n = sum(Combined_Vmax_at_lower_bound %in% TRUE),
      Most_Common_Model_Interpretation = {
        x <- Model_Interpretation[!is.na(Model_Interpretation)]
        if(length(x) == 0) NA_character_ else names(which.max(table(x)))
      },
      .groups = 'drop'
    )
  
  write.csv(table2, file.path(out_dir, "YEP_Model_Performance_Table.csv"), row.names = FALSE)
  message("Created model performance summary table")
  return(list(model_performance = table2))
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

message("\n=== YEP Analysis  ===")
message("Samples: ", length(samples))
message("Output directory: ", out_dir)
message("Trimming: First and last 2 minutes removed")
message("Mass/Volume normalization: Sample-specific from CSV")
message("Rate units: mg O2/kg dry sediment/h")

# Run analysis for all samples
all_results <- lapply(samples, function(sname) {
  dat_s <- df2 %>% filter(Sample_Name == sname)
  treatment_type <- unique(dat_s$Treatment_Type)[1]
  
  # Get conversion info for this sample
  conversion_info <- conversion_factors %>% filter(Sample_Name == sname)
  if(nrow(conversion_info) == 0) {
    conversion_info <- NULL
    message("No sediment mass/volume data for ", sname, "; sediment-normalized outputs will be NA")
  } else {
    conversion_info <- conversion_info[1, ]
  }
  
  fit_all_yep_models(dat_s, sname, treatment_type, conversion_info)
})

# Combine results
summary_df <- bind_rows(lapply(all_results, function(x) if(!is.null(x)) x$summary else NULL))

# Save and analyze results
if(nrow(summary_df) > 0) {
  enhanced_file <- file.path(out_dir, "YEP_Complete_Modeling_Outputs.csv")
  
  write.csv(summary_df, enhanced_file, row.names = FALSE)
  
  # Create all evaluation plots
  create_evaluation_plots(summary_df)
  
  # Create summary tables for PowerPoint
  tables <- create_summary_tables(summary_df)
  
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
      Mean_Rate_mg_L_h = round(mean(Linear_Rate_mg_per_L_per_h, na.rm = TRUE), 3),
      Mean_Rate_mg_kg_h = round(mean(Linear_Rate_mg_per_kg_per_h, na.rm = TRUE), 1),
      .groups = 'drop'
    )
  
  print("YEP Treatment Results:")
  print(treatment_results)
  
  # Model interpretation (exclude water-only/Unknown controls)
  model_preference <- summary_df %>%
    filter(Treatment_Type != "Unknown") %>%
    count(Model_Interpretation, sort = TRUE) %>%
    mutate(Percentage = round(n / sum(n) * 100, 1))
  
  print("\nModel interpretation:")
  print(model_preference)
  
  # Key findings
  fastest_treatment <- treatment_results$Treatment_Type[which.min(treatment_results$Mean_Time_h)]
  slowest_treatment <- treatment_results$Treatment_Type[which.max(treatment_results$Mean_Time_h)]
  most_common_model <- model_preference$Model_Interpretation[1]
  
  message("\n=== KEY FINDINGS ===")
  message("Fastest oxygen consumption: ", fastest_treatment)
  message("Slowest oxygen consumption: ", slowest_treatment)
  message("Most common model interpretation: ", most_common_model, " (", model_preference$Percentage[1], "% of sediment samples)")
  
} else {
  message("No successful analyses completed. Check data and sample naming.")
}



# ============================================================
# Km SENSITIVITY ANALYSIS
# ============================================================
# Primary Km = 10 mg O2/L following the apparent-parameter
# formulation used by Patel et al. (2024). Lower values test
# whether conclusions are robust to freshwater-sediment-relevant
# half-saturation assumptions. Km is not selected by lowest AIC.
# ============================================================

message("\n=== Km Sensitivity Analysis ===")

Km_original <- Km
save_sample_fit_plots_original <- save_sample_fit_plots
save_development_summary_plots_original <- save_development_summary_plots

save_sample_fit_plots <- FALSE
save_development_summary_plots <- FALSE

sensitivity_list <- list()

for(km_test in Km_sensitivity_values) {
  message("\nRunning Km = ", km_test, " mg/L")
  Km <- km_test
  
  km_results <- lapply(samples, function(sname) {
    dat_s <- df2 %>% filter(Sample_Name == sname)
    treatment_type <- unique(dat_s$Treatment_Type)[1]
    
    conversion_info <- conversion_factors %>% filter(Sample_Name == sname)
    if(nrow(conversion_info) == 0) {
      conversion_info <- NULL
    } else {
      conversion_info <- conversion_info[1, ]
    }
    
    fit_all_yep_models(dat_s, sname, treatment_type, conversion_info)
  })
  
  km_summary <- bind_rows(
    lapply(km_results, function(x) if(!is.null(x)) x$summary else NULL)
  )
  
  if(nrow(km_summary) > 0) {
    sensitivity_list[[as.character(km_test)]] <- km_summary
  }
}

# Restore primary settings
Km <- Km_original
save_sample_fit_plots <- save_sample_fit_plots_original
save_development_summary_plots <- save_development_summary_plots_original

Km_sensitivity_df <- bind_rows(sensitivity_list)

write.csv(
  Km_sensitivity_df,
  file.path(out_dir, "YEP_Km_Sensitivity_All_Samples.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Comparison to the primary Patel-style Km = 10 mg/L
# ------------------------------------------------------------
sensitivity_reference <- Km_sensitivity_df %>%
  filter(Km_mg_per_L == Km_primary) %>%
  select(
    Sample_Name,
    Reference_Combined_Vmax = Combined_Vmax_per_h,
    Reference_Combined_kL = Combined_kL_per_h,
    Reference_Combined_Rate = Combined_Rate_mg_per_kg_per_h,
    Reference_Model_Interpretation = Model_Interpretation
  )

sensitivity_comparison <- Km_sensitivity_df %>%
  left_join(sensitivity_reference, by = "Sample_Name") %>%
  mutate(
    Vmax_Percent_Change_From_Primary = ifelse(
      !is.na(Reference_Combined_Vmax) & Reference_Combined_Vmax != 0,
      100 * (Combined_Vmax_per_h - Reference_Combined_Vmax) / Reference_Combined_Vmax,
      NA_real_
    ),
    kL_Percent_Change_From_Primary = ifelse(
      !is.na(Reference_Combined_kL) & Reference_Combined_kL != 0,
      100 * (Combined_kL_per_h - Reference_Combined_kL) / Reference_Combined_kL,
      NA_real_
    ),
    Combined_Rate_Percent_Change_From_Primary = ifelse(
      !is.na(Reference_Combined_Rate) & Reference_Combined_Rate != 0,
      100 * (Combined_Rate_mg_per_kg_per_h - Reference_Combined_Rate) / Reference_Combined_Rate,
      NA_real_
    ),
    Model_Interpretation_Changed = Model_Interpretation != Reference_Model_Interpretation
  )

write.csv(
  sensitivity_comparison,
  file.path(out_dir, "YEP_Km_Sensitivity_Comparison_to_10.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Treatment-level sensitivity summary
# ------------------------------------------------------------
Km_sensitivity_treatment <- Km_sensitivity_df %>%
  filter(Treatment_Type != "Unknown") %>%
  group_by(Km_mg_per_L, Treatment_Type) %>%
  summarise(
    n = n(),
    Mean_Combined_Vmax = mean(Combined_Vmax_per_h, na.rm = TRUE),
    SD_Combined_Vmax = sd(Combined_Vmax_per_h, na.rm = TRUE),
    Mean_Combined_kL = mean(Combined_kL_per_h, na.rm = TRUE),
    SD_Combined_kL = sd(Combined_kL_per_h, na.rm = TRUE),
    Mean_Combined_Rate_mg_kg_h = mean(Combined_Rate_mg_per_kg_per_h, na.rm = TRUE),
    Mean_Combined_RMSE = mean(Combined_RMSE_mg_per_L, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  Km_sensitivity_treatment,
  file.path(out_dir, "YEP_Km_Sensitivity_Treatment_Summary.csv"),
  row.names = FALSE
)

# Model interpretation across Km values; exclude water-only controls
sensitivity_model_preference <- Km_sensitivity_df %>%
  filter(Treatment_Type != "Unknown") %>%
  count(Km_mg_per_L, Model_Interpretation) %>%
  group_by(Km_mg_per_L) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  sensitivity_model_preference,
  file.path(out_dir, "YEP_Km_Sensitivity_Model_Preference.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Sensitivity plots
# ------------------------------------------------------------
p_vmax <- ggplot(
  Km_sensitivity_treatment,
  aes(x = Km_mg_per_L, y = Mean_Combined_Vmax, group = Treatment_Type, linetype = Treatment_Type)
) +
  geom_line() +
  geom_point() +
  scale_x_log10(breaks = Km_sensitivity_values) +
  labs(
    x = expression(K[m]~"(mg O"[2]~L^{-1}*")"),
    y = expression("Mean combined-model " * V[max]~"(h"^{-1}*")"),
    title = "Sensitivity of fitted Vmax to fixed Km"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "YEP_Km_Sensitivity_Vmax.pdf"), p_vmax, width = 8, height = 6)

p_kl <- ggplot(
  Km_sensitivity_treatment,
  aes(x = Km_mg_per_L, y = Mean_Combined_kL, group = Treatment_Type, linetype = Treatment_Type)
) +
  geom_line() +
  geom_point() +
  scale_x_log10(breaks = Km_sensitivity_values) +
  labs(
    x = expression(K[m]~"(mg O"[2]~L^{-1}*")"),
    y = expression("Mean combined-model " * k[L]~"(h"^{-1}*")"),
    title = "Sensitivity of fitted kL to fixed Km"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "YEP_Km_Sensitivity_kL.pdf"), p_kl, width = 8, height = 6)

p_total <- ggplot(
  Km_sensitivity_treatment,
  aes(x = Km_mg_per_L, y = Mean_Combined_Rate_mg_kg_h, group = Treatment_Type, linetype = Treatment_Type)
) +
  geom_line() +
  geom_point() +
  scale_x_log10(breaks = Km_sensitivity_values) +
  labs(
    x = expression(K[m]~"(mg O"[2]~L^{-1}*")"),
    y = expression("Initial combined O"[2]~"consumption (mg O"[2]~kg^{-1}~h^{-1}*")"),
    title = "Sensitivity of total modeled O2 consumption to Km"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "YEP_Km_Sensitivity_Total_Rate.pdf"), p_total, width = 8, height = 6)

message("\n=== Km sensitivity analysis complete ===")
print(Km_sensitivity_treatment)
print(sensitivity_model_preference)
