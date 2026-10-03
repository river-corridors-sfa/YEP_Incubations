# ============================================================
# YEP dissolved oxygen modeling analysis
# Patel-style model fitting + Km sensitivity analysis
# ============================================================
#
# Primary modeling choices follow Patel et al. (2024):
#   1. Fit three model configurations separately by least squares:
#        - first-order (chemical/abiotic-like)
#        - biological Michaelis-Menten
#        - combined first-order + biological
#   2. Evaluate model fit using SSR (sum of squared residuals).
#   3. Use Km = 10 mg O2/L as the primary fixed apparent
#      half-saturation constant, as in Patel et al. (2024).
#   4. Initialize normalized active microbial biomass at CB = 1
#      and allow it to evolve during the simulation.
#
# A Km sensitivity analysis is included to address whether the
# primary soil-derived Patel value changes conclusions for these
# freshwater sediment incubations.
# ============================================================

rm(list = ls(all = TRUE))

library(deSolve)
library(FME)
library(ggplot2)
library(dplyr)
library(tidyr)


# ============================================================
# Output settings
# ============================================================

out_dir <- "modeling_outputs"

if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save_sample_fit_plots <- FALSE


# ============================================================
# Load data
# ============================================================

file_path <- file.path(
  "v2_data",
  "v2_YEP_Sample_Data",
  "YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv"
)

mapping_file_path <- file.path(
  "v2_data",
  "Merged_Firesting_Mapping_2025-10-24.csv"
)

mass_volume_file_path <- file.path(
  "v2_data",
  "v2_YEP_Sample_Data",
  "v2_YEP_Sediment_Water_Mass_Volume.csv"
)


mass_volume_df <- read.csv(
  mass_volume_file_path,
  stringsAsFactors = FALSE,
  skip = 2
) %>%
  filter(grepl("YEP", Sample_Name))


mapping_df <- read.csv(
  mapping_file_path,
  stringsAsFactors = FALSE
)


df <- read.csv(
  file_path,
  skip = 14,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(df) <- c(
  "Field_Name",
  "Sample_Name",
  "IGSN",
  "Material",
  "DateTime",
  "Elapsed_Seconds",
  "Temperature_degreesC",
  "DO_mg_per_L",
  "Firesting_Serial_Number",
  "Methods_Deviation"
)


df <- df %>%
  inner_join(
    mapping_df %>%
      select(
        Sample_Name,
        Start_Time_Elapsed_Seconds = Elapsed_Seconds_Start,
        End_Time_Elapsed_Seconds = Elapsed_Seconds_End
      ),
    by = "Sample_Name"
  ) %>%
  filter(
    Elapsed_Seconds >= Start_Time_Elapsed_Seconds,
    Elapsed_Seconds <= End_Time_Elapsed_Seconds
  ) %>%
  select(
    -Start_Time_Elapsed_Seconds,
    -End_Time_Elapsed_Seconds
  ) %>%
  mutate(
    Elapsed_Seconds = as.numeric(Elapsed_Seconds),
    DO_mg_per_L = as.numeric(DO_mg_per_L)
  ) %>%
  filter(
    !is.na(Elapsed_Seconds),
    !is.na(DO_mg_per_L),
    DO_mg_per_L > 0
  )


# ============================================================
# Fixed model assumptions
# ============================================================

# Primary value from Patel et al. (2024).
# Patel treated Km as a fixed apparent half-saturation parameter.
Km_primary <- 10       # mg O2/L
Km <- Km_primary

# Initial normalized active microbial biomass.
# CB is initialized at 1 and then evolves during the ODE solution.
CB <- 1.0

# Sensitivity values spanning freshwater sediment / hyporheic
# literature values and the Patel primary value.
Km_sensitivity_values <- c(
  0.2,
  0.5,
  1.0,
  5.28,
  10
)

# Optimization bounds retained from the original YEP analysis.
Vmax_lower_bound <- 0.01
Vmax_upper_bound <- 50
kL_lower_bound <- 0.001
kL_upper_bound <- 10


# ============================================================
# Treatment classification
# ============================================================

classify_yep_treatment <- function(sample_names) {
  case_when(
    grepl("YEP1.*S", sample_names, ignore.case = TRUE) ~ "Dry_Control",
    grepl("YEP1.*U", sample_names, ignore.case = TRUE) ~ "Dry_Unburned_DOC",
    grepl("YEP1.*H", sample_names, ignore.case = TRUE) ~ "Dry_HighBurn_DOC",
    grepl("YEP2.*S", sample_names, ignore.case = TRUE) ~ "Wet_Control",
    grepl("YEP2.*U", sample_names, ignore.case = TRUE) ~ "Wet_Unburned_DOC",
    grepl("YEP2.*H", sample_names, ignore.case = TRUE) ~ "Wet_HighBurn_DOC",
    TRUE ~ "Unknown"
  )
}


# ============================================================
# Sample-specific mass / volume conversion factors
# ============================================================

calculate_conversion_factors <- function(mass_volume_df) {
  mass_volume_df %>%
    filter(grepl("_INC-", Sample_Name)) %>%
    mutate(
      Dry_Sediment_Mass_g = as.numeric(Dry_Sediment_Mass_g),
      Water_Mass_g = as.numeric(Water_Mass_g),
      Water_Volume_L = Water_Mass_g / 1000,
      Conversion_Factor = Water_Volume_L / Dry_Sediment_Mass_g
    ) %>%
    select(
      Sample_Name,
      Dry_Sediment_Mass_g,
      Water_Volume_L,
      Conversion_Factor
    )
}

conversion_factors <- calculate_conversion_factors(mass_volume_df)


# ============================================================
# Trim first and last 2 minutes
# ============================================================

trim_time_series <- function(data_subset) {
  trim_hours <- 2 / 60
  
  data_subset %>%
    filter(
      time_hr >= min(time_hr) + trim_hours,
      time_hr <= max(time_hr) - trim_hours
    )
}


# ============================================================
# Model solver functions
# ============================================================

# ------------------------------------------------------------
# 1. First-order model
#
# dDO/dt = -kL * DO
# ------------------------------------------------------------

solveLinear <- function(pars, times, DO0) {
  
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      dDO <- -kL * DO
      list(c(dDO))
    })
  }
  
  state <- c(DO = DO0)
  
  tryCatch({
    out <- ode(
      y = state,
      times = times,
      func = derivs,
      parms = pars,
      method = "lsoda"
    )
    
    result <- as.data.frame(out[, c("time", "DO")])
    
    if(any(is.na(result$DO)) || any(result$DO < 0)) {
      stop("Invalid first-order solution")
    }
    
    result
    
  }, error = function(e) {
    data.frame(
      time = times,
      DO = NA_real_
    )
  })
}


# ------------------------------------------------------------
# 2. Biological Michaelis-Menten model
#
# r1 = Vmax * DO / (DO + Km)
# dDO/dt = -r1 * CB
# dCB/dt =  r1 * CB
#
# CB starts at 1 and evolves during the reaction, following
# the Patel model formulation.
# ------------------------------------------------------------

solveBiotic <- function(pars, times, DO0, CB, Km) {
  
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      
      r1 <- Vmax * DO / (DO + Km + 1e-10)
      
      dDO <- -r1 * CB1
      dCB1 <- r1 * CB1
      
      list(c(dDO, dCB1))
    })
  }
  
  state <- c(
    DO = DO0,
    CB1 = CB
  )
  
  tryCatch({
    out <- ode(
      y = state,
      times = times,
      func = derivs,
      parms = pars,
      method = "lsoda"
    )
    
    result <- as.data.frame(out[, c("time", "DO", "CB1")])
    
    if(any(is.na(result$DO)) || any(result$DO < 0)) {
      stop("Invalid biological solution")
    }
    
    result
    
  }, error = function(e) {
    data.frame(
      time = times,
      DO = NA_real_,
      CB1 = NA_real_
    )
  })
}


# ------------------------------------------------------------
# 3. Combined model
#
# r1 = Vmax * DO / (DO + Km)
# dDO/dt = -r1 * CB - kL * DO
# dCB/dt =  r1 * CB
# ------------------------------------------------------------

solveCombined <- function(pars, times, DO0, CB, Km) {
  
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      
      r1 <- Vmax * DO / (DO + Km + 1e-10)
      
      dDO <- -r1 * CB1 - kL * DO
      dCB1 <- r1 * CB1
      
      list(c(dDO, dCB1))
    })
  }
  
  state <- c(
    DO = DO0,
    CB1 = CB
  )
  
  tryCatch({
    out <- ode(
      y = state,
      times = times,
      func = derivs,
      parms = pars,
      method = "lsoda"
    )
    
    result <- as.data.frame(out[, c("time", "DO", "CB1")])
    
    if(any(is.na(result$DO)) || any(result$DO < 0)) {
      stop("Invalid combined solution")
    }
    
    result
    
  }, error = function(e) {
    data.frame(
      time = times,
      DO = NA_real_,
      CB1 = NA_real_
    )
  })
}


# ============================================================
# Helper functions
# ============================================================

calculate_time_to_anoxia <- function(time_vec, DO_vec, threshold = 0.1) {
  anoxic_idx <- which(DO_vec <= threshold)[1]
  
  if(is.na(anoxic_idx)) {
    return(max(time_vec))
  }
  
  time_vec[anoxic_idx]
}


# ============================================================
# Fit all three Patel model configurations
# ============================================================

fit_all_yep_models <- function(
    dat,
    sample_name,
    treatment_type = "Unknown",
    conversion_info = NULL
) {
  
  # ----------------------------------------------------------
  # Prepare observations
  # ----------------------------------------------------------
  
  dat_trimmed <- trim_time_series(dat)
  
  Data <- dat_trimmed %>%
    arrange(time_hr) %>%
    select(
      time = time_hr,
      DO = DO_mg_per_L
    ) %>%
    group_by(time) %>%
    summarise(
      DO = mean(DO, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      DO > 0,
      !is.na(DO),
      !is.na(time)
    ) %>%
    arrange(time) %>%
    as.data.frame()
  
  n <- nrow(Data)
  
  if(n < 5) {
    message(
      "Skipping ", sample_name,
      ": insufficient data after trimming (n = ", n, ")"
    )
    return(NULL)
  }
  
  Data$time <- Data$time - min(Data$time)
  
  DO0 <- Data$DO[1]
  time_to_anoxia_obs <- calculate_time_to_anoxia(Data$time, Data$DO)
  
  
  # ----------------------------------------------------------
  # Sample-specific water and sediment quantities
  # ----------------------------------------------------------
  
  if(!is.null(conversion_info)) {
    dry_mass_g <- conversion_info$Dry_Sediment_Mass_g
    water_vol_L <- conversion_info$Water_Volume_L
    conversion_factor <- conversion_info$Conversion_Factor
  } else {
    # Water-only controls have no sediment normalization.
    dry_mass_g <- NA_real_
    water_vol_L <- NA_real_
    conversion_factor <- NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Starting estimates
  # ----------------------------------------------------------
  
  linear_rate_est <- abs(
    coef(
      lm(DO ~ time, data = Data)
    )[["time"]]
  )
  
  if(
    is.na(linear_rate_est) ||
    !is.finite(linear_rate_est) ||
    linear_rate_est <= 0
  ) {
    linear_rate_est <- 0.1
  }
  
  start_kL <- linear_rate_est / DO0
  start_kL <- max(
    kL_lower_bound,
    min(start_kL, kL_upper_bound * 0.9)
  )
  
  start_Vmax <- max(
    Vmax_lower_bound,
    min(linear_rate_est, Vmax_upper_bound * 0.9)
  )
  
  
  # ----------------------------------------------------------
  # Objective functions: least-squares fitting using FME
  # ----------------------------------------------------------
  
  Objective_LIN <- function(x) {
    
    pars <- list(
      kL = x[["kL"]]
    )
    
    out <- solveLinear(
      pars,
      times = Data$time,
      DO0 = DO0
    )
    
    if(any(is.na(out$DO))) {
      return(1e10)
    }
    
    tryCatch(
      modCost(
        model = out[, c("time", "DO")],
        obs = Data
      ),
      error = function(e) 1e10
    )
  }
  
  
  Objective_BIO <- function(x) {
    
    pars <- list(
      Vmax = x[["Vmax"]]
    )
    
    out <- solveBiotic(
      pars,
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
    
    if(any(is.na(out$DO))) {
      return(1e10)
    }
    
    tryCatch(
      modCost(
        model = out[, c("time", "DO")],
        obs = Data
      ),
      error = function(e) 1e10
    )
  }
  
  
  Objective_COM <- function(x) {
    
    pars <- list(
      Vmax = x[["Vmax"]],
      kL = x[["kL"]]
    )
    
    out <- solveCombined(
      pars,
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
    
    if(any(is.na(out$DO))) {
      return(1e10)
    }
    
    tryCatch(
      modCost(
        model = out[, c("time", "DO")],
        obs = Data
      ),
      error = function(e) 1e10
    )
  }
  
  
  # ----------------------------------------------------------
  # Fit each model configuration separately
  # ----------------------------------------------------------
  
  fit_LIN <- NULL
  fit_BIO <- NULL
  fit_COM <- NULL
  
  
  try({
    fit_LIN <- modFit(
      p = c(kL = start_kL),
      f = Objective_LIN,
      lower = c(kL = kL_lower_bound),
      upper = c(kL = kL_upper_bound),
      control = list(maxiter = 500)
    )
  }, silent = TRUE)
  
  
  try({
    fit_BIO <- modFit(
      p = c(Vmax = start_Vmax),
      f = Objective_BIO,
      lower = c(Vmax = Vmax_lower_bound),
      upper = c(Vmax = Vmax_upper_bound),
      control = list(maxiter = 500)
    )
  }, silent = TRUE)
  
  
  try({
    fit_COM <- modFit(
      p = c(
        Vmax = max(start_Vmax / 2, Vmax_lower_bound),
        kL = max(start_kL / 2, kL_lower_bound)
      ),
      f = Objective_COM,
      lower = c(
        Vmax = Vmax_lower_bound,
        kL = kL_lower_bound
      ),
      upper = c(
        Vmax = Vmax_upper_bound,
        kL = kL_upper_bound
      ),
      control = list(maxiter = 500)
    )
  }, silent = TRUE)
  
  
  # ----------------------------------------------------------
  # Predictions
  # ----------------------------------------------------------
  
  model_LIN <- NULL
  model_BIO <- NULL
  model_COM <- NULL
  
  
  if(!is.null(fit_LIN)) {
    model_LIN <- solveLinear(
      list(kL = fit_LIN$par[["kL"]]),
      times = Data$time,
      DO0 = DO0
    )
  }
  
  
  if(!is.null(fit_BIO)) {
    model_BIO <- solveBiotic(
      list(Vmax = fit_BIO$par[["Vmax"]]),
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
  }
  
  
  if(!is.null(fit_COM)) {
    model_COM <- solveCombined(
      list(
        Vmax = fit_COM$par[["Vmax"]],
        kL = fit_COM$par[["kL"]]
      ),
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
  }
  
  
  # ----------------------------------------------------------
  # Initial oxygen-consumption rates
  # ----------------------------------------------------------
  
  DO_rate_linear <- if(!is.null(fit_LIN)) {
    fit_LIN$par[["kL"]] * DO0
  } else {
    NA_real_
  }
  
  
  DO_rate_biotic <- if(!is.null(fit_BIO)) {
    fit_BIO$par[["Vmax"]] *
      (DO0 / (DO0 + Km)) *
      CB
  } else {
    NA_real_
  }
  
  
  bio_rate_combined <- if(!is.null(fit_COM)) {
    fit_COM$par[["Vmax"]] *
      (DO0 / (DO0 + Km)) *
      CB
  } else {
    NA_real_
  }
  
  
  first_order_rate_combined <- if(!is.null(fit_COM)) {
    fit_COM$par[["kL"]] * DO0
  } else {
    NA_real_
  }
  
  
  DO_rate_combined <- if(!is.null(fit_COM)) {
    bio_rate_combined + first_order_rate_combined
  } else {
    NA_real_
  }
  
  
  # Relative contributions at initial DO
  Combined_Biotic_Fraction_Initial <- if(
    !is.na(DO_rate_combined) && DO_rate_combined > 0
  ) {
    bio_rate_combined / DO_rate_combined
  } else {
    NA_real_
  }
  
  
  Combined_FirstOrder_Fraction_Initial <- if(
    !is.na(DO_rate_combined) && DO_rate_combined > 0
  ) {
    first_order_rate_combined / DO_rate_combined
  } else {
    NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Mass-normalized rates
  # ----------------------------------------------------------
  
  DO_rate_linear_per_g <- DO_rate_linear * conversion_factor
  DO_rate_biotic_per_g <- DO_rate_biotic * conversion_factor
  DO_rate_combined_per_g <- DO_rate_combined * conversion_factor
  
  bio_rate_combined_per_g <- bio_rate_combined * conversion_factor
  first_order_rate_combined_per_g <- first_order_rate_combined * conversion_factor
  
  DO_rate_linear_per_kg <- DO_rate_linear_per_g * 1000
  DO_rate_biotic_per_kg <- DO_rate_biotic_per_g * 1000
  DO_rate_combined_per_kg <- DO_rate_combined_per_g * 1000
  
  bio_rate_combined_per_kg <- bio_rate_combined_per_g * 1000
  first_order_rate_combined_per_kg <- first_order_rate_combined_per_g * 1000
  
  
  # ----------------------------------------------------------
  # Optional sample fit plots
  # ----------------------------------------------------------
  
  if(save_sample_fit_plots) {
    
    pdf_file <- file.path(
      out_dir,
      paste0(
        "fit_",
        gsub("[^A-Za-z0-9_-]", "_", sample_name),
        ".pdf"
      )
    )
    
    pdf(pdf_file, width = 8, height = 6)
    
    plot(
      Data$time,
      Data$DO,
      pch = 16,
      xlab = "Time (h)",
      ylab = "DO (mg/L)",
      main = paste(sample_name, "Km =", Km, "mg/L")
    )
    
    if(!is.null(model_LIN)) {
      lines(model_LIN$time, model_LIN$DO, lwd = 2, lty = 1)
    }
    
    if(!is.null(model_BIO)) {
      lines(model_BIO$time, model_BIO$DO, lwd = 2, lty = 2)
    }
    
    if(!is.null(model_COM)) {
      lines(model_COM$time, model_COM$DO, lwd = 2, lty = 3)
    }
    
    legend(
      "topright",
      legend = c(
        "Observed",
        "First-order",
        "Biological",
        "Combined"
      ),
      pch = c(16, NA, NA, NA),
      lty = c(NA, 1, 2, 3),
      bty = "n"
    )
    
    dev.off()
  }
  
  
  # ----------------------------------------------------------
  # Summary row
  # ----------------------------------------------------------
  
  summary_row <- data.frame(
    
    Sample_Name = sample_name,
    Treatment_Type = treatment_type,
    n_points = n,
    Initial_DO_mg_per_L = round(DO0, 4),
    Time_to_Anoxia_h = round(time_to_anoxia_obs, 4),
    
    Km_mg_per_L = Km,
    CB_initial_dimensionless = CB,
    
    Dry_Sediment_Mass_g = round(dry_mass_g, 4),
    Water_Volume_mL = round(water_vol_L * 1000, 4),
    Conversion_Factor = round(conversion_factor, 6),
    
    # First-order model
    Linear_kL_per_h = ifelse(
      !is.null(fit_LIN),
      round(fit_LIN$par[["kL"]], 6),
      NA
    ),
    
    Linear_Rate_mg_per_L_per_h = round(DO_rate_linear, 6),
    Linear_Rate_mg_per_g_per_h = round(DO_rate_linear_per_g, 6),
    Linear_Rate_mg_per_kg_per_h = round(DO_rate_linear_per_kg, 4),
    
    Linear_SSR = ifelse(
      !is.null(fit_LIN),
      round(fit_LIN$ssr, 6),
      NA
    ),
    
    # Biological model
    Biotic_Vmax_per_h = ifelse(
      !is.null(fit_BIO),
      round(fit_BIO$par[["Vmax"]], 6),
      NA
    ),
    
    Biotic_Rate_mg_per_L_per_h = round(DO_rate_biotic, 6),
    Biotic_Rate_mg_per_g_per_h = round(DO_rate_biotic_per_g, 6),
    Biotic_Rate_mg_per_kg_per_h = round(DO_rate_biotic_per_kg, 4),
    
    Biotic_SSR = ifelse(
      !is.null(fit_BIO),
      round(fit_BIO$ssr, 6),
      NA
    ),
    
    # Combined model
    Combined_Vmax_per_h = ifelse(
      !is.null(fit_COM),
      round(fit_COM$par[["Vmax"]], 6),
      NA
    ),
    
    Combined_kL_per_h = ifelse(
      !is.null(fit_COM),
      round(fit_COM$par[["kL"]], 6),
      NA
    ),
    
    Combined_Biotic_Rate_mg_per_L_per_h = round(bio_rate_combined, 6),
    Combined_FirstOrder_Rate_mg_per_L_per_h = round(first_order_rate_combined, 6),
    Combined_Rate_mg_per_L_per_h = round(DO_rate_combined, 6),
    
    Combined_Biotic_Rate_mg_per_kg_per_h = round(bio_rate_combined_per_kg, 4),
    Combined_FirstOrder_Rate_mg_per_kg_per_h = round(first_order_rate_combined_per_kg, 4),
    Combined_Rate_mg_per_g_per_h = round(DO_rate_combined_per_g, 6),
    Combined_Rate_mg_per_kg_per_h = round(DO_rate_combined_per_kg, 4),
    
    Combined_Biotic_Fraction_Initial = round(Combined_Biotic_Fraction_Initial, 4),
    Combined_FirstOrder_Fraction_Initial = round(Combined_FirstOrder_Fraction_Initial, 4),
    
    Combined_SSR = ifelse(
      !is.null(fit_COM),
      round(fit_COM$ssr, 6),
      NA
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Best-fitting configuration following the Patel approach:
  # model fitness is evaluated using SSR.
  # ----------------------------------------------------------
  
  ssr_values <- c(
    summary_row$Linear_SSR,
    summary_row$Biotic_SSR,
    summary_row$Combined_SSR
  )
  
  ssr_names <- c(
    "FirstOrder",
    "Biotic",
    "Combined"
  )
  
  valid_ssr <- !is.na(ssr_values) & is.finite(ssr_values)
  
  if(any(valid_ssr)) {
    best_idx <- which.min(ssr_values[valid_ssr])
    
    summary_row$Best_Fit_Model <-
      ssr_names[valid_ssr][best_idx]
    
    summary_row$Best_SSR <-
      min(ssr_values[valid_ssr])
  } else {
    summary_row$Best_Fit_Model <- "None"
    summary_row$Best_SSR <- NA_real_
  }
  
  
  return(
    list(
      summary = summary_row,
      fits = list(
        first_order = fit_LIN,
        biotic = fit_BIO,
        combined = fit_COM
      ),
      data = Data
    )
  )
}


# ============================================================
# Prepare data for model fitting
# ============================================================

df2 <- df %>%
  group_by(Sample_Name) %>%
  arrange(Elapsed_Seconds, .by_group = TRUE) %>%
  mutate(
    time_hr =
      (Elapsed_Seconds - min(Elapsed_Seconds)) /
      3600
  ) %>%
  ungroup() %>%
  mutate(
    Treatment_Type = classify_yep_treatment(Sample_Name)
  )


samples <- unique(df2$Sample_Name)


# ============================================================
# Function to run all samples at the current Km
# ============================================================

run_all_samples <- function() {
  
  results <- lapply(
    samples,
    function(sname) {
      
      dat_s <- df2 %>%
        filter(Sample_Name == sname)
      
      treatment_type <- unique(dat_s$Treatment_Type)[1]
      
      conversion_info <- conversion_factors %>%
        filter(Sample_Name == sname)
      
      if(nrow(conversion_info) == 0) {
        conversion_info <- NULL
      } else {
        conversion_info <- conversion_info[1, ]
      }
      
      fit_all_yep_models(
        dat = dat_s,
        sample_name = sname,
        treatment_type = treatment_type,
        conversion_info = conversion_info
      )
    }
  )
  
  bind_rows(
    lapply(
      results,
      function(x) {
        if(!is.null(x)) x$summary else NULL
      }
    )
  )
}


# ============================================================
# PRIMARY ANALYSIS
# Km = 10 mg O2/L, following Patel et al. (2024)
# ============================================================

Km <- Km_primary

message("\n============================================")
message("PRIMARY YEP ANALYSIS")
message("Km = ", Km, " mg O2/L")
message("Initial CB = ", CB)
message("Model fitness evaluated using SSR")
message("============================================")


summary_df <- run_all_samples()


write.csv(
  summary_df,
  file.path(
    out_dir,
    "YEP_Complete_Modeling_Outputs.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Primary model performance summary
# ------------------------------------------------------------

model_performance <- summary_df %>%
  filter(Treatment_Type != "Unknown") %>%
  group_by(Treatment_Type) %>%
  summarise(
    n = n(),
    FirstOrder_SSR = mean(Linear_SSR, na.rm = TRUE),
    Biotic_SSR = mean(Biotic_SSR, na.rm = TRUE),
    Combined_SSR = mean(Combined_SSR, na.rm = TRUE),
    Most_Common_Best_Fit_Model = {
      x <- Best_Fit_Model[!is.na(Best_Fit_Model)]
      if(length(x) == 0) {
        NA_character_
      } else {
        names(which.max(table(x)))
      }
    },
    .groups = "drop"
  )

write.csv(
  model_performance,
  file.path(
    out_dir,
    "YEP_Model_Performance_Table.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Primary treatment summary
# ------------------------------------------------------------

primary_treatment_summary <- summary_df %>%
  filter(Treatment_Type != "Unknown") %>%
  group_by(Treatment_Type) %>%
  summarise(
    n = n(),
    Mean_Combined_Vmax = mean(Combined_Vmax_per_h, na.rm = TRUE),
    SD_Combined_Vmax = sd(Combined_Vmax_per_h, na.rm = TRUE),
    Mean_Combined_kL = mean(Combined_kL_per_h, na.rm = TRUE),
    SD_Combined_kL = sd(Combined_kL_per_h, na.rm = TRUE),
    Mean_Combined_Rate_mg_kg_h = mean(
      Combined_Rate_mg_per_kg_per_h,
      na.rm = TRUE
    ),
    Mean_Biotic_Fraction_Initial = mean(
      Combined_Biotic_Fraction_Initial,
      na.rm = TRUE
    ),
    Mean_FirstOrder_Fraction_Initial = mean(
      Combined_FirstOrder_Fraction_Initial,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  primary_treatment_summary,
  file.path(
    out_dir,
    "YEP_Treatment_Model_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Primary best-fit model counts
# ------------------------------------------------------------

primary_model_preference <- summary_df %>%
  filter(Treatment_Type != "Unknown") %>%
  count(Best_Fit_Model, sort = TRUE) %>%
  mutate(
    Percentage = 100 * n / sum(n)
  )

write.csv(
  primary_model_preference,
  file.path(
    out_dir,
    "YEP_Best_Fit_Model_Preference.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Primary model performance plot
# ------------------------------------------------------------

plot_performance <- model_performance %>%
  select(
    Treatment_Type,
    FirstOrder_SSR,
    Biotic_SSR,
    Combined_SSR
  ) %>%
  pivot_longer(
    cols = c(
      FirstOrder_SSR,
      Biotic_SSR,
      Combined_SSR
    ),
    names_to = "Model",
    values_to = "Mean_SSR"
  )

p_performance <- ggplot(
  plot_performance,
  aes(
    x = Treatment_Type,
    y = Mean_SSR,
    fill = Model
  )
) +
  geom_col(position = "dodge") +
  labs(
    x = NULL,
    y = "Mean SSR",
    title = "Model fit by treatment"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


ggsave(
  file.path(
    out_dir,
    "Model_Performance_Summary.pdf"
  ),
  p_performance,
  width = 10,
  height = 6
)


# ============================================================
# Km SENSITIVITY ANALYSIS
# ============================================================
#
# Rationale:
# Km = 10 mg/L is retained as the primary value to reproduce
# the Patel et al. (2024) apparent process-model formulation.
# Lower Km values test whether the conclusions are sensitive to
# freshwater sediment / hyporheic values reported in the literature.
#
# Km is not chosen by whichever sensitivity run has the lowest SSR.
# The purpose of the sensitivity analysis is to test robustness of
# conclusions to the fixed Km assumption.
# ============================================================

message("\n============================================")
message("BEGINNING Km SENSITIVITY ANALYSIS")
message("============================================")


save_sample_fit_plots_primary <- save_sample_fit_plots
save_sample_fit_plots <- FALSE


sensitivity_list <- list()


for(km_test in Km_sensitivity_values) {
  
  message("Running Km = ", km_test, " mg O2/L")
  
  Km <- km_test
  
  sensitivity_list[[as.character(km_test)]] <-
    run_all_samples()
}


# Restore primary model setting
Km <- Km_primary
save_sample_fit_plots <- save_sample_fit_plots_primary


Km_sensitivity_df <- bind_rows(
  sensitivity_list
)


write.csv(
  Km_sensitivity_df,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_All_Samples.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Compare each Km against primary Km = 10 mg/L
# ------------------------------------------------------------

sensitivity_reference <- Km_sensitivity_df %>%
  filter(Km_mg_per_L == Km_primary) %>%
  select(
    Sample_Name,
    Reference_Combined_Vmax = Combined_Vmax_per_h,
    Reference_Combined_kL = Combined_kL_per_h,
    Reference_Combined_Rate = Combined_Rate_mg_per_kg_per_h,
    Reference_Best_Fit_Model = Best_Fit_Model
  )


sensitivity_comparison <- Km_sensitivity_df %>%
  left_join(
    sensitivity_reference,
    by = "Sample_Name"
  ) %>%
  mutate(
    
    Vmax_Percent_Change_From_Primary = ifelse(
      !is.na(Reference_Combined_Vmax) &
        Reference_Combined_Vmax != 0,
      100 *
        (Combined_Vmax_per_h - Reference_Combined_Vmax) /
        Reference_Combined_Vmax,
      NA_real_
    ),
    
    kL_Percent_Change_From_Primary = ifelse(
      !is.na(Reference_Combined_kL) &
        Reference_Combined_kL != 0,
      100 *
        (Combined_kL_per_h - Reference_Combined_kL) /
        Reference_Combined_kL,
      NA_real_
    ),
    
    Combined_Rate_Percent_Change_From_Primary = ifelse(
      !is.na(Reference_Combined_Rate) &
        Reference_Combined_Rate != 0,
      100 *
        (Combined_Rate_mg_per_kg_per_h - Reference_Combined_Rate) /
        Reference_Combined_Rate,
      NA_real_
    ),
    
    Best_Fit_Model_Changed =
      Best_Fit_Model != Reference_Best_Fit_Model
  )


write.csv(
  sensitivity_comparison,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Comparison_to_10.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Treatment-level sensitivity summary
# ------------------------------------------------------------

Km_sensitivity_treatment <- Km_sensitivity_df %>%
  filter(Treatment_Type != "Unknown") %>%
  group_by(
    Km_mg_per_L,
    Treatment_Type
  ) %>%
  summarise(
    n = n(),
    
    Mean_Combined_Vmax = mean(
      Combined_Vmax_per_h,
      na.rm = TRUE
    ),
    
    SD_Combined_Vmax = sd(
      Combined_Vmax_per_h,
      na.rm = TRUE
    ),
    
    Mean_Combined_kL = mean(
      Combined_kL_per_h,
      na.rm = TRUE
    ),
    
    SD_Combined_kL = sd(
      Combined_kL_per_h,
      na.rm = TRUE
    ),
    
    Mean_Combined_Rate_mg_kg_h = mean(
      Combined_Rate_mg_per_kg_per_h,
      na.rm = TRUE
    ),
    
    Mean_FirstOrder_SSR = mean(
      Linear_SSR,
      na.rm = TRUE
    ),
    
    Mean_Biotic_SSR = mean(
      Biotic_SSR,
      na.rm = TRUE
    ),
    
    Mean_Combined_SSR = mean(
      Combined_SSR,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


write.csv(
  Km_sensitivity_treatment,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Treatment_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Best-fit model across Km values
# ------------------------------------------------------------

sensitivity_model_preference <- Km_sensitivity_df %>%
  filter(Treatment_Type != "Unknown") %>%
  count(
    Km_mg_per_L,
    Best_Fit_Model
  ) %>%
  group_by(Km_mg_per_L) %>%
  mutate(
    Percentage = 100 * n / sum(n)
  ) %>%
  ungroup()


write.csv(
  sensitivity_model_preference,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Model_Preference.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Sensitivity plots
# ------------------------------------------------------------

p_vmax <- ggplot(
  Km_sensitivity_treatment,
  aes(
    x = Km_mg_per_L,
    y = Mean_Combined_Vmax,
    group = Treatment_Type,
    linetype = Treatment_Type
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(
    breaks = Km_sensitivity_values
  ) +
  labs(
    x = expression(K[m]~"(mg O"[2]~L^{-1}*")"),
    y = expression("Mean combined-model " * V[max]~"(h"^{-1}*")"),
    title = "Sensitivity of fitted Vmax to fixed Km"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Vmax.pdf"
  ),
  p_vmax,
  width = 8,
  height = 6
)


p_kl <- ggplot(
  Km_sensitivity_treatment,
  aes(
    x = Km_mg_per_L,
    y = Mean_Combined_kL,
    group = Treatment_Type,
    linetype = Treatment_Type
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(
    breaks = Km_sensitivity_values
  ) +
  labs(
    x = expression(K[m]~"(mg O"[2]~L^{-1}*")"),
    y = expression("Mean combined-model " * k[L]~"(h"^{-1}*")"),
    title = "Sensitivity of fitted kL to fixed Km"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_kL.pdf"
  ),
  p_kl,
  width = 8,
  height = 6
)


p_total <- ggplot(
  Km_sensitivity_treatment,
  aes(
    x = Km_mg_per_L,
    y = Mean_Combined_Rate_mg_kg_h,
    group = Treatment_Type,
    linetype = Treatment_Type
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(
    breaks = Km_sensitivity_values
  ) +
  labs(
    x = expression(K[m]~"(mg O"[2]~L^{-1}*")"),
    y = "Mean initial combined O2 consumption rate\n(mg O2 kg dry sediment^-1 h^-1)",
    title = "Sensitivity of total modeled O2 consumption to Km"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Total_Rate.pdf"
  ),
  p_total,
  width = 8,
  height = 6
)


# ============================================================
# Final messages
# ============================================================

message("\n============================================")
message("YEP MODELING COMPLETE")
message("Primary Km = ", Km_primary, " mg O2/L")
message("Initial normalized CB = ", CB)
message("Model fitness evaluated using SSR")
message(
  "Km sensitivity values = ",
  paste(Km_sensitivity_values, collapse = ", "),
  " mg O2/L"
)
message("Output directory: ", out_dir)
message("============================================")

print(model_performance)
print(primary_model_preference)
print(Km_sensitivity_treatment)
print(sensitivity_model_preference)
