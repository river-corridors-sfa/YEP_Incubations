# ============================================================
# Complete YEP dissolved oxygen modeling analysis
# ============================================================

rm(list = ls(all = TRUE))

library(deSolve)
library(FME)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)


# ============================================================
# Output settings
# ============================================================

out_dir <- "modeling_outputs"

if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
}

save_sample_fit_plots <- FALSE
save_development_summary_plots <- FALSE


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


# Mass/volume data
mass_volume_df <- read.csv(
  mass_volume_file_path,
  stringsAsFactors = FALSE,
  skip = 2
) %>%
  filter(grepl("YEP", Sample_Name))


# Mapping data
mapping_df <- read.csv(
  mapping_file_path,
  stringsAsFactors = FALSE
)


# Raw dissolved oxygen data
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


# Apply start/end mapping
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
  )


# Clean numeric data
df <- df %>%
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
# Fixed kinetic parameters
# ============================================================

# Primary literature-constrained oxygen half-saturation constant
Km_primary <- 0.2       # mg O2/L

# Km currently used by model
Km <- Km_primary

# Fixed dimensionless biomass normalization factor
CB <- 1.0

# Values used for sensitivity analysis
Km_sensitivity_values <- c(
  0.1,
  0.2,
  0.5,
  1.0
)

# Parameter bounds used during optimization
# These are also used to flag when the combined model collapses
# onto one of the simpler component models.
Vmax_lower_bound <- 0.01
Vmax_upper_bound <- 50
kL_lower_bound <- 0.001
kL_upper_bound <- 10

# Numerical tolerance used only to identify estimates that are
# effectively sitting on an optimization bound.
bound_tolerance <- 1e-5


# ============================================================
# Treatment classification
# ============================================================

classify_yep_treatment <- function(sample_names) {
  
  treatment_types <- case_when(
    
    grepl(
      "YEP1.*S",
      sample_names,
      ignore.case = TRUE
    ) ~ "Dry_Control",
    
    grepl(
      "YEP1.*U",
      sample_names,
      ignore.case = TRUE
    ) ~ "Dry_Unburned_DOC",
    
    grepl(
      "YEP1.*H",
      sample_names,
      ignore.case = TRUE
    ) ~ "Dry_HighBurn_DOC",
    
    grepl(
      "YEP2.*S",
      sample_names,
      ignore.case = TRUE
    ) ~ "Wet_Control",
    
    grepl(
      "YEP2.*U",
      sample_names,
      ignore.case = TRUE
    ) ~ "Wet_Unburned_DOC",
    
    grepl(
      "YEP2.*H",
      sample_names,
      ignore.case = TRUE
    ) ~ "Wet_HighBurn_DOC",
    
    TRUE ~ "Unknown"
  )
  
  return(treatment_types)
}


# ============================================================
# Sample-specific water/sediment conversion factors
# ============================================================

calculate_conversion_factors <- function(mass_volume_df) {
  
  incubation_samples <- mass_volume_df %>%
    filter(
      grepl("_INC-", Sample_Name)
    ) %>%
    mutate(
      
      Dry_Sediment_Mass_g =
        as.numeric(Dry_Sediment_Mass_g),
      
      Water_Mass_g =
        as.numeric(Water_Mass_g),
      
      # Density of water approximated as 1 g/mL
      Water_Volume_L =
        Water_Mass_g / 1000,
      
      # L water per g dry sediment
      Conversion_Factor =
        Water_Volume_L /
        Dry_Sediment_Mass_g
    ) %>%
    select(
      Sample_Name,
      Dry_Sediment_Mass_g,
      Water_Volume_L,
      Conversion_Factor
    )
  
  return(incubation_samples)
}


conversion_factors <- calculate_conversion_factors(
  mass_volume_df
)


# ============================================================
# Trim first and last 2 minutes
# ============================================================

trim_time_series <- function(data_subset) {
  
  trim_hours <- 2 / 60
  
  data_trimmed <- data_subset %>%
    filter(
      time_hr >= min(time_hr) + trim_hours,
      time_hr <= max(time_hr) - trim_hours
    )
  
  return(data_trimmed)
}


# ============================================================
# Model solver functions
#
# Adapted from Patel et al. (2024)
#
# Important difference:
# CB is held constant at 1 as a dimensionless normalization
# factor rather than modeled as a growing biomass state.
# ============================================================


# ------------------------------------------------------------
# 1. First-order O2 consumption model
#
# dDO/dt = -kL * DO
#
# kL units = h^-1
# ------------------------------------------------------------

solveLinear <- function(pars, times, DO0) {
  
  derivs <- function(t, state, pars) {
    
    with(as.list(c(state, pars)), {
      
      dDO <- -kL * DO
      
      list(c(dDO))
    })
  }
  
  
  state <- c(
    DO = DO0
  )
  
  
  tryCatch({
    
    out <- ode(
      y = state,
      times = times,
      func = derivs,
      parms = pars,
      method = "lsoda"
    )
    
    
    result <- as.data.frame(
      out[, c("time", "DO")]
    )
    
    
    if(
      any(is.na(result$DO)) ||
      any(result$DO < 0)
    ) {
      stop("Invalid solution")
    }
    
    
    return(result)
    
    
  }, error = function(e) {
    
    # Exact analytical solution for first-order decay
    result <- data.frame(
      
      time = times,
      
      DO = pmax(
        0,
        DO0 * exp(
          -pars$kL * times
        )
      )
    )
    
    return(result)
  })
}


# ------------------------------------------------------------
# 2. Biological Michaelis-Menten model
#
# dDO/dt =
# -(Vmax * DO/(DO + Km) * CB)
#
# Vmax units = mg O2/L/h
# Km units   = mg O2/L
# CB         = fixed dimensionless normalization factor
# ------------------------------------------------------------

solveBiotic <- function(
    pars,
    times,
    DO0,
    CB,
    Km
) {
  
  derivs <- function(t, state, pars) {
    
    with(as.list(c(state, pars)), {
      
      # Biological O2 consumption
      r_bio <-
        Vmax *
        DO /
        (DO + Km)
      
      
      # CB is constant and dimensionless
      dDO <-
        -r_bio *
        CB
      
      
      list(
        c(dDO)
      )
    })
  }
  
  
  # DO is the only state variable
  state <- c(
    DO = DO0
  )
  
  
  tryCatch({
    
    out <- ode(
      y = state,
      times = times,
      func = derivs,
      parms = pars,
      method = "lsoda"
    )
    
    
    result <- as.data.frame(
      out[, c("time", "DO")]
    )
    
    
    if(
      any(is.na(result$DO)) ||
      any(result$DO < 0)
    ) {
      stop("Invalid solution")
    }
    
    
    return(result)
    
    
  }, error = function(e) {
    
    # Do not replace the nonlinear equation with
    # a different approximation if the solver fails
    data.frame(
      time = times,
      DO = NA_real_
    )
  })
}


# ------------------------------------------------------------
# 3. Combined biological + first-order model
#
# dDO/dt =
# -(Vmax * DO/(DO + Km) * CB)
# - kL * DO
#
# Vmax units = mg O2/L/h
# Km units   = mg O2/L
# CB         = fixed dimensionless normalization factor
# kL units   = h^-1
# ------------------------------------------------------------

solveCombined <- function(
    pars,
    times,
    DO0,
    CB,
    Km
) {
  
  derivs <- function(t, state, pars) {
    
    with(as.list(c(state, pars)), {
      
      # Biological component
      r_bio <-
        Vmax *
        DO /
        (DO + Km)
      
      
      # First-order component
      r_first_order <-
        kL *
        DO
      
      
      # Total change in dissolved oxygen
      dDO <-
        -(r_bio * CB) -
        r_first_order
      
      
      list(
        c(dDO)
      )
    })
  }
  
  
  state <- c(
    DO = DO0
  )
  
  
  tryCatch({
    
    out <- ode(
      y = state,
      times = times,
      func = derivs,
      parms = pars,
      method = "lsoda"
    )
    
    
    result <- as.data.frame(
      out[, c("time", "DO")]
    )
    
    
    if(
      any(is.na(result$DO)) ||
      any(result$DO < 0)
    ) {
      stop("Invalid solution")
    }
    
    
    return(result)
    
    
  }, error = function(e) {
    
    data.frame(
      time = times,
      DO = NA_real_
    )
  })
}


# ============================================================
# Information criteria
# ============================================================

compute_info <- function(
    SSR,
    k,
    n
) {
  
  if(
    SSR <= 0 ||
    n <= k
  ) {
    
    return(
      list(
        AIC = Inf,
        BIC = Inf
      )
    )
  }
  
  
  AIC <-
    n *
    log(SSR / n) +
    2 * k
  
  
  BIC <-
    n *
    log(SSR / n) +
    k * log(n)
  
  
  list(
    AIC = AIC,
    BIC = BIC
  )
}


# ============================================================
# Model diagnostic helper functions
# ============================================================

calc_rmse <- function(obs, pred) {
  if(is.null(pred) || length(obs) != length(pred)) {
    return(NA_real_)
  }
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}


calc_lag1_residual_ac <- function(obs, pred) {
  if(is.null(pred) || length(obs) != length(pred)) {
    return(NA_real_)
  }
  residuals <- obs - pred
  residuals <- residuals[is.finite(residuals)]
  if(length(residuals) < 3 || sd(residuals) == 0) {
    return(NA_real_)
  }
  cor(
    residuals[-length(residuals)],
    residuals[-1],
    use = "complete.obs"
  )
}


at_lower_bound <- function(value, lower_bound, tolerance = bound_tolerance) {
  if(length(value) == 0 || is.na(value) || !is.finite(value)) {
    return(NA)
  }
  value <= (lower_bound + tolerance)
}


at_upper_bound <- function(value, upper_bound, tolerance = bound_tolerance) {
  if(length(value) == 0 || is.na(value) || !is.finite(value)) {
    return(NA)
  }
  value >= (upper_bound - tolerance)
}


# ============================================================
# Complete model fitting function
# ============================================================

fit_all_yep_models <- function(
    dat,
    sample_name,
    treatment_type = "Unknown",
    conversion_info = NULL
) {
  
  
  # ----------------------------------------------------------
  # Trim data
  # ----------------------------------------------------------
  
  dat_trimmed <- trim_time_series(
    dat
  )
  
  
  Data <- dat_trimmed %>%
    arrange(
      time_hr
    ) %>%
    select(
      time = time_hr,
      DO = DO_mg_per_L
    ) %>%
    group_by(
      time
    ) %>%
    summarise(
      DO = mean(
        DO,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    filter(
      DO > 0,
      !is.na(DO),
      !is.na(time)
    ) %>%
    arrange(
      time
    ) %>%
    as.data.frame()
  
  
  n <- nrow(Data)
  
  
  if(n < 5) {
    
    message(
      "Skipping ",
      sample_name,
      ": insufficient data after trimming (n=",
      n,
      ")"
    )
    
    return(NULL)
  }
  
  
  # Rebase time to zero
  Data$time <-
    Data$time -
    min(Data$time)
  
  
  DO0 <- Data$DO[1]
  
  
  # ----------------------------------------------------------
  # Sample-specific water and sediment quantities
  # ----------------------------------------------------------
  if(!is.null(conversion_info)) {
    
    dry_mass_g <- conversion_info$Dry_Sediment_Mass_g
    water_vol_L <- conversion_info$Water_Volume_L
    conversion_factor <- conversion_info$Conversion_Factor
    
  } else {
    
    dry_mass_g <- NA_real_
    water_vol_L <- NA_real_
    conversion_factor <- NA_real_
    
    message(
      "No sediment mass/volume normalization for ",
      sample_name
    )
  }
  
  # ----------------------------------------------------------
  # Starting parameter estimate
  # ----------------------------------------------------------
  
  linear_rate_est <- abs(
    coef(
      lm(
        DO ~ time,
        data = Data
      )
    )[["time"]]
  )
  
  
  if(
    is.na(linear_rate_est) ||
    !is.finite(linear_rate_est) ||
    linear_rate_est <= 0
  ) {
    
    linear_rate_est <- 0.1
  }
  
  
  # ==========================================================
  # Objective functions
  # ==========================================================
  
  
  # ----------------------------------------------------------
  # First-order model
  # ----------------------------------------------------------
  
  Objective_LIN <- function(x) {
    
    if(
      x[["kL"]] <= 0 ||
      x[["kL"]] > kL_upper_bound
    ) {
      
      return(1e10)
    }
    
    
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
        model = out,
        obs = Data
      ),
      
      error = function(e) 1e10
    )
  }
  
  
  # ----------------------------------------------------------
  # Biological model
  # ----------------------------------------------------------
  
  Objective_BIO <- function(x) {
    
    if(
      x[["Vmax"]] <= 0 ||
      x[["Vmax"]] > Vmax_upper_bound
    ) {
      
      return(1e10)
    }
    
    
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
        model = out,
        obs = Data
      ),
      
      error = function(e) 1e10
    )
  }
  
  
  # ----------------------------------------------------------
  # Combined model
  # ----------------------------------------------------------
  
  Objective_COM <- function(x) {
    
    if(
      x[["Vmax"]] <= 0 ||
      x[["Vmax"]] > Vmax_upper_bound ||
      x[["kL"]] <= 0 ||
      x[["kL"]] > kL_upper_bound
    ) {
      
      return(1e10)
    }
    
    
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
        model = out,
        obs = Data
      ),
      
      error = function(e) 1e10
    )
  }
  
  
  # ==========================================================
  # Fit models
  # ==========================================================
  
  fit_LIN <- NULL
  fit_BIO <- NULL
  fit_COM <- NULL
  
  
  # Starting estimates
  start_kL <-
    linear_rate_est /
    DO0
  
  start_kL <-
    max(
      kL_lower_bound,
      min(
        start_kL,
        kL_upper_bound * 0.9
      )
    )
  
  
  start_Vmax <-
    max(
      Vmax_lower_bound,
      min(
        linear_rate_est,
        Vmax_upper_bound * 0.9
      )
    )
  
  
  # First-order fit
  try({
    
    fit_LIN <- modFit(
      
      p = c(
        kL = start_kL
      ),
      
      f = Objective_LIN,
      
      lower = c(
        kL = kL_lower_bound
      ),
      
      upper = c(
        kL = kL_upper_bound
      ),
      
      control = list(
        maxiter = 500
      )
    )
    
  }, silent = TRUE)
  
  
  # Biological fit
  try({
    
    fit_BIO <- modFit(
      
      p = c(
        Vmax = start_Vmax
      ),
      
      f = Objective_BIO,
      
      lower = c(
        Vmax = Vmax_lower_bound
      ),
      
      upper = c(
        Vmax = Vmax_upper_bound
      ),
      
      control = list(
        maxiter = 500
      )
    )
    
  }, silent = TRUE)
  
  
  # Combined fit
  try({
    
    fit_COM <- modFit(
      
      p = c(
        Vmax = max(
          start_Vmax / 2,
          Vmax_lower_bound
        ),
        
        kL = max(
          start_kL / 2,
          kL_lower_bound
        )
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
      
      control = list(
        maxiter = 500
      )
    )
    
  }, silent = TRUE)
  
  
  # ==========================================================
  # Model predictions
  # ==========================================================
  
  model_LIN <- NULL
  model_BIO <- NULL
  model_COM <- NULL
  
  
  if(!is.null(fit_LIN)) {
    
    model_LIN <- solveLinear(
      
      list(
        kL =
          fit_LIN$par[["kL"]]
      ),
      
      times = Data$time,
      DO0 = DO0
    )
  }
  
  
  if(!is.null(fit_BIO)) {
    
    model_BIO <- solveBiotic(
      
      list(
        Vmax =
          fit_BIO$par[["Vmax"]]
      ),
      
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
  }
  
  
  if(!is.null(fit_COM)) {
    
    model_COM <- solveCombined(
      
      list(
        Vmax =
          fit_COM$par[["Vmax"]],
        
        kL =
          fit_COM$par[["kL"]]
      ),
      
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
  }
  
  
  # ==========================================================
  # Model diagnostics
  # ==========================================================
  #
  # AIC below is retained as a descriptive diagnostic only.
  # The 5-second observations are temporally autocorrelated, so
  # raw AIC is NOT used by itself to declare a best model.
  #
  # RMSE is reported in mg O2/L. Residual lag-1 autocorrelation
  # describes remaining temporal structure in model residuals;
  # values closer to zero indicate less serial structure.
  # ==========================================================
  
  Linear_RMSE <- if(!is.null(model_LIN)) {
    calc_rmse(Data$DO, model_LIN$DO)
  } else {
    NA_real_
  }
  
  Biotic_RMSE <- if(!is.null(model_BIO)) {
    calc_rmse(Data$DO, model_BIO$DO)
  } else {
    NA_real_
  }
  
  Combined_RMSE <- if(!is.null(model_COM)) {
    calc_rmse(Data$DO, model_COM$DO)
  } else {
    NA_real_
  }
  
  Linear_Residual_Lag1_AC <- if(!is.null(model_LIN)) {
    calc_lag1_residual_ac(Data$DO, model_LIN$DO)
  } else {
    NA_real_
  }
  
  Biotic_Residual_Lag1_AC <- if(!is.null(model_BIO)) {
    calc_lag1_residual_ac(Data$DO, model_BIO$DO)
  } else {
    NA_real_
  }
  
  Combined_Residual_Lag1_AC <- if(!is.null(model_COM)) {
    calc_lag1_residual_ac(Data$DO, model_COM$DO)
  } else {
    NA_real_
  }
  
  Combined_kL_at_lower_bound <- if(!is.null(fit_COM)) {
    at_lower_bound(
      fit_COM$par[["kL"]],
      kL_lower_bound
    )
  } else {
    NA
  }
  
  Combined_Vmax_at_lower_bound <- if(!is.null(fit_COM)) {
    at_lower_bound(
      fit_COM$par[["Vmax"]],
      Vmax_lower_bound
    )
  } else {
    NA
  }
  
  Combined_kL_at_upper_bound <- if(!is.null(fit_COM)) {
    at_upper_bound(
      fit_COM$par[["kL"]],
      kL_upper_bound
    )
  } else {
    NA
  }
  
  Combined_Vmax_at_upper_bound <- if(!is.null(fit_COM)) {
    at_upper_bound(
      fit_COM$par[["Vmax"]],
      Vmax_upper_bound
    )
  } else {
    NA
  }
  
  Best_Single_RMSE <- suppressWarnings(
    min(
      c(Linear_RMSE, Biotic_RMSE),
      na.rm = TRUE
    )
  )
  
  if(!is.finite(Best_Single_RMSE)) {
    Best_Single_RMSE <- NA_real_
  }
  
  Combined_RMSE_Improvement_vs_BestSingle_pct <-
    if(
      !is.na(Combined_RMSE) &&
      !is.na(Best_Single_RMSE) &&
      Best_Single_RMSE > 0
    ) {
      100 *
        (Best_Single_RMSE - Combined_RMSE) /
        Best_Single_RMSE
    } else {
      NA_real_
    }
  
  
  # ==========================================================
  # Optional individual fit plots
  # ==========================================================
  
  if(save_sample_fit_plots) {
    
    tryCatch({
      
      pdf_file <- file.path(
        out_dir,
        paste0(
          gsub(
            "[^A-Za-z0-9_-]",
            "_",
            sample_name
          ),
          ".pdf"
        )
      )
      
      
      pdf(
        pdf_file,
        width = 14,
        height = 10
      )
      
      
      par(
        mfrow = c(2, 2),
        mar = c(4, 4, 3, 2)
      )
      
      
      # All models
      plot(
        Data$time,
        Data$DO,
        pch = 16,
        col = "orange",
        cex = 1.3,
        xlab = "Time (h)",
        ylab = "DO (mg/L)",
        main = paste(
          sample_name,
          "\nTreatment:",
          treatment_type
        )
      )
      
      
      if(!is.null(model_LIN)) {
        lines(
          model_LIN$time,
          model_LIN$DO,
          col = "blue",
          lwd = 3
        )
      }
      
      
      if(!is.null(model_BIO)) {
        lines(
          model_BIO$time,
          model_BIO$DO,
          col = "red",
          lwd = 3,
          lty = 2
        )
      }
      
      
      if(!is.null(model_COM)) {
        lines(
          model_COM$time,
          model_COM$DO,
          col = "green",
          lwd = 3,
          lty = 3
        )
      }
      
      
      legend(
        "topright",
        c(
          "Observed",
          "First-order",
          "Biological",
          "Combined"
        ),
        pch = c(
          16,
          NA,
          NA,
          NA
        ),
        lty = c(
          NA,
          1,
          2,
          3
        ),
        col = c(
          "orange",
          "blue",
          "red",
          "green"
        ),
        lwd = c(
          NA,
          3,
          3,
          3
        ),
        bty = "n"
      )
      
      
      models_list <- list(
        
        list(
          fit = fit_LIN,
          model = model_LIN,
          title = "First-order",
          color = "blue"
        ),
        
        list(
          fit = fit_BIO,
          model = model_BIO,
          title = "Biological",
          color = "red"
        ),
        
        list(
          fit = fit_COM,
          model = model_COM,
          title = "Combined",
          color = "green"
        )
      )
      
      
      for(i in 1:3) {
        
        if(!is.null(models_list[[i]]$fit)) {
          
          plot(
            Data$time,
            Data$DO,
            pch = 16,
            col = "orange",
            xlab = "Time (h)",
            ylab = "DO (mg/L)",
            main = paste(
              models_list[[i]]$title,
              "\nSSR =",
              round(
                models_list[[i]]$fit$ssr,
                2
              )
            )
          )
          
          
          lines(
            models_list[[i]]$model$time,
            models_list[[i]]$model$DO,
            col = models_list[[i]]$color,
            lwd = 3
          )
          
          
          if(i == 1) {
            
            param_text <- paste(
              "kL =",
              round(
                fit_LIN$par[["kL"]],
                3
              ),
              "h^-1"
            )
            
          } else if(i == 2) {
            
            param_text <- paste(
              "Vmax =",
              round(
                fit_BIO$par[["Vmax"]],
                3
              ),
              "mg O2 L^-1 h^-1"
            )
            
          } else {
            
            param_text <- paste(
              "Vmax =",
              round(
                fit_COM$par[["Vmax"]],
                3
              ),
              "mg O2 L^-1 h^-1\n",
              "kL =",
              round(
                fit_COM$par[["kL"]],
                3
              ),
              "h^-1"
            )
          }
          
          
          text(
            x = min(Data$time),
            y = max(Data$DO) * 0.85,
            labels = param_text,
            adj = c(0, 1)
          )
          
        } else {
          
          plot(
            1,
            1,
            type = "n",
            xlab = "",
            ylab = "",
            main = "Fit failed"
          )
        }
      }
      
      
      dev.off()
      
      
    }, error = function(e) {
      
      message(
        "Error creating plot for ",
        sample_name,
        ": ",
        e$message
      )
    })
  }
  
  
  # ==========================================================
  # Initial oxygen-consumption rates
  # ==========================================================
  
  
  # ----------------------------------------------------------
  # First-order model
  # mg O2/L/h
  # ----------------------------------------------------------
  
  DO_rate_linear <- if(!is.null(fit_LIN)) {
    
    fit_LIN$par[["kL"]] *
      DO0
    
  } else {
    
    NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Biological model
  # mg O2/L/h
  # ----------------------------------------------------------
  
  DO_rate_biotic <- if(!is.null(fit_BIO)) {
    
    fit_BIO$par[["Vmax"]] *
      (
        DO0 /
          (DO0 + Km)
      ) *
      CB
    
  } else {
    
    NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Combined model components
  # mg O2/L/h
  # ----------------------------------------------------------
  
  bio_rate_combined <- if(!is.null(fit_COM)) {
    
    fit_COM$par[["Vmax"]] *
      (
        DO0 /
          (DO0 + Km)
      ) *
      CB
    
  } else {
    
    NA_real_
  }
  
  
  first_order_rate_combined <- if(!is.null(fit_COM)) {
    
    fit_COM$par[["kL"]] *
      DO0
    
  } else {
    
    NA_real_
  }
  
  
  DO_rate_combined <- if(!is.null(fit_COM)) {
    
    bio_rate_combined +
      first_order_rate_combined
    
  } else {
    
    NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Relative model-attributed contributions
  # at initial DO
  # ----------------------------------------------------------
  
  Combined_Biotic_Fraction_Initial <-
    if(
      !is.null(fit_COM) &&
      !is.na(DO_rate_combined) &&
      DO_rate_combined > 0
    ) {
      
      bio_rate_combined /
        DO_rate_combined
      
    } else {
      
      NA_real_
    }
  
  
  Combined_FirstOrder_Fraction_Initial <-
    if(
      !is.null(fit_COM) &&
      !is.na(DO_rate_combined) &&
      DO_rate_combined > 0
    ) {
      
      first_order_rate_combined /
        DO_rate_combined
      
    } else {
      
      NA_real_
    }
  
  
  # ==========================================================
  # Mass-normalize instantaneous rates
  #
  # mg/L/h
  # x L/g dry sediment
  # = mg/g dry sediment/h
  # ==========================================================
  
  DO_rate_linear_per_g <-
    DO_rate_linear *
    conversion_factor
  
  
  DO_rate_biotic_per_g <-
    DO_rate_biotic *
    conversion_factor
  
  
  DO_rate_combined_per_g <-
    DO_rate_combined *
    conversion_factor
  
  
  bio_rate_combined_per_g <-
    bio_rate_combined *
    conversion_factor
  
  
  first_order_rate_combined_per_g <-
    first_order_rate_combined *
    conversion_factor
  
  
  # Convert g sediment denominator to kg sediment denominator
  
  DO_rate_linear_per_kg <-
    DO_rate_linear_per_g *
    1000
  
  
  DO_rate_biotic_per_kg <-
    DO_rate_biotic_per_g *
    1000
  
  
  DO_rate_combined_per_kg <-
    DO_rate_combined_per_g *
    1000
  
  
  bio_rate_combined_per_kg <-
    bio_rate_combined_per_g *
    1000
  
  
  first_order_rate_combined_per_kg <-
    first_order_rate_combined_per_g *
    1000
  
  
  # ==========================================================
  # Mass-normalize Vmax
  #
  # Fitted Vmax remains mg O2/L/h in the ODE.
  #
  # These additional values are for comparison among
  # sediment samples with different water:sediment ratios.
  # ==========================================================
  
  Biotic_Vmax_mg_per_g_per_h <-
    if(!is.null(fit_BIO)) {
      
      fit_BIO$par[["Vmax"]] *
        conversion_factor
      
    } else {
      
      NA_real_
    }
  
  
  Biotic_Vmax_mg_per_kg_per_h <-
    Biotic_Vmax_mg_per_g_per_h *
    1000
  
  
  Combined_Vmax_mg_per_g_per_h <-
    if(!is.null(fit_COM)) {
      
      fit_COM$par[["Vmax"]] *
        conversion_factor
      
    } else {
      
      NA_real_
    }
  
  
  Combined_Vmax_mg_per_kg_per_h <-
    Combined_Vmax_mg_per_g_per_h *
    1000
  
  
  # ==========================================================
  # Summary row
  # ==========================================================
  
  summary_row <- data.frame(
    
    Sample_Name =
      sample_name,
    
    Treatment_Type =
      treatment_type,
    
    n_points =
      n,
    
    Initial_DO_mg_per_L =
      round(
        DO0,
        4
      ),
    
    Km_mg_per_L =
      Km,
    
    CB_dimensionless =
      CB,
    
    Dry_Sediment_Mass_g =
      round(
        dry_mass_g,
        4
      ),
    
    Water_Volume_mL =
      round(
        water_vol_L * 1000,
        4
      ),
    
    Conversion_Factor_L_per_g =
      round(
        conversion_factor,
        6
      ),
    
    
    # --------------------------------------------------------
    # First-order model
    # --------------------------------------------------------
    
    Linear_kL_per_h =
      ifelse(
        !is.null(fit_LIN),
        round(
          fit_LIN$par[["kL"]],
          6
        ),
        NA
      ),
    
    Linear_Rate_mg_per_L_per_h =
      round(
        DO_rate_linear,
        6
      ),
    
    Linear_Rate_mg_per_g_per_h =
      round(
        DO_rate_linear_per_g,
        6
      ),
    
    Linear_Rate_mg_per_kg_per_h =
      round(
        DO_rate_linear_per_kg,
        4
      ),
    
    Linear_SSR =
      ifelse(
        !is.null(fit_LIN),
        round(
          fit_LIN$ssr,
          4
        ),
        NA
      ),
    
    Linear_AIC =
      ifelse(
        !is.null(fit_LIN),
        round(
          compute_info(
            fit_LIN$ssr,
            1,
            n
          )$AIC,
          4
        ),
        NA
      ),
    
    Linear_RMSE_mg_per_L =
      round(Linear_RMSE, 6),
    
    Linear_Residual_Lag1_AC =
      round(Linear_Residual_Lag1_AC, 6),
    
    
    # --------------------------------------------------------
    # Biological model
    # --------------------------------------------------------
    
    Biotic_Vmax_mg_per_L_per_h =
      ifelse(
        !is.null(fit_BIO),
        round(
          fit_BIO$par[["Vmax"]],
          6
        ),
        NA
      ),
    
    Biotic_Vmax_mg_per_g_per_h =
      round(
        Biotic_Vmax_mg_per_g_per_h,
        6
      ),
    
    Biotic_Vmax_mg_per_kg_per_h =
      round(
        Biotic_Vmax_mg_per_kg_per_h,
        4
      ),
    
    Biotic_Rate_mg_per_L_per_h =
      round(
        DO_rate_biotic,
        6
      ),
    
    Biotic_Rate_mg_per_g_per_h =
      round(
        DO_rate_biotic_per_g,
        6
      ),
    
    Biotic_Rate_mg_per_kg_per_h =
      round(
        DO_rate_biotic_per_kg,
        4
      ),
    
    Biotic_SSR =
      ifelse(
        !is.null(fit_BIO),
        round(
          fit_BIO$ssr,
          4
        ),
        NA
      ),
    
    Biotic_AIC =
      ifelse(
        !is.null(fit_BIO),
        round(
          compute_info(
            fit_BIO$ssr,
            1,
            n
          )$AIC,
          4
        ),
        NA
      ),
    
    Biotic_RMSE_mg_per_L =
      round(Biotic_RMSE, 6),
    
    Biotic_Residual_Lag1_AC =
      round(Biotic_Residual_Lag1_AC, 6),
    
    
    # --------------------------------------------------------
    # Combined model
    # --------------------------------------------------------
    
    Combined_Vmax_mg_per_L_per_h =
      ifelse(
        !is.null(fit_COM),
        round(
          fit_COM$par[["Vmax"]],
          6
        ),
        NA
      ),
    
    Combined_Vmax_mg_per_g_per_h =
      round(
        Combined_Vmax_mg_per_g_per_h,
        6
      ),
    
    Combined_Vmax_mg_per_kg_per_h =
      round(
        Combined_Vmax_mg_per_kg_per_h,
        4
      ),
    
    Combined_kL_per_h =
      ifelse(
        !is.null(fit_COM),
        round(
          fit_COM$par[["kL"]],
          6
        ),
        NA
      ),
    
    Combined_Biotic_Rate_mg_per_L_per_h =
      round(
        bio_rate_combined,
        6
      ),
    
    Combined_FirstOrder_Rate_mg_per_L_per_h =
      round(
        first_order_rate_combined,
        6
      ),
    
    Combined_Rate_mg_per_L_per_h =
      round(
        DO_rate_combined,
        6
      ),
    
    Combined_Biotic_Rate_mg_per_kg_per_h =
      round(
        bio_rate_combined_per_kg,
        4
      ),
    
    Combined_FirstOrder_Rate_mg_per_kg_per_h =
      round(
        first_order_rate_combined_per_kg,
        4
      ),
    
    Combined_Rate_mg_per_g_per_h =
      round(
        DO_rate_combined_per_g,
        6
      ),
    
    Combined_Rate_mg_per_kg_per_h =
      round(
        DO_rate_combined_per_kg,
        4
      ),
    
    Combined_Biotic_Fraction_Initial =
      round(
        Combined_Biotic_Fraction_Initial,
        4
      ),
    
    Combined_FirstOrder_Fraction_Initial =
      round(
        Combined_FirstOrder_Fraction_Initial,
        4
      ),
    
    Combined_SSR =
      ifelse(
        !is.null(fit_COM),
        round(
          fit_COM$ssr,
          4
        ),
        NA
      ),
    
    Combined_AIC =
      ifelse(
        !is.null(fit_COM),
        round(
          compute_info(
            fit_COM$ssr,
            2,
            n
          )$AIC,
          4
        ),
        NA
      ),
    
    Combined_RMSE_mg_per_L =
      round(Combined_RMSE, 6),
    
    Combined_Residual_Lag1_AC =
      round(Combined_Residual_Lag1_AC, 6),
    
    Combined_RMSE_Improvement_vs_BestSingle_pct =
      round(
        Combined_RMSE_Improvement_vs_BestSingle_pct,
        4
      ),
    
    Combined_kL_at_lower_bound =
      Combined_kL_at_lower_bound,
    
    Combined_Vmax_at_lower_bound =
      Combined_Vmax_at_lower_bound,
    
    Combined_kL_at_upper_bound =
      Combined_kL_at_upper_bound,
    
    Combined_Vmax_at_upper_bound =
      Combined_Vmax_at_upper_bound,
    
    stringsAsFactors = FALSE
  )
  
  
  # ==========================================================
  # Model interpretation
  # ==========================================================
  #
  # Raw AIC is retained as a diagnostic, but because the DO
  # measurements occur every 5 seconds and are temporally
  # autocorrelated, it is not treated as the sole model-selection
  # criterion.
  #
  # Interpretation hierarchy:
  # 1. If one combined-model parameter is effectively at its lower
  #    bound, the combined model has collapsed onto the simpler model.
  # 2. If both combined parameters are away from their bounds,
  #    support for the combined model requires BOTH lower RMSE and
  #    lower absolute lag-1 residual autocorrelation than each
  #    single-component model.
  # 3. Conflicting diagnostics are flagged for manual review.
  # ==========================================================
  
  aic_values <- c(
    summary_row$Linear_AIC,
    summary_row$Biotic_AIC,
    summary_row$Combined_AIC
  )
  
  aic_names <- c(
    "FirstOrder",
    "Biotic",
    "Combined"
  )
  
  valid_aic <-
    !is.na(aic_values) &
    is.finite(aic_values)
  
  if(any(valid_aic)) {
    lowest_aic_idx <- which.min(aic_values[valid_aic])
    summary_row$Lowest_AIC_Model <-
      aic_names[valid_aic][lowest_aic_idx]
    summary_row$Lowest_AIC <-
      min(aic_values[valid_aic])
  } else {
    summary_row$Lowest_AIC_Model <- "None"
    summary_row$Lowest_AIC <- NA_real_
  }
  
  
  # Best single-component model by RMSE
  if(
    !is.na(Linear_RMSE) &&
    !is.na(Biotic_RMSE)
  ) {
    Best_Single_Model_RMSE <-
      ifelse(
        Linear_RMSE <= Biotic_RMSE,
        "FirstOrder",
        "Biotic"
      )
  } else if(!is.na(Linear_RMSE)) {
    Best_Single_Model_RMSE <- "FirstOrder"
  } else if(!is.na(Biotic_RMSE)) {
    Best_Single_Model_RMSE <- "Biotic"
  } else {
    Best_Single_Model_RMSE <- "None"
  }
  
  summary_row$Best_Single_Model_RMSE <-
    Best_Single_Model_RMSE
  
  
  combined_lower_rmse_than_both <-
    !is.na(Combined_RMSE) &&
    !is.na(Linear_RMSE) &&
    !is.na(Biotic_RMSE) &&
    Combined_RMSE < Linear_RMSE &&
    Combined_RMSE < Biotic_RMSE
  
  combined_lower_abs_ac_than_both <-
    !is.na(Combined_Residual_Lag1_AC) &&
    !is.na(Linear_Residual_Lag1_AC) &&
    !is.na(Biotic_Residual_Lag1_AC) &&
    abs(Combined_Residual_Lag1_AC) <
    abs(Linear_Residual_Lag1_AC) &&
    abs(Combined_Residual_Lag1_AC) <
    abs(Biotic_Residual_Lag1_AC)
  
  summary_row$Combined_Lower_RMSE_Than_Both_Singles <-
    combined_lower_rmse_than_both
  
  summary_row$Combined_Lower_Abs_Lag1_AC_Than_Both_Singles <-
    combined_lower_abs_ac_than_both
  
  
  if(is.null(fit_COM)) {
    
    summary_row$Model_Interpretation <-
      Best_Single_Model_RMSE
    
    summary_row$Model_Interpretation_Basis <-
      "Combined fit failed; selected lower-RMSE single-component model"
    
  } else if(
    isTRUE(Combined_kL_at_lower_bound) &&
    !isTRUE(Combined_Vmax_at_lower_bound)
  ) {
    
    summary_row$Model_Interpretation <- "Biotic"
    
    summary_row$Model_Interpretation_Basis <-
      "Combined kL is at its lower bound; combined model collapses to biotic model"
    
  } else if(
    isTRUE(Combined_Vmax_at_lower_bound) &&
    !isTRUE(Combined_kL_at_lower_bound)
  ) {
    
    summary_row$Model_Interpretation <- "FirstOrder"
    
    summary_row$Model_Interpretation_Basis <-
      "Combined Vmax is at its lower bound; combined model collapses to first-order model"
    
  } else if(
    isTRUE(Combined_kL_at_lower_bound) &&
    isTRUE(Combined_Vmax_at_lower_bound)
  ) {
    
    summary_row$Model_Interpretation <- "ManualReview"
    
    summary_row$Model_Interpretation_Basis <-
      "Both combined-model parameters are at lower bounds"
    
  } else if(
    isTRUE(combined_lower_rmse_than_both) &&
    isTRUE(combined_lower_abs_ac_than_both)
  ) {
    
    summary_row$Model_Interpretation <- "Combined"
    
    summary_row$Model_Interpretation_Basis <-
      "Both parameters are interior; combined model has lower RMSE and lower absolute lag-1 residual autocorrelation than both single models"
    
  } else {
    
    summary_row$Model_Interpretation <- "ManualReview"
    
    summary_row$Model_Interpretation_Basis <-
      "Combined parameters are interior but RMSE and residual-autocorrelation diagnostics do not both favor the combined model"
  }
  
  
  return(
    
    list(
      
      summary =
        summary_row,
      
      fits = list(
        first_order = fit_LIN,
        biotic = fit_BIO,
        combined = fit_COM
      ),
      
      data =
        Data
    )
  )
}


# ============================================================
# Evaluation plots
# ============================================================

create_evaluation_plots <- function(summary_df) {
  
  plot_data <- summary_df %>%
    filter(
      Treatment_Type != "Unknown"
    )
  
  
  # ----------------------------------------------------------
  # Model performance
  # ----------------------------------------------------------
  
  pdf(
    file.path(
      out_dir,
      "Model_Performance_Summary.pdf"
    ),
    width = 12,
    height = 8
  )
  
  
  rmse_summary <- plot_data %>%
    group_by(
      Treatment_Type
    ) %>%
    summarise(
      
      FirstOrder_RMSE =
        mean(
          Linear_RMSE_mg_per_L,
          na.rm = TRUE
        ),
      
      Biotic_RMSE =
        mean(
          Biotic_RMSE_mg_per_L,
          na.rm = TRUE
        ),
      
      Combined_RMSE =
        mean(
          Combined_RMSE_mg_per_L,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  rmse_matrix <- as.matrix(
    rmse_summary[
      ,
      c(
        "FirstOrder_RMSE",
        "Biotic_RMSE",
        "Combined_RMSE"
      )
    ]
  )
  
  
  rownames(rmse_matrix) <-
    rmse_summary$Treatment_Type
  
  
  barplot(
    t(rmse_matrix),
    beside = TRUE,
    col = c(
      "lightblue",
      "lightcoral",
      "lightgreen"
    ),
    main =
      "Model Performance by Treatment\nLower RMSE = Better Fit",
    ylab =
      "RMSE (mg O2/L)",
    legend.text = c(
      "First-order",
      "Biotic",
      "Combined"
    ),
    args.legend = list(
      x = "topright",
      bty = "n"
    ),
    las = 2
  )
  
  
  dev.off()
  
  
  # ----------------------------------------------------------
  # Optional development plots
  # ----------------------------------------------------------
  
  if(save_development_summary_plots) {
    
    
    # --------------------------------------------------------
    # Combined-model fitted parameters
    # --------------------------------------------------------
    
    parameter_summary <- plot_data %>%
      group_by(
        Treatment_Type
      ) %>%
      summarise(
        
        n = n(),
        
        Mean_Vmax =
          mean(
            Combined_Vmax_mg_per_kg_per_h,
            na.rm = TRUE
          ),
        
        SE_Vmax =
          sd(
            Combined_Vmax_mg_per_kg_per_h,
            na.rm = TRUE
          ) /
          sqrt(n()),
        
        Mean_kL =
          mean(
            Combined_kL_per_h,
            na.rm = TRUE
          ),
        
        SE_kL =
          sd(
            Combined_kL_per_h,
            na.rm = TRUE
          ) /
          sqrt(n()),
        
        .groups = "drop"
      )
    
    
    pdf(
      file.path(
        out_dir,
        "Combined_Model_Parameters.pdf"
      ),
      width = 12,
      height = 6
    )
    
    
    par(
      mfrow = c(1, 2)
    )
    
    
    bp1 <- barplot(
      parameter_summary$Mean_Vmax,
      names.arg =
        gsub(
          "_",
          "\n",
          parameter_summary$Treatment_Type
        ),
      main =
        "Combined-model Vmax",
      ylab =
        "Vmax (mg O2/kg dry sediment/h)",
      las = 2
    )
    
    
    arrows(
      bp1,
      parameter_summary$Mean_Vmax -
        parameter_summary$SE_Vmax,
      bp1,
      parameter_summary$Mean_Vmax +
        parameter_summary$SE_Vmax,
      angle = 90,
      code = 3,
      length = 0.05
    )
    
    
    bp2 <- barplot(
      parameter_summary$Mean_kL,
      names.arg =
        gsub(
          "_",
          "\n",
          parameter_summary$Treatment_Type
        ),
      main =
        "Combined-model kL",
      ylab =
        "kL (h^-1)",
      las = 2
    )
    
    
    arrows(
      bp2,
      parameter_summary$Mean_kL -
        parameter_summary$SE_kL,
      bp2,
      parameter_summary$Mean_kL +
        parameter_summary$SE_kL,
      angle = 90,
      code = 3,
      length = 0.05
    )
    
    
    dev.off()
    
    
    # --------------------------------------------------------
    # Mass-normalized initial rates
    # --------------------------------------------------------
    
    mass_rate_summary <- plot_data %>%
      group_by(
        Treatment_Type
      ) %>%
      summarise(
        
        n = n(),
        
        FirstOrder_Rate =
          mean(
            Linear_Rate_mg_per_kg_per_h,
            na.rm = TRUE
          ),
        
        Biotic_Rate =
          mean(
            Biotic_Rate_mg_per_kg_per_h,
            na.rm = TRUE
          ),
        
        Combined_Rate =
          mean(
            Combined_Rate_mg_per_kg_per_h,
            na.rm = TRUE
          ),
        
        .groups = "drop"
      )
    
    
    mass_matrix <- as.matrix(
      mass_rate_summary[
        ,
        c(
          "FirstOrder_Rate",
          "Biotic_Rate",
          "Combined_Rate"
        )
      ]
    )
    
    
    rownames(mass_matrix) <-
      mass_rate_summary$Treatment_Type
    
    
    pdf(
      file.path(
        out_dir,
        "Mass_Normalized_Initial_Rates.pdf"
      ),
      width = 12,
      height = 7
    )
    
    
    barplot(
      t(mass_matrix),
      beside = TRUE,
      col = c(
        "lightblue",
        "lightcoral",
        "lightgreen"
      ),
      ylab =
        "Initial O2 consumption rate (mg O2/kg dry sediment/h)",
      las = 2,
      legend.text = c(
        "First-order",
        "Biotic",
        "Combined"
      ),
      args.legend = list(
        x = "topright",
        bty = "n"
      )
    )
    
    
    dev.off()
  }
}


# ============================================================
# Model-comparison summary table
# ============================================================

create_summary_tables <- function(summary_df) {
  
  table2 <- summary_df %>%
    filter(
      Treatment_Type != "Unknown"
    ) %>%
    group_by(
      Treatment_Type
    ) %>%
    summarise(
      
      n = n(),
      
      FirstOrder_SSR =
        round(mean(Linear_SSR, na.rm = TRUE), 2),
      
      Biotic_SSR =
        round(mean(Biotic_SSR, na.rm = TRUE), 2),
      
      Combined_SSR =
        round(mean(Combined_SSR, na.rm = TRUE), 2),
      
      FirstOrder_RMSE =
        round(mean(Linear_RMSE_mg_per_L, na.rm = TRUE), 4),
      
      Biotic_RMSE =
        round(mean(Biotic_RMSE_mg_per_L, na.rm = TRUE), 4),
      
      Combined_RMSE =
        round(mean(Combined_RMSE_mg_per_L, na.rm = TRUE), 4),
      
      FirstOrder_Mean_Abs_Lag1_AC =
        round(mean(abs(Linear_Residual_Lag1_AC), na.rm = TRUE), 4),
      
      Biotic_Mean_Abs_Lag1_AC =
        round(mean(abs(Biotic_Residual_Lag1_AC), na.rm = TRUE), 4),
      
      Combined_Mean_Abs_Lag1_AC =
        round(mean(abs(Combined_Residual_Lag1_AC), na.rm = TRUE), 4),
      
      Combined_kL_at_lower_bound_n =
        sum(Combined_kL_at_lower_bound %in% TRUE),
      
      Combined_Vmax_at_lower_bound_n =
        sum(Combined_Vmax_at_lower_bound %in% TRUE),
      
      Most_Common_Model_Interpretation = {
        x <- Model_Interpretation[!is.na(Model_Interpretation)]
        if(length(x) == 0) {
          NA_character_
        } else {
          names(which.max(table(x)))
        }
      },
      
      .groups = "drop"
    )
  
  
  write.csv(
    table2,
    file.path(
      out_dir,
      "YEP_Model_Performance_Table.csv"
    ),
    row.names = FALSE
  )
  
  return(table2)
}


# ============================================================
# Prepare data for model fitting
# ============================================================

df2 <- df %>%
  group_by(
    Sample_Name
  ) %>%
  arrange(
    Elapsed_Seconds,
    .by_group = TRUE
  ) %>%
  mutate(
    
    time_hr =
      (
        Elapsed_Seconds -
          min(Elapsed_Seconds)
      ) /
      3600
  ) %>%
  ungroup() %>%
  mutate(
    
    Treatment_Type =
      classify_yep_treatment(
        Sample_Name
      )
  )


# Treatment check
treatment_summary <- df2 %>%
  count(
    Treatment_Type,
    sort = TRUE
  )

print(
  treatment_summary
)


# Conversion-factor check
print(
  head(
    conversion_factors
  )
)


samples <- unique(
  df2$Sample_Name
)


# ============================================================
# PRIMARY ANALYSIS
#
# Km = 0.2 mg O2/L
# CB = 1 dimensionless and constant
# ============================================================

Km <- Km_primary


message(
  "\n============================================"
)

message(
  "PRIMARY YEP ANALYSIS"
)

message(
  "Km = ",
  Km,
  " mg O2/L"
)

message(
  "CB = ",
  CB,
  " dimensionless and constant"
)

message(
  "============================================"
)


all_results <- lapply(
  samples,
  function(sname) {
    
    dat_s <- df2 %>%
      filter(
        Sample_Name == sname
      )
    
    
    treatment_type <-
      unique(
        dat_s$Treatment_Type
      )[1]
    
    
    conversion_info <-
      conversion_factors %>%
      filter(
        Sample_Name == sname
      )
    
    
    if(nrow(conversion_info) == 0) {
      
      conversion_info <- NULL
      
      message(
        "No sediment mass/volume normalization for ",
        sname,
        "; mass-normalized outputs will be NA."
      )
      
    } else {
      
      conversion_info <-
        conversion_info[1, ]
    }
    
    
    fit_all_yep_models(
      dat = dat_s,
      sample_name = sname,
      treatment_type = treatment_type,
      conversion_info = conversion_info
    )
  }
)


# Combine primary results
summary_df <- bind_rows(
  lapply(
    all_results,
    function(x) {
      
      if(!is.null(x)) {
        x$summary
      } else {
        NULL
      }
    }
  )
)


# ============================================================
# Save primary analysis
# ============================================================

if(nrow(summary_df) > 0) {
  
  
  write.csv(
    summary_df,
    file.path(
      out_dir,
      "YEP_Complete_Modeling_Outputs.csv"
    ),
    row.names = FALSE
  )
  
  
  create_evaluation_plots(
    summary_df
  )
  
  
  model_table <- create_summary_tables(
    summary_df
  )
  
  
  # ----------------------------------------------------------
  # Treatment-level summary
  # ----------------------------------------------------------
  
  treatment_results <- summary_df %>%
    filter(
      Treatment_Type != "Unknown"
    ) %>%
    group_by(
      Treatment_Type
    ) %>%
    summarise(
      
      n = n(),
      
      Mean_Combined_Vmax_mg_per_kg_per_h =
        mean(
          Combined_Vmax_mg_per_kg_per_h,
          na.rm = TRUE
        ),
      
      SD_Combined_Vmax_mg_per_kg_per_h =
        sd(
          Combined_Vmax_mg_per_kg_per_h,
          na.rm = TRUE
        ),
      
      Mean_Combined_kL_per_h =
        mean(
          Combined_kL_per_h,
          na.rm = TRUE
        ),
      
      SD_Combined_kL_per_h =
        sd(
          Combined_kL_per_h,
          na.rm = TRUE
        ),
      
      Mean_Combined_Rate_mg_per_kg_per_h =
        mean(
          Combined_Rate_mg_per_kg_per_h,
          na.rm = TRUE
        ),
      
      Mean_Biotic_Fraction_Initial =
        mean(
          Combined_Biotic_Fraction_Initial,
          na.rm = TRUE
        ),
      
      Mean_FirstOrder_Fraction_Initial =
        mean(
          Combined_FirstOrder_Fraction_Initial,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  write.csv(
    treatment_results,
    file.path(
      out_dir,
      "YEP_Treatment_Model_Summary.csv"
    ),
    row.names = FALSE
  )
  
  
  print(
    treatment_results
  )
  
  
  # ----------------------------------------------------------
  # Overall model interpretation
  # ----------------------------------------------------------
  
  model_preference <- summary_df %>%
    count(
      Model_Interpretation,
      sort = TRUE
    ) %>%
    mutate(
      
      Percentage =
        round(
          n /
            sum(n) *
            100,
          1
        )
    )
  
  
  print(
    model_preference
  )
  
  
} else {
  
  message(
    "No successful primary analyses completed."
  )
}


# ============================================================
# Km SENSITIVITY ANALYSIS
#
# 0.1 mg/L = lower bracketing scenario
# 0.2 mg/L = primary literature-constrained value
# 0.5 mg/L = moderate scenario
# 1.0 mg/L = upper literature-informed scenario
#
# Purpose:
# Determine whether fitted Vmax, kL, model interpretation,
# and treatment patterns depend strongly on fixed Km.
#
# Km is NOT selected by choosing the value with lowest AIC.
# ============================================================


message(
  "\n============================================"
)

message(
  "BEGINNING Km SENSITIVITY ANALYSIS"
)

message(
  "============================================"
)


# Save plotting settings
save_sample_fit_plots_original <-
  save_sample_fit_plots

save_development_summary_plots_original <-
  save_development_summary_plots


# Do not make sample PDFs repeatedly
save_sample_fit_plots <- FALSE
save_development_summary_plots <- FALSE


sensitivity_results_list <- list()


for(km_test in Km_sensitivity_values) {
  
  
  message(
    "\n--------------------------------------------"
  )
  
  message(
    "Km = ",
    km_test,
    " mg O2/L"
  )
  
  message(
    "--------------------------------------------"
  )
  
  
  # Set Km for this sensitivity run
  Km <- km_test
  
  
  sensitivity_results_km <- lapply(
    samples,
    function(sname) {
      
      dat_s <- df2 %>%
        filter(
          Sample_Name == sname
        )
      
      
      treatment_type <-
        unique(
          dat_s$Treatment_Type
        )[1]
      
      
      conversion_info <-
        conversion_factors %>%
        filter(
          Sample_Name == sname
        )
      
      
      if(nrow(conversion_info) == 0) {
        
        conversion_info <- NULL
        
      } else {
        
        conversion_info <-
          conversion_info[1, ]
      }
      
      
      fit_all_yep_models(
        dat = dat_s,
        sample_name = sname,
        treatment_type = treatment_type,
        conversion_info = conversion_info
      )
    }
  )
  
  
  sensitivity_summary_km <-
    bind_rows(
      
      lapply(
        sensitivity_results_km,
        function(x) {
          
          if(!is.null(x)) {
            x$summary
          } else {
            NULL
          }
        }
      )
    )
  
  
  if(nrow(sensitivity_summary_km) > 0) {
    
    sensitivity_results_list[[as.character(km_test)]] <- sensitivity_summary_km
    
  }
}


# Restore primary Km
Km <- Km_primary


# Restore plotting options
save_sample_fit_plots <-
  save_sample_fit_plots_original

save_development_summary_plots <-
  save_development_summary_plots_original


# ============================================================
# Combine sensitivity results
# ============================================================

sensitivity_df <- bind_rows(
  sensitivity_results_list
)


write.csv(
  sensitivity_df,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_All_Samples.csv"
  ),
  row.names = FALSE
)


# ============================================================
# Compare each Km against primary Km = 0.2
# ============================================================

sensitivity_reference <- sensitivity_df %>%
  filter(
    Km_mg_per_L == Km_primary
  ) %>%
  select(
    
    Sample_Name,
    
    Reference_Combined_Vmax_mg_per_kg_per_h =
      Combined_Vmax_mg_per_kg_per_h,
    
    Reference_Combined_kL_per_h =
      Combined_kL_per_h,
    
    Reference_Combined_Rate_mg_per_kg_per_h =
      Combined_Rate_mg_per_kg_per_h,
    
    Reference_Model_Interpretation =
      Model_Interpretation
  )


sensitivity_comparison <- sensitivity_df %>%
  left_join(
    sensitivity_reference,
    by = "Sample_Name"
  ) %>%
  mutate(
    
    Vmax_Percent_Change_From_Primary =
      ifelse(
        
        !is.na(
          Reference_Combined_Vmax_mg_per_kg_per_h
        ) &
          Reference_Combined_Vmax_mg_per_kg_per_h != 0,
        
        100 *
          (
            Combined_Vmax_mg_per_kg_per_h -
              Reference_Combined_Vmax_mg_per_kg_per_h
          ) /
          Reference_Combined_Vmax_mg_per_kg_per_h,
        
        NA
      ),
    
    
    kL_Percent_Change_From_Primary =
      ifelse(
        
        !is.na(
          Reference_Combined_kL_per_h
        ) &
          Reference_Combined_kL_per_h != 0,
        
        100 *
          (
            Combined_kL_per_h -
              Reference_Combined_kL_per_h
          ) /
          Reference_Combined_kL_per_h,
        
        NA
      ),
    
    
    Combined_Rate_Percent_Change_From_Primary =
      ifelse(
        
        !is.na(
          Reference_Combined_Rate_mg_per_kg_per_h
        ) &
          Reference_Combined_Rate_mg_per_kg_per_h != 0,
        
        100 *
          (
            Combined_Rate_mg_per_kg_per_h -
              Reference_Combined_Rate_mg_per_kg_per_h
          ) /
          Reference_Combined_Rate_mg_per_kg_per_h,
        
        NA
      ),
    
    
    Model_Interpretation_Changed =
      Model_Interpretation !=
      Reference_Model_Interpretation
  )


write.csv(
  sensitivity_comparison,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Comparison_to_0.2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# Treatment-level sensitivity summary
# ============================================================

sensitivity_treatment_summary <- sensitivity_df %>%
  filter(
    Treatment_Type != "Unknown"
  ) %>%
  group_by(
    Km_mg_per_L,
    Treatment_Type
  ) %>%
  summarise(
    
    n = n(),
    
    Mean_Combined_Vmax_mg_per_kg_per_h =
      mean(
        Combined_Vmax_mg_per_kg_per_h,
        na.rm = TRUE
      ),
    
    SD_Combined_Vmax_mg_per_kg_per_h =
      sd(
        Combined_Vmax_mg_per_kg_per_h,
        na.rm = TRUE
      ),
    
    Mean_Combined_kL_per_h =
      mean(
        Combined_kL_per_h,
        na.rm = TRUE
      ),
    
    SD_Combined_kL_per_h =
      sd(
        Combined_kL_per_h,
        na.rm = TRUE
      ),
    
    Mean_Combined_Rate_mg_per_kg_per_h =
      mean(
        Combined_Rate_mg_per_kg_per_h,
        na.rm = TRUE
      ),
    
    Mean_Biotic_Fraction_Initial =
      mean(
        Combined_Biotic_Fraction_Initial,
        na.rm = TRUE
      ),
    
    Mean_FirstOrder_Fraction_Initial =
      mean(
        Combined_FirstOrder_Fraction_Initial,
        na.rm = TRUE
      ),
    
    Mean_Combined_AIC =
      mean(
        Combined_AIC,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


write.csv(
  sensitivity_treatment_summary,
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Treatment_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# Model interpretation across Km values
# ============================================================

sensitivity_model_preference <- sensitivity_df %>%
  count(
    Km_mg_per_L,
    Model_Interpretation
  ) %>%
  group_by(
    Km_mg_per_L
  ) %>%
  mutate(
    
    Percentage =
      100 *
      n /
      sum(n)
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


# ============================================================
# Sensitivity plots
# ============================================================


# ------------------------------------------------------------
# Vmax sensitivity
# ------------------------------------------------------------

p_sens_vmax <- sensitivity_df %>%
  filter(
    Treatment_Type != "Unknown"
  ) %>%
  mutate(
    Km_mg_per_L =
      factor(
        Km_mg_per_L,
        levels =
          Km_sensitivity_values
      )
  ) %>%
  ggplot(
    aes(
      x = Km_mg_per_L,
      y = Combined_Vmax_mg_per_kg_per_h,
      group = Treatment_Type
    )
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(
      linetype = Treatment_Type
    )
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(
      shape = Treatment_Type
    ),
    size = 3
  ) +
  labs(
    x =
      expression(
        K[m]~"(mg O"[2]~L^-1*")"
      ),
    y =
      expression(
        V[max]~
          "(mg O"[2]~
          kg^-1~
          dry~sediment~
          h^-1*")"
      ),
    title =
      "Sensitivity of fitted Vmax to fixed Km"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_Vmax.pdf"
  ),
  p_sens_vmax,
  width = 8,
  height = 6
)


# ------------------------------------------------------------
# kL sensitivity
# ------------------------------------------------------------

p_sens_kL <- sensitivity_df %>%
  filter(
    Treatment_Type != "Unknown"
  ) %>%
  mutate(
    Km_mg_per_L =
      factor(
        Km_mg_per_L,
        levels =
          Km_sensitivity_values
      )
  ) %>%
  ggplot(
    aes(
      x = Km_mg_per_L,
      y = Combined_kL_per_h,
      group = Treatment_Type
    )
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(
      linetype = Treatment_Type
    )
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(
      shape = Treatment_Type
    ),
    size = 3
  ) +
  labs(
    x =
      expression(
        K[m]~"(mg O"[2]~L^-1*")"
      ),
    y =
      expression(
        k[L]~"(h"^-1*")"
      ),
    title =
      "Sensitivity of fitted kL to fixed Km"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "YEP_Km_Sensitivity_kL.pdf"
  ),
  p_sens_kL,
  width = 8,
  height = 6
)


# ============================================================
# Final messages
# ============================================================

message(
  "\n============================================"
)

message(
  "YEP MODELING COMPLETE"
)

message(
  "Primary Km = ",
  Km_primary,
  " mg O2/L"
)

message(
  "CB = ",
  CB,
  " (fixed dimensionless normalization factor)"
)

message(
  "Km sensitivity values = ",
  paste(
    Km_sensitivity_values,
    collapse = ", "
  ),
  " mg O2/L"
)

message(
  "Output directory: ",
  out_dir
)

message(
  "============================================"
)


print(
  sensitivity_treatment_summary
)

print(
  sensitivity_model_preference
)