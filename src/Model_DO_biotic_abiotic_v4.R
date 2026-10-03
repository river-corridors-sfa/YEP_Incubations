# ============================================================
# Complete YEP dissolved oxygen modeling analysis - v3
#
# v3 formulation retained:
#   - primary Km = 0.2 mg O2/L
#   - Km sensitivity = 0.1, 0.2, 0.5, 1.0 mg O2/L
#   - CB = 1, fixed dimensionless normalization factor
#   - first-order, biological, and combined models fit separately
#
# Changes in this version:
#   1. Best-fitting model selected by minimum SSR, following the
#      model-fit comparison used by Patel et al. (2024).
#   2. RMSE, AIC, residual autocorrelation, and bound flags are
#      retained as diagnostics only; they do not select the model.
#   3. Km sensitivity plots use treatment colors and sediment-history
#      line types/shapes.
#   4. A Figure-1-style dry-vs-wet plot is generated for ALL Km values,
#      with a free y-axis for every individual panel.
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
fig_dir <- file.path("Figures", "Km_sensitivity")

if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
}

if(!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figure directory: ", fig_dir)
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

# ORIGINAL v3 sensitivity range
Km_sensitivity_values <- c(
  0.13,
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.7,
  0.8,
  0.9,
  1.0,
  1.1,
  1.22
)

# Parameter bounds used during optimization
Vmax_lower_bound <- 0.01
Vmax_upper_bound <- 50
kL_lower_bound <- 0.001
kL_upper_bound <- 10

# Numerical tolerance used only to identify estimates sitting on a bound
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
# IMPORTANT v3 difference from Patel:
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
# dDO/dt = -(Vmax * DO/(DO + Km) * CB)
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
      
      r_bio <-
        Vmax *
        DO /
        (DO + Km)
      
      dDO <-
        -r_bio *
        CB
      
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
    
    # Do not substitute a different approximation if solver fails
    data.frame(
      time = times,
      DO = NA_real_
    )
  })
}


# ------------------------------------------------------------
# 3. Combined biological + first-order model
#
# dDO/dt = -(Vmax * DO/(DO + Km) * CB) - kL * DO
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
      
      r_bio <-
        Vmax *
        DO /
        (DO + Km)
      
      r_first_order <-
        kL *
        DO
      
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
# Information criteria - retained as diagnostics only
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
  
  if(
    is.null(pred) ||
    length(obs) != length(pred) ||
    all(!is.finite(pred))
  ) {
    return(NA_real_)
  }
  
  sqrt(
    mean(
      (obs - pred)^2,
      na.rm = TRUE
    )
  )
}


calc_lag1_residual_ac <- function(obs, pred) {
  
  if(
    is.null(pred) ||
    length(obs) != length(pred)
  ) {
    return(NA_real_)
  }
  
  residuals <- obs - pred
  residuals <- residuals[is.finite(residuals)]
  
  if(
    length(residuals) < 3 ||
    sd(residuals) == 0
  ) {
    return(NA_real_)
  }
  
  cor(
    residuals[-length(residuals)],
    residuals[-1],
    use = "complete.obs"
  )
}


at_lower_bound <- function(
    value,
    lower_bound,
    tolerance = bound_tolerance
) {
  
  if(
    length(value) == 0 ||
    is.na(value) ||
    !is.finite(value)
  ) {
    return(NA)
  }
  
  value <= (
    lower_bound + tolerance
  )
}


at_upper_bound <- function(
    value,
    upper_bound,
    tolerance = bound_tolerance
) {
  
  if(
    length(value) == 0 ||
    is.na(value) ||
    !is.finite(value)
  ) {
    return(NA)
  }
  
  value >= (
    upper_bound - tolerance
  )
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
        kL = fit_LIN$par[["kL"]]
      ),
      times = Data$time,
      DO0 = DO0
    )
  }
  
  if(!is.null(fit_BIO)) {
    
    model_BIO <- solveBiotic(
      list(
        Vmax = fit_BIO$par[["Vmax"]]
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
        Vmax = fit_COM$par[["Vmax"]],
        kL = fit_COM$par[["kL"]]
      ),
      times = Data$time,
      DO0 = DO0,
      CB = CB,
      Km = Km
    )
  }
  
  
  # ==========================================================
  # Model diagnostics
  #
  # These are retained for inspection but DO NOT select the model.
  # Model selection below is minimum SSR, following Patel et al.
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
      c(
        Linear_RMSE,
        Biotic_RMSE
      ),
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
  
  # First-order model: mg O2/L/h
  DO_rate_linear <- if(!is.null(fit_LIN)) {
    
    fit_LIN$par[["kL"]] *
      DO0
    
  } else {
    
    NA_real_
  }
  
  
  # Biological model: mg O2/L/h
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
  
  
  # Combined model components: mg O2/L/h
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
  
  
  # Relative model-attributed contributions at initial DO
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
  # mg/L/h x L/g dry sediment = mg/g dry sediment/h
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
  # Best-fitting model by minimum SSR - Patel-style comparison
  #
  # IMPORTANT:
  # SSR is used here only to identify which fitted configuration
  # reproduces the observed trajectory with the smallest residual
  # sum of squares. It does not penalize the combined model for
  # its additional parameter.
  # ==========================================================
  
  ssr_values <- c(
    FirstOrder = if(!is.null(fit_LIN)) fit_LIN$ssr else NA_real_,
    Biotic = if(!is.null(fit_BIO)) fit_BIO$ssr else NA_real_,
    Combined = if(!is.null(fit_COM)) fit_COM$ssr else NA_real_
  )
  
  valid_ssr <-
    !is.na(ssr_values) &
    is.finite(ssr_values)
  
  if(any(valid_ssr)) {
    
    Best_Model <-
      names(ssr_values[valid_ssr])[
        which.min(
          ssr_values[valid_ssr]
        )
      ]
    
    Best_SSR <-
      min(
        ssr_values[valid_ssr]
      )
    
  } else {
    
    Best_Model <- "None"
    Best_SSR <- NA_real_
  }
  
  
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
    
    
    # First-order model
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
          6
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
      round(
        Linear_RMSE,
        6
      ),
    
    Linear_Residual_Lag1_AC =
      round(
        Linear_Residual_Lag1_AC,
        6
      ),
    
    
    # Biological model
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
          6
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
      round(
        Biotic_RMSE,
        6
      ),
    
    Biotic_Residual_Lag1_AC =
      round(
        Biotic_Residual_Lag1_AC,
        6
      ),
    
    
    # Combined model
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
          6
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
      round(
        Combined_RMSE,
        6
      ),
    
    Combined_Residual_Lag1_AC =
      round(
        Combined_Residual_Lag1_AC,
        6
      ),
    
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
    
    Best_Model =
      Best_Model,
    
    Best_SSR =
      round(
        Best_SSR,
        6
      ),
    
    stringsAsFactors = FALSE
  )
  
  
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
# Evaluation plots
# ============================================================

create_evaluation_plots <- function(summary_df) {
  
  plot_data <- summary_df %>%
    filter(
      Treatment_Type != "Unknown"
    )
  
  
  # ----------------------------------------------------------
  # Model performance: SSR, matching Patel-style comparison
  # ----------------------------------------------------------
  
  pdf(
    file.path(
      out_dir,
      "Model_Performance_Summary.pdf"
    ),
    width = 12,
    height = 8
  )
  
  ssr_summary <- plot_data %>%
    group_by(
      Treatment_Type
    ) %>%
    summarise(
      
      FirstOrder_SSR =
        mean(
          Linear_SSR,
          na.rm = TRUE
        ),
      
      Biotic_SSR =
        mean(
          Biotic_SSR,
          na.rm = TRUE
        ),
      
      Combined_SSR =
        mean(
          Combined_SSR,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  ssr_matrix <- as.matrix(
    ssr_summary[
      ,
      c(
        "FirstOrder_SSR",
        "Biotic_SSR",
        "Combined_SSR"
      )
    ]
  )
  
  rownames(ssr_matrix) <-
    ssr_summary$Treatment_Type
  
  barplot(
    t(ssr_matrix),
    beside = TRUE,
    col = c(
      "lightblue",
      "lightcoral",
      "lightgreen"
    ),
    main =
      "Model Performance by Treatment\nLower SSR = Better Fit",
    ylab =
      "Sum of Squared Residuals",
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
        round(
          mean(
            Linear_SSR,
            na.rm = TRUE
          ),
          3
        ),
      
      Biotic_SSR =
        round(
          mean(
            Biotic_SSR,
            na.rm = TRUE
          ),
          3
        ),
      
      Combined_SSR =
        round(
          mean(
            Combined_SSR,
            na.rm = TRUE
          ),
          3
        ),
      
      FirstOrder_RMSE =
        round(
          mean(
            Linear_RMSE_mg_per_L,
            na.rm = TRUE
          ),
          4
        ),
      
      Biotic_RMSE =
        round(
          mean(
            Biotic_RMSE_mg_per_L,
            na.rm = TRUE
          ),
          4
        ),
      
      Combined_RMSE =
        round(
          mean(
            Combined_RMSE_mg_per_L,
            na.rm = TRUE
          ),
          4
        ),
      
      Combined_kL_at_lower_bound_n =
        sum(
          Combined_kL_at_lower_bound %in% TRUE
        ),
      
      Combined_Vmax_at_lower_bound_n =
        sum(
          Combined_Vmax_at_lower_bound %in% TRUE
        ),
      
      Most_Common_Best_Model = {
        
        x <- Best_Model[
          !is.na(Best_Model) &
            Best_Model != "None"
        ]
        
        if(length(x) == 0) {
          
          NA_character_
          
        } else {
          
          names(
            which.max(
              table(x)
            )
          )
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
  # Overall best-fitting model by minimum SSR
  # ----------------------------------------------------------
  
  model_preference <- summary_df %>%
    filter(
      Treatment_Type != "Unknown"
    ) %>%
    count(
      Best_Model,
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
  
  write.csv(
    model_preference,
    file.path(
      out_dir,
      "YEP_Best_Fit_Model_Preference.csv"
    ),
    row.names = FALSE
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
# ORIGINAL v3 RANGE:
# 0.1 mg/L = lower bracketing scenario
# 0.2 mg/L = primary value
# 0.5 mg/L = moderate scenario
# 1.0 mg/L = upper scenario
#
# Purpose:
# Refit the entire model at each fixed Km and assess sensitivity
# of Vmax, kL, overall fitted rate, and best-fitting configuration.
#
# Km is NOT selected by whichever value gives the lowest SSR.
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
    
    sensitivity_results_list[[
      as.character(km_test)
    ]] <- sensitivity_summary_km
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
    
    Reference_Best_Model =
      Best_Model
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
        
        NA_real_
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
        
        NA_real_
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
        
        NA_real_
      ),
    
    Best_Model_Changed =
      Best_Model !=
      Reference_Best_Model
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
    
    Mean_Combined_SSR =
      mean(
        Combined_SSR,
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
# Best-fitting configuration across Km values
# ============================================================

sensitivity_model_preference <- sensitivity_df %>%
  filter(
    Treatment_Type != "Unknown"
  ) %>%
  count(
    Km_mg_per_L,
    Best_Model
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
# Prepare sensitivity data for plots
# ============================================================

# Treatment colors match the manuscript treatment plots
treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

# Sediment-history colors for Figure 1 style plot
sediment_colors <- c(
  "Dry Sediments" = "#56B4E9",
  "Wet Sediments" = "#1B4F72"
)

# Factor ordering
treatment_order <- c(
  "Control",
  "Unburned + DOC",
  "High Burn + DOC"
)

sediment_order <- c(
  "Dry Sediments",
  "Wet Sediments"
)

km_label_order <- paste0(
  "Km = ",
  Km_sensitivity_values
)


sensitivity_plot_df <- sensitivity_df %>%
  filter(
    Treatment_Type != "Unknown"
  ) %>%
  mutate(
    
    Treatment = case_when(
      grepl(
        "Control",
        Treatment_Type
      ) ~ "Control",
      grepl(
        "Unburned",
        Treatment_Type
      ) ~ "Unburned + DOC",
      grepl(
        "HighBurn",
        Treatment_Type
      ) ~ "High Burn + DOC",
      TRUE ~ NA_character_
    ),
    
    Sediment_Type = case_when(
      grepl(
        "^Dry_",
        Treatment_Type
      ) ~ "Dry Sediments",
      grepl(
        "^Wet_",
        Treatment_Type
      ) ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    
    Treatment = factor(
      Treatment,
      levels = treatment_order
    ),
    
    Sediment_Type = factor(
      Sediment_Type,
      levels = sediment_order
    ),
    
    Km_Factor = factor(
      Km_mg_per_L,
      levels = Km_sensitivity_values
    ),
    
    Km_Label = factor(
      paste0(
        "Km = ",
        Km_mg_per_L
      ),
      levels = km_label_order
    )
  ) %>%
  filter(
    !is.na(Treatment),
    !is.na(Sediment_Type)
  )


# ============================================================
# Prepare summary statistics for sensitivity plots
#
# Raw points = individual replicates
# Large points/lines = group means
# Error bars = 95% CI of the mean across replicates
# ============================================================

vmax_sensitivity_summary <- sensitivity_plot_df %>%
  group_by(
    Km_Factor,
    Treatment,
    Sediment_Type
  ) %>%
  summarise(
    n = sum(!is.na(Combined_Vmax_mg_per_kg_per_h)),
    
    mean = mean(
      Combined_Vmax_mg_per_kg_per_h,
      na.rm = TRUE
    ),
    
    sd = sd(
      Combined_Vmax_mg_per_kg_per_h,
      na.rm = TRUE
    ),
    
    se = sd / sqrt(n),
    
    t_crit = ifelse(
      n > 1,
      qt(0.975, df = n - 1),
      NA_real_
    ),
    
    ci_lower = mean - t_crit * se,
    ci_upper = mean + t_crit * se,
    
    .groups = "drop"
  )


kl_sensitivity_summary <- sensitivity_plot_df %>%
  group_by(
    Km_Factor,
    Treatment,
    Sediment_Type
  ) %>%
  summarise(
    n = sum(!is.na(Combined_kL_per_h)),
    
    mean = mean(
      Combined_kL_per_h,
      na.rm = TRUE
    ),
    
    sd = sd(
      Combined_kL_per_h,
      na.rm = TRUE
    ),
    
    se = sd / sqrt(n),
    
    t_crit = ifelse(
      n > 1,
      qt(0.975, df = n - 1),
      NA_real_
    ),
    
    ci_lower = mean - t_crit * se,
    ci_upper = mean + t_crit * se,
    
    .groups = "drop"
  )


# Same horizontal dodge for means, lines, and error bars
summary_dodge <- position_dodge(
  width = 0.45
)


# ============================================================
# Sensitivity plot 1: Vmax across Km
#
# Small points = all individual replicates
# Large points = mean
# Error bars = 95% confidence interval
#
# Color = treatment
# Line type / shape = dry vs wet history
# ============================================================

p_sens_vmax <- ggplot() +
  
  # ----------------------------------------------------------
# Individual replicate points
# ----------------------------------------------------------

geom_point(
  data = sensitivity_plot_df,
  aes(
    x = Km_Factor,
    y = Combined_Vmax_mg_per_kg_per_h,
    color = Treatment,
    shape = Sediment_Type,
    group = interaction(
      Treatment,
      Sediment_Type
    )
  ),
  position = position_jitterdodge(
    jitter.width = 0.06,
    jitter.height = 0,
    dodge.width = 0.45
  ),
  size = 1.8,
  alpha = 0.45
) +
  
  
  # ----------------------------------------------------------
# Mean lines
# ----------------------------------------------------------

geom_line(
  data = vmax_sensitivity_summary,
  aes(
    x = Km_Factor,
    y = mean,
    group = interaction(
      Treatment,
      Sediment_Type
    ),
    color = Treatment,
    linetype = Sediment_Type
  ),
  position = summary_dodge,
  linewidth = 1
) +
  
  
  # ----------------------------------------------------------
# 95% confidence intervals
# ----------------------------------------------------------

geom_errorbar(
  data = vmax_sensitivity_summary,
  aes(
    x = Km_Factor,
    ymin = ci_lower,
    ymax = ci_upper,
    color = Treatment,
    group = interaction(
      Treatment,
      Sediment_Type
    )
  ),
  position = summary_dodge,
  width = 0.12,
  linewidth = 0.7
) +
  
  
  # ----------------------------------------------------------
# Mean points
# ----------------------------------------------------------

geom_point(
  data = vmax_sensitivity_summary,
  aes(
    x = Km_Factor,
    y = mean,
    color = Treatment,
    shape = Sediment_Type,
    group = interaction(
      Treatment,
      Sediment_Type
    )
  ),
  position = summary_dodge,
  size = 3
) +
  
  
  scale_color_manual(
    values = treatment_colors
  ) +
  
  scale_linetype_manual(
    values = c(
      "Dry Sediments" = "solid",
      "Wet Sediments" = "dashed"
    )
  ) +
  
  labs(
    x = expression(
      K[m]~"(mg O"[2]~L^-1*")"
    ),
    y = expression(
      V[max]~
        "(mg O"[2]~
        kg^-1~
        dry~sediment~
        h^-1*")"
    ),
    color = "Treatment",
    linetype = "Sediment Type",
    shape = "Sediment Type"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


print(
  p_sens_vmax
)


ggsave(
  file.path(
    fig_dir,
    "YEP_Km_Sensitivity_Vmax.png"
  ),
  p_sens_vmax,
  width = 9,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    fig_dir,
    "YEP_Km_Sensitivity_Vmax.pdf"
  ),
  p_sens_vmax,
  width = 9,
  height = 6
)


# ============================================================
# Sensitivity plot 2: kL across Km
#
# Small points = all individual replicates
# Large points = mean
# Error bars = 95% confidence interval
# ============================================================

p_sens_kL <- ggplot() +
  
  # ----------------------------------------------------------
# Individual replicate points
# ----------------------------------------------------------

geom_point(
  data = sensitivity_plot_df,
  aes(
    x = Km_Factor,
    y = Combined_kL_per_h,
    color = Treatment,
    shape = Sediment_Type,
    group = interaction(
      Treatment,
      Sediment_Type
    )
  ),
  position = position_jitterdodge(
    jitter.width = 0.06,
    jitter.height = 0,
    dodge.width = 0.45
  ),
  size = 1.8,
  alpha = 0.45
) +
  
  
  # ----------------------------------------------------------
# Mean lines
# ----------------------------------------------------------

geom_line(
  data = kl_sensitivity_summary,
  aes(
    x = Km_Factor,
    y = mean,
    group = interaction(
      Treatment,
      Sediment_Type
    ),
    color = Treatment,
    linetype = Sediment_Type
  ),
  position = summary_dodge,
  linewidth = 1
) +
  
  
  # ----------------------------------------------------------
# 95% confidence intervals
# ----------------------------------------------------------

geom_errorbar(
  data = kl_sensitivity_summary,
  aes(
    x = Km_Factor,
    ymin = ci_lower,
    ymax = ci_upper,
    color = Treatment,
    group = interaction(
      Treatment,
      Sediment_Type
    )
  ),
  position = summary_dodge,
  width = 0.12,
  linewidth = 0.7
) +
  
  
  # ----------------------------------------------------------
# Mean points
# ----------------------------------------------------------

geom_point(
  data = kl_sensitivity_summary,
  aes(
    x = Km_Factor,
    y = mean,
    color = Treatment,
    shape = Sediment_Type,
    group = interaction(
      Treatment,
      Sediment_Type
    )
  ),
  position = summary_dodge,
  size = 3
) +
  
  
  scale_color_manual(
    values = treatment_colors
  ) +
  
  scale_linetype_manual(
    values = c(
      "Dry Sediments" = "solid",
      "Wet Sediments" = "dashed"
    )
  ) +
  
  labs(
    x = expression(
      K[m]~"(mg O"[2]~L^-1*")"
    ),
    y = expression(
      k[L]~"(h"^-1*")"
    ),
    color = "Treatment",
    linetype = "Sediment Type",
    shape = "Sediment Type"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


print(
  p_sens_kL
)


ggsave(
  file.path(
    fig_dir,
    "YEP_Km_Sensitivity_kL.png"
  ),
  p_sens_kL,
  width = 9,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    fig_dir,
    "YEP_Km_Sensitivity_kL.pdf"
  ),
  p_sens_kL,
  width = 9,
  height = 6
)
# ============================================================
# FIGURE 1-STYLE PLOT FOR ALL Km VALUES
#
# Dry vs wet sediments for Vmax and kL.
#
# IMPORTANT:
# facet_wrap(..., scales = "free_y") gives every Km/parameter
# panel an independent y-axis range.
# ============================================================

fig1_all_km_data <- sensitivity_plot_df %>%
  select(
    Sample_Name,
    Km_mg_per_L,
    Km_Label,
    Treatment,
    Sediment_Type,
    Combined_Vmax_mg_per_kg_per_h,
    Combined_kL_per_h
  ) %>%
  pivot_longer(
    cols = c(
      Combined_Vmax_mg_per_kg_per_h,
      Combined_kL_per_h
    ),
    names_to = "Parameter",
    values_to = "Fitted_Value"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      Combined_Vmax_mg_per_kg_per_h = "Vmax",
      Combined_kL_per_h = "kL"
    ),
    Parameter = factor(
      Parameter,
      levels = c(
        "Vmax",
        "kL"
      )
    )
  ) %>%
  filter(
    !is.na(Fitted_Value)
  )


# ------------------------------------------------------------
# Wilcoxon dry vs wet at each Km and for each fitted parameter
# ------------------------------------------------------------

fig1_all_km_stats <- fig1_all_km_data %>%
  group_by(
    Km_mg_per_L,
    Km_Label,
    Parameter
  ) %>%
  summarise(
    
    n_dry =
      sum(
        Sediment_Type ==
          "Dry Sediments"
      ),
    
    n_wet =
      sum(
        Sediment_Type ==
          "Wet Sediments"
      ),
    
    p = if(
      n_dry > 0 &&
      n_wet > 0
    ) {
      
      wilcox.test(
        Fitted_Value ~
          Sediment_Type,
        exact = FALSE
      )$p.value
      
    } else {
      
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p) ~ "Wilcox p = NA",
      p < 0.001 ~ "Wilcox p < 0.001",
      TRUE ~ sprintf(
        "Wilcox p = %.3f",
        p
      )
    )
  )


# Position p-value independently within every facet
fig1_all_km_positions <- fig1_all_km_data %>%
  group_by(
    Km_mg_per_L,
    Km_Label,
    Parameter
  ) %>%
  summarise(
    y_min = min(
      Fitted_Value,
      na.rm = TRUE
    ),
    y_max = max(
      Fitted_Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    y_range =
      y_max -
      y_min,
    y_range = ifelse(
      y_range == 0,
      pmax(
        abs(y_max),
        1
      ),
      y_range
    ),
    y =
      y_max +
      0.18 *
      y_range,
    x = 1.5
  )


fig1_all_km_stats_plot <- fig1_all_km_stats %>%
  left_join(
    fig1_all_km_positions,
    by = c(
      "Km_mg_per_L",
      "Km_Label",
      "Parameter"
    )
  )


write.csv(
  fig1_all_km_stats,
  file.path(
    fig_dir,
    "Figure1_All_Km_DryWet_Wilcox_stats.csv"
  ),
  row.names = FALSE
)


# Parameter strip labels
parameter_labels <- c(
  "Vmax" = "Biotic O2 consumption at saturation (Vmax)",
  "kL" = "First-order O2-consumption coefficient (kL)"
)


fig1_all_km <- ggplot(
  fig1_all_km_data,
  aes(
    x = Sediment_Type,
    y = Fitted_Value,
    fill = Sediment_Type
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(
      width = 0.10
    ),
    size = 1.8,
    alpha = 0.75
  ) +
  geom_blank(
    data = fig1_all_km_stats_plot,
    aes(
      x = x,
      y = y
    ),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = fig1_all_km_stats_plot,
    aes(
      x = x,
      y = y,
      label = p_label
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.2
  ) +
  scale_fill_manual(
    values = sediment_colors
  ) +
  scale_x_discrete(
    labels = c(
      "Dry Sediments" = "Dry",
      "Wet Sediments" = "Wet"
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Sediment Type"
  ) +
  facet_wrap(
    vars(
      Parameter,
      Km_Label
    ),
    nrow = 2,
    scales = "free_y",
    labeller = labeller(
      Parameter = parameter_labels
    )
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      size = 10
    ),
    panel.spacing = unit(
      0.8,
      "lines"
    )
  )


print(
  fig1_all_km
)


ggsave(
  file.path(
    fig_dir,
    "Figure1_All_Km_DryWet_freeY.png"
  ),
  fig1_all_km,
  width = 13,
  height = 7,
  dpi = 300
)


ggsave(
  file.path(
    fig_dir,
    "Figure1_All_Km_DryWet_freeY.pdf"
  ),
  fig1_all_km,
  width = 13,
  height = 7
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
  "Best-fitting model = minimum SSR"
)

message(
  "Model output directory: ",
  out_dir
)

message(
  "Sensitivity figure directory: ",
  fig_dir
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
