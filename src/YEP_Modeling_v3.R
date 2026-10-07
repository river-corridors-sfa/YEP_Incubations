# =============================================================================
# YEP sediment O2 consumption: combined biotic + abiotic model
# Version 3 (revision) -- Km sensitivity analysis
# -----------------------------------------------------------------------------
# Model (single combined model; abiotic-only and biotic-only are nested cases):
#
#   dDO/dt = -kL * DO  -  Vmax * DO / (Km + DO)
#
#   kL   : first-order abiotic rate constant (h^-1)
#   Vmax : maximum biotic O2 consumption rate (mg O2 L^-1 h^-1)
#   Km   : half-saturation constant for O2 (mg L^-1), FIXED (not estimated)
#   DO0  : initial DO (mg L^-1), ESTIMATED with kL and Vmax
#
# Changes from v2 (Model_DO_biotic_abiotic_v2.R):
#   1. Km = 0.2 mg/L base (Gu et al. 2007; within Sihi et al. 2020 posterior,
#      0.13-1.22 mg/L), replacing Km = 10 mg/L. Km varied in two grids.
#   2. Biomass state variable CB removed. In v2, CB grew by the O2 consumed
#      (dCB/dt = r1*CB), reaching 2-8x its initial value within each run.
#      Vmax is now in mg O2 L^-1 h^-1.
#   3. Reads the corrected merged DO file (includes H1/U1 start-time fixes).
#   4. Anoxic plateau removed: series cut before the first DO <= 0.1 mg/L.
#   5. End-of-run DO rise removed: series truncated at the DO minimum.
#   6. No silent fallback solutions when the ODE solver fails.
#   7. Three starting guesses per fit; agreement among them is recorded.
#   8. Bounds start at 0 so kL or Vmax can reach zero; bound hits are flagged.
#   9. DO0 estimated rather than fixed to the first observation.
#  10. Water-only controls (YEP_INC-W*) excluded explicitly.
#  11. Only the combined model is fit (no AIC-based model selection).
#      AIC/BIC are reported as descriptive metrics only (5-s data are
#      autocorrelated; see lag-1 residual autocorrelation column).
#  12. First 2 min and last 2 min trimmed (unchanged from v2).
# =============================================================================

rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  library(deSolve)
  library(minpack.lm)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
do_file <- file.path("v2_data", "Merged_DO_Firesting_Data_2026-10-07.csv")
mass_volume_file <- file.path("v2_data", "v2_YEP_Sample_Data",
                              "v2_YEP_Sediment_Water_Mass_Volume.csv")
out_dir <- "modeling_outputs_v3"

Km_base  <- 0.2                                                   # mg/L
Km_gridA <- c(0.05, 0.1, 0.13, 0.2, 0.3, 0.51, 0.75, 1.0, 1.22)   # mg/L
Km_gridB <- c(Km_gridA, 3.2, 5.28, 10)                            # mg/L

trim_start_min   <- 2      # minutes removed at start
trim_end_min     <- 2      # minutes removed at end
anoxic_threshold <- 0.1    # mg/L; series cut before first value <= this
DO_reference     <- 5      # mg/L; rates also reported at this common DO

bounds_lower <- c(kL = 0,  Vmax = 0,   DO0 = 0)
bounds_upper <- c(kL = 50, Vmax = 500, DO0 = 15)

save_fit_qc_plot <- TRUE   # observed vs fitted DO for all vials at Km_base

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load and label data
# -----------------------------------------------------------------------------
df_raw <- read.csv(do_file, stringsAsFactors = FALSE)

df <- df_raw %>%
  filter(!grepl("^YEP_INC-W", Sample_Name)) %>%          # water-only controls
  mutate(Elapsed_Seconds = as.numeric(Elapsed_Seconds),
         DO_mg_per_L     = as.numeric(DO_mg_per_L)) %>%
  filter(!is.na(Elapsed_Seconds), !is.na(DO_mg_per_L)) %>%
  mutate(
    Moisture = case_when(grepl("^YEP1", Sample_Name) ~ "Dry",
                         grepl("^YEP2", Sample_Name) ~ "Wet",
                         TRUE ~ NA_character_),
    Treatment_Code = sub(".*_INC-([A-Za-z]).*", "\\1", Sample_Name),
    Treatment = recode(Treatment_Code,
                       S = "Control", U = "Unburned DOM", H = "Burned DOM",
                       .default = NA_character_),
    Incubation_Set = sub("^YEP[12]([A-Z]).*", "\\1", Sample_Name)
  )

if (any(is.na(df$Moisture)) || any(is.na(df$Treatment))) {
  stop("Unrecognized sample names: ",
       paste(unique(df$Sample_Name[is.na(df$Moisture) | is.na(df$Treatment)]),
             collapse = ", "))
}
df$Treatment <- factor(df$Treatment,
                       levels = c("Control", "Unburned DOM", "Burned DOM"))

# Mass/volume conversion (L water per g dry sediment); optional
if (file.exists(mass_volume_file)) {
  conversion_factors <- read.csv(mass_volume_file, stringsAsFactors = FALSE,
                                 skip = 2) %>%
    filter(grepl("_INC-", Sample_Name)) %>%
    mutate(Dry_Sediment_Mass_g = as.numeric(Dry_Sediment_Mass_g),
           Water_Volume_L      = as.numeric(Water_Mass_g) / 1000,
           Conversion_Factor   = Water_Volume_L / Dry_Sediment_Mass_g) %>%
    select(Sample_Name, Dry_Sediment_Mass_g, Water_Volume_L, Conversion_Factor)
} else {
  warning("Mass/volume file not found; mass-normalized rates will be NA.")
  conversion_factors <- data.frame(Sample_Name = character(),
                                   Dry_Sediment_Mass_g = numeric(),
                                   Water_Volume_L = numeric(),
                                   Conversion_Factor = numeric())
}

# -----------------------------------------------------------------------------
# Data preparation per vial
# -----------------------------------------------------------------------------
prepare_vial <- function(d) {
  d <- d %>% arrange(Elapsed_Seconds)
  t_hr <- (d$Elapsed_Seconds - min(d$Elapsed_Seconds)) / 3600
  y    <- d$DO_mg_per_L
  n_raw <- length(y)
  
  # 1. Trim first and last 2 min
  keep <- t_hr >= trim_start_min / 60 & t_hr <= max(t_hr) - trim_end_min / 60
  t_hr <- t_hr[keep]; y <- y[keep]
  n_after_trim <- length(y)
  
  # 2. Truncate at the DO minimum (removes end-of-run DO rise)
  i_min <- which.min(y)
  n_removed_after_min <- length(y) - i_min
  t_hr <- t_hr[seq_len(i_min)]; y <- y[seq_len(i_min)]
  
  # 3. Cut before first DO <= anoxic threshold (removes calibration plateau)
  i_anox <- which(y <= anoxic_threshold)[1]
  reached_anoxia <- !is.na(i_anox)
  n_removed_anoxic <- 0
  if (reached_anoxia) {
    n_removed_anoxic <- length(y) - (i_anox - 1)
    t_hr <- t_hr[seq_len(i_anox - 1)]; y <- y[seq_len(i_anox - 1)]
  }
  
  # Rebase time to 0
  t_hr <- t_hr - t_hr[1]
  
  list(
    data = data.frame(time = t_hr, DO = y),
    log = data.frame(
      n_raw = n_raw,
      n_after_trim = n_after_trim,
      n_removed_after_min = n_removed_after_min,
      reached_anoxia = reached_anoxia,
      n_removed_anoxic = n_removed_anoxic,
      n_fit = length(y),
      duration_fit_min = round(60 * max(t_hr), 2),
      DO_first_obs = y[1],
      DO_last_obs = y[length(y)]
    )
  )
}

# -----------------------------------------------------------------------------
# Model
# -----------------------------------------------------------------------------
derivs <- function(t, state, p) {
  DO <- max(state[["DO"]], 0)
  list(-p[["kL"]] * DO - p[["Vmax"]] * DO / (p[["Km"]] + DO))
}

solve_model <- function(kL, Vmax, DO0, Km, times) {
  out <- ode(y = c(DO = DO0), times = times, func = derivs,
             parms = c(kL = kL, Vmax = Vmax, Km = Km),
             method = "lsoda", rtol = 1e-8, atol = 1e-10)
  if (nrow(out) != length(times) || any(!is.finite(out[, "DO"]))) {
    return(rep(NA_real_, length(times)))
  }
  out[, "DO"]
}

# -----------------------------------------------------------------------------
# Fit one vial at one Km
# -----------------------------------------------------------------------------
fit_vial <- function(dat, Km) {
  times <- dat$time
  obs   <- dat$DO
  
  resid_fn <- function(p) {
    pred <- solve_model(p[["kL"]], p[["Vmax"]], p[["DO0"]], Km, times)
    r <- obs - pred
    r[!is.finite(r)] <- 1e3      # penalize solver failure; not hidden
    r
  }
  
  # Starting guesses
  slope <- abs(coef(lm(obs ~ times))[[2]])   # mg/L/h
  DOm   <- mean(obs)
  starts <- list(
    biotic_dominant  = c(kL = 0.01,          Vmax = slope,     DO0 = obs[1]),
    abiotic_dominant = c(kL = slope / DOm,   Vmax = 0.1,       DO0 = obs[1]),
    mixed            = c(kL = slope / DOm / 2, Vmax = slope / 2, DO0 = obs[1])
  )
  
  run_lm <- function(s) {
    s <- pmin(pmax(s, bounds_lower + 1e-6), bounds_upper - 1e-6)
    tryCatch(
      nls.lm(par = s, fn = resid_fn,
             lower = bounds_lower, upper = bounds_upper,
             control = nls.lm.control(maxiter = 300, ftol = 1e-12,
                                      ptol = 1e-12)),
      error = function(e) NULL)
  }
  
  # Each start is fit, then refit once from its own result ("polish").
  # Convergence slows when a parameter sits on its bound (e.g. kL = 0 in
  # dry vials); the polish step removes that small stopping error.
  fits <- lapply(starts, function(s) {
    f1 <- run_lm(s)
    if (is.null(f1)) return(NULL)
    f2 <- run_lm(f1$par)
    if (!is.null(f2) && sum(f2$fvec^2) <= sum(f1$fvec^2)) f2 else f1
  })
  
  ssr <- sapply(fits, function(f) if (is.null(f)) Inf else sum(f$fvec^2))
  if (all(!is.finite(ssr))) {
    return(data.frame(fit_ok = FALSE))
  }
  best <- fits[[which.min(ssr)]]
  p <- best$par
  
  # Agreement among starts: SSR within 1% of the best, and kL and Vmax
  # within 2% (or 0.05 absolute, for values near zero) of the best fit
  agree <- sapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    if (is.null(f) || !is.finite(ssr[i])) return(FALSE)
    abs(ssr[i] - min(ssr)) <= 0.01 * max(min(ssr), 1e-12) &&
      all(abs(f$par[c("kL", "Vmax")] - p[c("kL", "Vmax")]) <=
            pmax(0.02 * abs(p[c("kL", "Vmax")]), 0.05))
  })
  
  pred <- solve_model(p[["kL"]], p[["Vmax"]], p[["DO0"]], Km, times)
  res  <- obs - pred
  n <- length(obs); k <- 3
  SSR <- sum(res^2)
  lag1 <- if (n > 2) cor(res[-1], res[-n]) else NA_real_
  
  data.frame(
    fit_ok = TRUE,
    kL_per_h = p[["kL"]],
    Vmax_mg_L_h = p[["Vmax"]],
    DO0_est = p[["DO0"]],
    # error metrics
    n_points = n,
    SSR = SSR,
    RMSE = sqrt(SSR / n),
    MAE = mean(abs(res)),
    R2 = 1 - SSR / sum((obs - mean(obs))^2),
    Mean_Bias = mean(pred - obs),
    Max_Abs_Resid = max(abs(res)),
    Resid_Lag1_Autocorr = lag1,
    AIC_descriptive = n * log(SSR / n) + 2 * k,
    BIC_descriptive = n * log(SSR / n) + k * log(n),
    # diagnostics
    n_starts_agree = sum(agree),
    n_starts = length(starts),
    kL_at_lower = p[["kL"]] < 1e-3,
    kL_at_upper = p[["kL"]] > 0.99 * bounds_upper[["kL"]],
    Vmax_at_lower = p[["Vmax"]] < 1e-2,
    Vmax_at_upper = p[["Vmax"]] > 0.99 * bounds_upper[["Vmax"]],
    nls_info = best$info
  )
}

# -----------------------------------------------------------------------------
# Run all vials x all Km
# -----------------------------------------------------------------------------
samples <- sort(unique(df$Sample_Name))
Km_all  <- sort(unique(Km_gridB))

prep <- lapply(samples, function(s) prepare_vial(df %>% filter(Sample_Name == s)))
names(prep) <- samples

meta <- df %>% distinct(Sample_Name, Moisture, Treatment, Incubation_Set)

prep_log <- bind_rows(lapply(samples, function(s)
  cbind(Sample_Name = s, prep[[s]]$log))) %>%
  left_join(meta, by = "Sample_Name") %>%
  relocate(Moisture, Treatment, Incubation_Set, .after = Sample_Name)
write.csv(prep_log, file.path(out_dir, "YEP_DataPrep_Log.csv"), row.names = FALSE)

message("Fitting ", length(samples), " vials x ", length(Km_all), " Km values")
results <- list()
for (Km in Km_all) {
  message("  Km = ", Km)
  for (s in samples) {
    fr <- fit_vial(prep[[s]]$data, Km)
    results[[length(results) + 1]] <- cbind(Sample_Name = s, Km_mg_L = Km, fr)
  }
}

per_vial <- bind_rows(results) %>%
  left_join(meta, by = "Sample_Name") %>%
  left_join(conversion_factors, by = "Sample_Name") %>%
  mutate(
    # Rates at the fitted initial DO
    Abiotic_Rate_DO0_mg_L_h = kL_per_h * DO0_est,
    Biotic_Rate_DO0_mg_L_h  = Vmax_mg_L_h * DO0_est / (Km_mg_L + DO0_est),
    Total_Rate_DO0_mg_L_h   = Abiotic_Rate_DO0_mg_L_h + Biotic_Rate_DO0_mg_L_h,
    Abiotic_Fraction_DO0    = Abiotic_Rate_DO0_mg_L_h / Total_Rate_DO0_mg_L_h,
    # Rates at a common reference DO (removes differences in starting DO)
    Abiotic_Rate_DOref_mg_L_h = kL_per_h * DO_reference,
    Biotic_Rate_DOref_mg_L_h  = Vmax_mg_L_h * DO_reference / (Km_mg_L + DO_reference),
    Total_Rate_DOref_mg_L_h   = Abiotic_Rate_DOref_mg_L_h + Biotic_Rate_DOref_mg_L_h,
    # Mass-normalized (mg O2 per kg dry sediment per h), at fitted initial DO
    Abiotic_Rate_DO0_mg_kg_h = Abiotic_Rate_DO0_mg_L_h * Conversion_Factor * 1000,
    Biotic_Rate_DO0_mg_kg_h  = Biotic_Rate_DO0_mg_L_h  * Conversion_Factor * 1000,
    Total_Rate_DO0_mg_kg_h   = Total_Rate_DO0_mg_L_h   * Conversion_Factor * 1000,
    In_GridA = Km_mg_L %in% Km_gridA,
    In_GridB = Km_mg_L %in% Km_gridB
  ) %>%
  relocate(Moisture, Treatment, Incubation_Set, .after = Sample_Name)

write.csv(per_vial, file.path(out_dir, "YEP_KmSensitivity_PerVial_AllKm.csv"),
          row.names = FALSE)
write.csv(per_vial %>% filter(Km_mg_L == Km_base),
          file.path(out_dir, paste0("YEP_CombinedModel_Km", Km_base, "_PerVial.csv")),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------
summarise_groups <- function(x) {
  x %>%
    group_by(Km_mg_L, Moisture, Treatment) %>%
    summarise(
      n = n(),
      kL_mean = mean(kL_per_h), kL_sd = sd(kL_per_h),
      Vmax_mean = mean(Vmax_mg_L_h), Vmax_sd = sd(Vmax_mg_L_h),
      Total_Rate_DO0_mean = mean(Total_Rate_DO0_mg_L_h),
      Abiotic_Fraction_DO0_mean = mean(Abiotic_Fraction_DO0),
      .groups = "drop")
}

summarise_errors <- function(x) {
  x %>%
    group_by(Km_mg_L, Moisture) %>%
    summarise(
      n_vials = n(),
      SSR_mean = mean(SSR), SSR_median = median(SSR),
      RMSE_mean = mean(RMSE), RMSE_max = max(RMSE),
      MAE_mean = mean(MAE),
      R2_mean = mean(R2), R2_min = min(R2),
      Mean_Bias_mean = mean(Mean_Bias),
      Max_Abs_Resid_max = max(Max_Abs_Resid),
      Resid_Lag1_Autocorr_mean = mean(Resid_Lag1_Autocorr),
      AIC_descriptive_mean = mean(AIC_descriptive),
      n_kL_at_lower = sum(kL_at_lower),
      n_Vmax_at_lower = sum(Vmax_at_lower),
      n_any_upper_bound = sum(kL_at_upper | Vmax_at_upper),
      n_starts_disagree = sum(n_starts_agree < n_starts),
      .groups = "drop")
}

# Dry vs wet comparison per Km (vials as replicates; Wilcoxon rank-sum)
compare_dry_wet <- function(x) {
  x %>%
    group_by(Km_mg_L) %>%
    summarise(
      kL_dry = mean(kL_per_h[Moisture == "Dry"]),
      kL_wet = mean(kL_per_h[Moisture == "Wet"]),
      kL_p = wilcox.test(kL_per_h ~ Moisture, exact = FALSE)$p.value,
      Vmax_dry = mean(Vmax_mg_L_h[Moisture == "Dry"]),
      Vmax_wet = mean(Vmax_mg_L_h[Moisture == "Wet"]),
      Vmax_p = wilcox.test(Vmax_mg_L_h ~ Moisture, exact = FALSE)$p.value,
      TotalRate_DOref_dry = mean(Total_Rate_DOref_mg_L_h[Moisture == "Dry"]),
      TotalRate_DOref_wet = mean(Total_Rate_DOref_mg_L_h[Moisture == "Wet"]),
      TotalRate_DOref_p = wilcox.test(Total_Rate_DOref_mg_L_h ~ Moisture,
                                      exact = FALSE)$p.value,
      .groups = "drop") %>%
    mutate(kL_direction   = ifelse(kL_wet > kL_dry, "Wet > Dry", "Dry > Wet"),
           Vmax_direction = ifelse(Vmax_wet > Vmax_dry, "Wet > Dry", "Dry > Wet"))
}

for (g in c("A", "B")) {
  x <- per_vial %>% filter(fit_ok, if (g == "A") In_GridA else In_GridB)
  write.csv(summarise_groups(x),
            file.path(out_dir, paste0("YEP_KmSensitivity_GroupSummary_Grid", g, ".csv")),
            row.names = FALSE)
  write.csv(summarise_errors(x),
            file.path(out_dir, paste0("YEP_KmSensitivity_ErrorMetrics_Grid", g, ".csv")),
            row.names = FALSE)
  write.csv(compare_dry_wet(x),
            file.path(out_dir, paste0("YEP_KmSensitivity_DryVsWet_Grid", g, ".csv")),
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
trt_cols <- c("Control" = "grey40", "Unburned DOM" = "#1b9e77",
              "Burned DOM" = "#d95f02")
moist_shapes <- c("Dry" = 17, "Wet" = 16)

plot_parameters_vs_Km <- function(x, title_suffix) {
  long <- x %>%
    select(Sample_Name, Moisture, Treatment, Km_mg_L, kL_per_h, Vmax_mg_L_h) %>%
    pivot_longer(c(kL_per_h, Vmax_mg_L_h), names_to = "Parameter",
                 values_to = "Value") %>%
    mutate(Parameter = recode(Parameter,
                              kL_per_h = "kL (h^-1)",
                              Vmax_mg_L_h = "Vmax (mg O2 L^-1 h^-1)"))
  means <- long %>%
    group_by(Parameter, Km_mg_L, Moisture, Treatment) %>%
    summarise(Value = mean(Value), .groups = "drop")
  
  ggplot(long, aes(Km_mg_L, Value, colour = Treatment, shape = Moisture)) +
    annotate("rect", xmin = 0.13, xmax = 1.22, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_vline(xintercept = Km_base, linetype = "dashed", colour = "grey30") +
    geom_point(aes(group = interaction(Moisture, Treatment)),
               position = position_dodge(width = 0.05), alpha = 0.35, size = 1.6) +
    geom_line(data = means,
              aes(group = interaction(Moisture, Treatment), linetype = Moisture),
              linewidth = 0.7) +
    geom_point(data = means, size = 3) +
    scale_x_log10(breaks = sort(unique(x$Km_mg_L)),
                  labels = function(b) format(b, drop0trailing = TRUE)) +
    scale_colour_manual(values = trt_cols) +
    scale_shape_manual(values = moist_shapes) +
    scale_linetype_manual(values = c("Dry" = "dashed", "Wet" = "solid")) +
    facet_wrap(~Parameter, ncol = 1, scales = "free_y") +
    labs(x = "Km (mg L^-1, log scale)", y = NULL,
         title = paste("Combined-model parameters vs. Km", title_suffix),
         subtitle = paste0("Small points: individual vials; large points and lines: ",
                           "group means.\nShaded: Sihi et al. (2020) posterior ",
                           "range (0.13-1.22 mg/L); dashed line: base Km = ",
                           Km_base, " mg/L")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

plot_errors_vs_Km <- function(x, title_suffix) {
  ggplot(x, aes(Km_mg_L, RMSE, colour = Treatment, shape = Moisture)) +
    annotate("rect", xmin = 0.13, xmax = 1.22, ymin = 0, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_vline(xintercept = Km_base, linetype = "dashed", colour = "grey30") +
    geom_point(aes(group = interaction(Moisture, Treatment)),
               position = position_dodge(width = 0.05), alpha = 0.7, size = 2) +
    stat_summary(aes(group = Moisture, linetype = Moisture), fun = mean,
                 geom = "line", colour = "black") +
    scale_x_log10(breaks = sort(unique(x$Km_mg_L)),
                  labels = function(b) format(b, drop0trailing = TRUE)) +
    scale_colour_manual(values = trt_cols) +
    scale_shape_manual(values = moist_shapes) +
    scale_linetype_manual(values = c("Dry" = "dashed", "Wet" = "solid")) +
    labs(x = "Km (mg L^-1, log scale)", y = "RMSE (mg O2 L^-1)",
         title = paste("Fit error vs. Km", title_suffix),
         subtitle = "Black lines: mean RMSE by inundation history") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (g in c("A", "B")) {
  x <- per_vial %>% filter(fit_ok, if (g == "A") In_GridA else In_GridB)
  lab <- if (g == "A") "(Grid A: 0.05-1.22 mg/L)" else "(Grid B: 0.05-10 mg/L)"
  p1 <- plot_parameters_vs_Km(x, lab)
  p2 <- plot_errors_vs_Km(x, lab)
  ggsave(file.path(out_dir, paste0("Fig_KmSensitivity_Parameters_Grid", g, ".pdf")),
         p1, width = 8, height = 8)
  ggsave(file.path(out_dir, paste0("Fig_KmSensitivity_Parameters_Grid", g, ".png")),
         p1, width = 8, height = 8, dpi = 300)
  ggsave(file.path(out_dir, paste0("Fig_KmSensitivity_RMSE_Grid", g, ".pdf")),
         p2, width = 8, height = 5)
  ggsave(file.path(out_dir, paste0("Fig_KmSensitivity_RMSE_Grid", g, ".png")),
         p2, width = 8, height = 5, dpi = 300)
}

# QC: observed vs fitted DO for every vial at the base Km
if (save_fit_qc_plot) {
  base <- per_vial %>% filter(Km_mg_L == Km_base, fit_ok)
  qc <- bind_rows(lapply(seq_len(nrow(base)), function(i) {
    r <- base[i, ]
    d <- prep[[r$Sample_Name]]$data
    d$Fitted <- solve_model(r$kL_per_h, r$Vmax_mg_L_h, r$DO0_est, Km_base, d$time)
    d$Sample_Name <- r$Sample_Name
    d$Label <- sprintf("%s\nkL=%.2f, Vmax=%.1f, RMSE=%.3f",
                       r$Sample_Name, r$kL_per_h, r$Vmax_mg_L_h, r$RMSE)
    d
  }))
  p_qc <- ggplot(qc, aes(time * 60)) +
    geom_point(aes(y = DO), size = 0.4, colour = "grey50") +
    geom_line(aes(y = Fitted), colour = "red", linewidth = 0.6) +
    facet_wrap(~Label, ncol = 6, scales = "free") +
    labs(x = "Time since start of fitted window (min)", y = "DO (mg/L)",
         title = paste0("Combined model fits at Km = ", Km_base, " mg/L")) +
    theme_bw(base_size = 7)
  ggsave(file.path(out_dir, paste0("QC_ModelFits_Km", Km_base, ".pdf")),
         p_qc, width = 16, height = 11)
}

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
message("\nFits failed: ", sum(!per_vial$fit_ok))
message("Fits where starting guesses disagreed: ",
        sum(per_vial$n_starts_agree < per_vial$n_starts, na.rm = TRUE))
print(compare_dry_wet(per_vial %>% filter(fit_ok)) %>%
        select(Km_mg_L, kL_dry, kL_wet, Vmax_dry, Vmax_wet, Vmax_direction) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))))
message("Outputs written to: ", out_dir)