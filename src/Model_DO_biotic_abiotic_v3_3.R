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
#  13. Diagnostics: DO0_minus_first_obs per fit, and the number of vials whose
#      fitted DO0 differs from the first observation by more than
#      do0_gap_flag, per Km and sediment type.
#  14. Equifinality diagnostics: parameter correlations for every fit; at
#      Km_base, a kL profile (SSR as kL is fixed at a series of values) and
#      moving-block bootstrap confidence intervals for kL, Vmax and DO0.
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
out_dir <- "modeling_outputs_v3"

Km_base  <- 0.2                                                   # mg/L
Km_gridA <- c(0.05, 0.1, 0.13, 0.2, 0.3, 0.51, 0.75, 1.0, 1.22)   # mg/L
Km_gridB <- c(Km_gridA, 3.2, 5.28, 10)                            # mg/L

trim_start_min   <- 2      # minutes removed at start
trim_end_min     <- 2      # minutes removed at end
anoxic_threshold <- 0.1    # mg/L; series cut before first value <= this
do0_gap_flag     <- 0.3    # mg/L; flag |DO0_est - first observation| above this

bounds_lower <- c(kL = 0,  Vmax = 0,   DO0 = 0)
bounds_upper <- c(kL = 50, Vmax = 500, DO0 = 15)

save_fit_qc_plot <- TRUE   # observed vs fitted DO for all vials at Km_base

# Equifinality diagnostics (run at Km_base only)
run_kL_profile  <- TRUE    # refit with kL fixed at a series of values
profile_ssr_tol <- 0.10    # kL range = values with SSR within 10% of the best
run_bootstrap   <- TRUE    # moving-block residual bootstrap
n_boot          <- 200     # bootstrap replicates per vial
boot_acf_cutoff <- 0.2     # block length = first lag where residual ACF < this
boot_seed       <- 42

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
                       S = "Control", U = "Unburned + DOC", H = "High Burn + DOC",
                       .default = NA_character_),
    Sediment_Type = paste(Moisture, "Sediments"),
    Incubation_Set = sub("^YEP[12]([A-Z]).*", "\\1", Sample_Name)
  )

if (any(is.na(df$Moisture)) || any(is.na(df$Treatment))) {
  stop("Unrecognized sample names: ",
       paste(unique(df$Sample_Name[is.na(df$Moisture) | is.na(df$Treatment)]),
             collapse = ", "))
}
df$Treatment <- factor(df$Treatment,
                       levels = c("Control", "Unburned + DOC", "High Burn + DOC"))
df$Sediment_Type <- factor(df$Sediment_Type,
                           levels = c("Dry Sediments", "Wet Sediments"))

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
# Parameter correlations from the Jacobian of the model at the best fit.
# When kL sits at its lower bound (0), only the Vmax-DO0 correlation is
# meaningful, and the kL correlations are returned as NA.
# -----------------------------------------------------------------------------
param_correlations <- function(p, Km, times) {
  out <- c(r_kL_Vmax = NA_real_, r_kL_DO0 = NA_real_, r_Vmax_DO0 = NA_real_)
  nm <- if (p[["kL"]] < 1e-3) c("Vmax", "DO0") else c("kL", "Vmax", "DO0")
  J <- sapply(nm, function(k) {
    h  <- 1e-4 * max(abs(p[[k]]), 0.01)
    up <- p; dn <- p
    up[[k]] <- p[[k]] + h
    dn[[k]] <- max(p[[k]] - h, 0)
    (solve_model(up[["kL"]], up[["Vmax"]], up[["DO0"]], Km, times) -
       solve_model(dn[["kL"]], dn[["Vmax"]], dn[["DO0"]], Km, times)) /
      (up[[k]] - dn[[k]])
  })
  cm <- tryCatch(cov2cor(solve(crossprod(J))), error = function(e) NULL)
  if (is.null(cm) || any(!is.finite(cm))) return(out)
  if (length(nm) == 2) {
    out[["r_Vmax_DO0"]] <- cm[1, 2]
  } else {
    out[["r_kL_Vmax"]]  <- cm[1, 2]
    out[["r_kL_DO0"]]   <- cm[1, 3]
    out[["r_Vmax_DO0"]] <- cm[2, 3]
  }
  out
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
  cors <- param_correlations(p, Km, times)
  n <- length(obs); k <- 3
  SSR <- sum(res^2)
  lag1 <- if (n > 2) cor(res[-1], res[-n]) else NA_real_

  data.frame(
    fit_ok = TRUE,
    kL_per_h = p[["kL"]],
    Vmax_mg_L_h = p[["Vmax"]],
    DO0_est = p[["DO0"]],
    DO_first_obs = obs[1],
    DO0_minus_first_obs = p[["DO0"]] - obs[1],
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
    r_kL_Vmax = cors[["r_kL_Vmax"]],
    r_kL_DO0 = cors[["r_kL_DO0"]],
    r_Vmax_DO0 = cors[["r_Vmax_DO0"]],
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

meta <- df %>% distinct(Sample_Name, Moisture, Sediment_Type, Treatment, Incubation_Set)

prep_log <- bind_rows(lapply(samples, function(s)
  cbind(Sample_Name = s, prep[[s]]$log))) %>%
  left_join(meta, by = "Sample_Name") %>%
  relocate(Moisture, Sediment_Type, Treatment, Incubation_Set, .after = Sample_Name)
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
  mutate(
    In_GridA = Km_mg_L %in% Km_gridA,
    In_GridB = Km_mg_L %in% Km_gridB
  ) %>%
  relocate(Moisture, Sediment_Type, Treatment, Incubation_Set, .after = Sample_Name)

write.csv(per_vial, file.path(out_dir, "YEP_KmSensitivity_PerVial_AllKm.csv"),
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
      r_kL_Vmax_median = median(r_kL_Vmax, na.rm = TRUE),
      r_kL_Vmax_strongest = if (all(is.na(r_kL_Vmax))) NA_real_ else
        r_kL_Vmax[which.max(abs(r_kL_Vmax))],
      DO0_gap_abs_max = max(abs(DO0_minus_first_obs)),
      !!paste0("vials_DO0_gap_gt_", do0_gap_flag) :=
        sum(abs(DO0_minus_first_obs) > do0_gap_flag),
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
# Plots (colors, labels and theme follow the Figure 1 script,
# Plot_fitted_Vmax_kL_DryWet_Wilcox.R)
# -----------------------------------------------------------------------------
sediment_colors <- c(
  "Dry Sediments" = "#56B4E9",
  "Wet Sediments" = "#1B4F72"
)
treatment_shapes    <- c("Control" = 16, "Unburned + DOC" = 17,
                         "High Burn + DOC" = 15)
treatment_linetypes <- c("Control" = "solid", "Unburned + DOC" = "dashed",
                         "High Burn + DOC" = "dotted")
parameter_labels <- c(
  "Vmax" = "Biotic O\u2082 consumption at saturation,\nVmax (mg O\u2082 L\u207b\u00b9 h\u207b\u00b9)",
  "kL"   = "Chemical reaction rate constant,\nkL (h\u207b\u00b9)"
)

pub_theme <- theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey85", color = "grey20"),
    strip.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 11, hjust = 0),
    plot.margin = margin(10, 24, 10, 10)
  )

km_axis <- function(x) {
  scale_x_log10(breaks = sort(unique(x$Km_mg_L)),
                labels = function(b) format(b, drop0trailing = TRUE))
}

km_caption <- paste0(
  "Small points: individual microcosms; large points and lines: group means. ",
  "Shaded: Sihi et al. (2020)\nposterior range (0.13-1.22 mg/L); ",
  "dashed line: Km = ", Km_base, " mg/L used in the main analysis.")

plot_parameters_vs_Km <- function(x) {
  long <- x %>%
    select(Sample_Name, Sediment_Type, Treatment, Km_mg_L,
           kL_per_h, Vmax_mg_L_h) %>%
    pivot_longer(c(Vmax_mg_L_h, kL_per_h), names_to = "Parameter",
                 values_to = "Value") %>%
    mutate(Parameter = factor(recode(Parameter, Vmax_mg_L_h = "Vmax",
                                     kL_per_h = "kL"),
                              levels = c("Vmax", "kL")))
  means <- long %>%
    group_by(Parameter, Km_mg_L, Sediment_Type, Treatment) %>%
    summarise(Value = mean(Value), .groups = "drop")

  ggplot(long, aes(Km_mg_L, Value, colour = Sediment_Type, shape = Treatment)) +
    annotate("rect", xmin = 0.13, xmax = 1.22, ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = Km_base, linetype = "dashed", colour = "grey30") +
    geom_point(aes(group = interaction(Sediment_Type, Treatment)),
               position = position_dodge(width = 0.05), alpha = 0.35, size = 1.8) +
    geom_line(data = means,
              aes(group = interaction(Sediment_Type, Treatment),
                  linetype = Treatment),
              linewidth = 0.8) +
    geom_point(data = means, size = 3.2) +
    km_axis(x) +
    scale_colour_manual(values = sediment_colors) +
    scale_shape_manual(values = treatment_shapes) +
    scale_linetype_manual(values = treatment_linetypes) +
    facet_wrap(~Parameter, ncol = 1, scales = "free_y",
               labeller = labeller(Parameter = parameter_labels)) +
    labs(x = expression(K[m] ~ "(mg" ~ L^{-1} * ", log scale)"), y = NULL,
         colour = "Sediment Type", shape = "Treatment", linetype = "Treatment",
         caption = km_caption) +
    pub_theme
}

plot_errors_vs_Km <- function(x) {
  means <- x %>%
    group_by(Km_mg_L, Sediment_Type) %>%
    summarise(RMSE = mean(RMSE), .groups = "drop")

  ggplot(x, aes(Km_mg_L, RMSE, colour = Sediment_Type)) +
    annotate("rect", xmin = 0.13, xmax = 1.22, ymin = 0, ymax = Inf,
             fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = Km_base, linetype = "dashed", colour = "grey30") +
    geom_point(aes(shape = Treatment,
                   group = interaction(Sediment_Type, Treatment)),
               position = position_dodge(width = 0.05), alpha = 0.6, size = 2.2) +
    geom_line(data = means, aes(group = Sediment_Type), linewidth = 1) +
    km_axis(x) +
    scale_colour_manual(values = sediment_colors) +
    scale_shape_manual(values = treatment_shapes) +
    labs(x = expression(K[m] ~ "(mg" ~ L^{-1} * ", log scale)"),
         y = expression("RMSE (mg O"[2] ~ L^{-1} * ")"),
         colour = "Sediment Type", shape = "Treatment",
         caption = paste0("Points: individual microcosms; lines: mean RMSE by ",
                          "sediment type. Shaded and dashed line as in the ",
                          "parameter figure.")) +
    pub_theme
}

save_plot <- function(p, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
  tryCatch(
    ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
           width = width, height = height, device = cairo_pdf),
    error = function(e) message("PDF not saved for ", name, ": ", conditionMessage(e)))
}

for (g in c("A", "B")) {
  x <- per_vial %>% filter(fit_ok, if (g == "A") In_GridA else In_GridB)
  save_plot(plot_parameters_vs_Km(x),
            paste0("Fig_KmSensitivity_Parameters_Grid", g), 10, 11)
  save_plot(plot_errors_vs_Km(x),
            paste0("Fig_KmSensitivity_RMSE_Grid", g), 10, 7)
}

# -----------------------------------------------------------------------------
# Equifinality diagnostics at the base Km
# -----------------------------------------------------------------------------
base_fits <- per_vial %>% filter(Km_mg_L == Km_base, fit_ok)

# Least-squares fit from one starting point, with optional fixed kL.
# Used by the profile (kL fixed) and the bootstrap (all three free).
fit_from <- function(times, obs, Km, start, kL_fixed = NULL) {
  free <- if (is.null(kL_fixed)) c("kL", "Vmax", "DO0") else c("Vmax", "DO0")
  rf <- function(q) {
    kL <- if (is.null(kL_fixed)) q[["kL"]] else kL_fixed
    r <- obs - solve_model(kL, q[["Vmax"]], q[["DO0"]], Km, times)
    r[!is.finite(r)] <- 1e3
    r
  }
  lo <- bounds_lower[free]; up <- bounds_upper[free]
  run <- function(s) tryCatch(
    nls.lm(par = pmin(pmax(s, lo + 1e-6), up - 1e-6), fn = rf,
           lower = lo, upper = up,
           control = nls.lm.control(maxiter = 300, ftol = 1e-12, ptol = 1e-12)),
    error = function(e) NULL)
  f <- run(start[free])
  if (is.null(f)) return(NULL)
  f2 <- run(f$par)                                   # polish
  if (!is.null(f2) && sum(f2$fvec^2) <= sum(f$fvec^2)) f <- f2
  list(par = f$par, SSR = sum(f$fvec^2))
}

# ---- kL profile -------------------------------------------------------------
# For each vial, kL is fixed at a series of values around its best estimate,
# Vmax and DO0 are refit, and SSR is recorded. A sharp V-shaped profile means
# kL is well determined; a flat profile means it is not.
profile_grid <- function(kL_hat) {
  g <- if (kL_hat < 0.5) {         # near zero: fixed absolute grid
    c(0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1)
  } else {
    kL_hat * (1 + c(-1, -0.6, -0.4, -0.25, -0.15, -0.08, 0,
                    0.08, 0.15, 0.25, 0.4, 0.6, 1))
  }
  sort(unique(round(pmax(c(g, kL_hat), 0), 6)))
}

# kL values whose SSR stays within (1 + tol) x best SSR, by interpolation
profile_range <- function(kL, rel, kL_hat, tol) {
  thr <- 1 + tol
  i0 <- which.min(abs(kL - kL_hat))
  left <- seq_len(i0); right <- i0:length(kL)
  find_cross <- function(idx, rev = FALSE) {
    if (rev) idx <- rev(idx)
    above <- which(rel[idx] > thr)
    if (length(above) == 0) return(list(value = kL[idx[length(idx)]], open = TRUE))
    j <- above[1]
    if (j == 1) return(list(value = kL[idx[1]], open = FALSE))
    a <- idx[j - 1]; b <- idx[j]
    list(value = approx(rel[c(a, b)], kL[c(a, b)], xout = thr)$y, open = FALSE)
  }
  lw <- find_cross(left, rev = TRUE)
  rt <- find_cross(right)
  data.frame(kL_profile_lower = lw$value, kL_profile_upper = rt$value,
             kL_profile_lower_open = lw$open && kL[1] > 0,
             kL_profile_upper_open = rt$open)
}

if (run_kL_profile) {
  message("\nkL profile at Km = ", Km_base)
  prof <- bind_rows(lapply(seq_len(nrow(base_fits)), function(i) {
    r <- base_fits[i, ]
    d <- prep[[r$Sample_Name]]$data
    start <- c(kL = r$kL_per_h, Vmax = r$Vmax_mg_L_h, DO0 = r$DO0_est)
    bind_rows(lapply(profile_grid(r$kL_per_h), function(k) {
      f <- fit_from(d$time, d$DO, Km_base, start, kL_fixed = k)
      data.frame(Sample_Name = r$Sample_Name, kL_fixed = k,
                 Vmax_refit = if (is.null(f)) NA else f$par[["Vmax"]],
                 DO0_refit  = if (is.null(f)) NA else f$par[["DO0"]],
                 SSR = if (is.null(f)) NA else f$SSR)
    }))
  })) %>%
    left_join(base_fits %>% select(Sample_Name, Sediment_Type, Treatment,
                                   kL_best = kL_per_h, SSR_best = SSR),
              by = "Sample_Name") %>%
    mutate(SSR_relative = SSR / pmin(SSR_best, ave(SSR, Sample_Name,
                                                   FUN = function(z) min(z, na.rm = TRUE))))

  prof_range <- prof %>%
    filter(is.finite(SSR_relative)) %>%
    group_by(Sample_Name) %>%
    group_modify(~ profile_range(.x$kL_fixed, .x$SSR_relative,
                                 .x$kL_best[1], profile_ssr_tol)) %>%
    ungroup()

  write.csv(prof, file.path(out_dir, paste0("YEP_kL_Profile_Km", Km_base, ".csv")),
            row.names = FALSE)

  p_prof <- ggplot(prof, aes(kL_fixed, SSR_relative, colour = Sediment_Type)) +
    geom_hline(yintercept = 1 + profile_ssr_tol, linetype = "dashed",
               colour = "grey40") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    geom_vline(aes(xintercept = kL_best), colour = "grey30", linewidth = 0.3) +
    scale_colour_manual(values = sediment_colors) +
    scale_y_log10() +
    facet_wrap(~Sample_Name, ncol = 6, scales = "free_x") +
    labs(x = expression(k[L] ~ "fixed at (" * h^{-1} * ")"),
         y = "SSR relative to best fit (log scale)",
         colour = "Sediment Type",
         caption = paste0("Vmax and initial DO refit at each fixed kL (Km = ",
                          Km_base, " mg/L). Vertical line: best-fit kL. ",
                          "Dashed line: SSR ", round(100 * profile_ssr_tol),
                          "% above the best fit.")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8),
          plot.caption = element_text(hjust = 0))
  save_plot(p_prof, paste0("Fig_kL_Profile_Km", Km_base), 14, 11)
}

# ---- Moving-block residual bootstrap -----------------------------------------
# Residuals from the best fit are resampled in contiguous blocks, so the
# autocorrelation of the 5-s readings is kept. Block length for each vial is
# the first lag at which its residual autocorrelation falls below
# boot_acf_cutoff (at most a quarter of the series). Each bootstrap series
# (fitted curve + resampled residuals) is refit, and the 2.5th and 97.5th
# percentiles give 95% intervals.
if (run_bootstrap) {
  message("Bootstrap at Km = ", Km_base, " (", n_boot, " replicates per vial)")
  set.seed(boot_seed)
  boot_draws <- bind_rows(lapply(seq_len(nrow(base_fits)), function(i) {
    r <- base_fits[i, ]
    d <- prep[[r$Sample_Name]]$data
    p_hat <- c(kL = r$kL_per_h, Vmax = r$Vmax_mg_L_h, DO0 = r$DO0_est)
    pred <- solve_model(p_hat[["kL"]], p_hat[["Vmax"]], p_hat[["DO0"]],
                        Km_base, d$time)
    res <- d$DO - pred
    res <- res - mean(res)
    n <- length(res)
    max_L <- max(floor(n / 4), 2)
    a <- acf(res, lag.max = max_L, plot = FALSE)$acf[-1]
    L <- which(a < boot_acf_cutoff)[1]
    if (is.na(L)) L <- max_L
    L <- min(max(L, 2), max_L)
    message("  ", r$Sample_Name, ": block length ", L, " readings")
    bind_rows(lapply(seq_len(n_boot), function(b) {
      starts <- sample.int(n - L + 1, ceiling(n / L), replace = TRUE)
      idx <- unlist(lapply(starts, function(s0) s0:(s0 + L - 1)))[seq_len(n)]
      f <- fit_from(d$time, pred + res[idx], Km_base, p_hat)
      data.frame(Sample_Name = r$Sample_Name, replicate = b,
                 block_length = L,
                 kL = if (is.null(f)) NA else f$par[["kL"]],
                 Vmax = if (is.null(f)) NA else f$par[["Vmax"]],
                 DO0 = if (is.null(f)) NA else f$par[["DO0"]])
    }))
  }))
  write.csv(boot_draws,
            file.path(out_dir, paste0("YEP_Bootstrap_Km", Km_base, "_Draws.csv")),
            row.names = FALSE)

  boot_ci <- boot_draws %>%
    group_by(Sample_Name) %>%
    summarise(
      boot_n_ok = sum(!is.na(kL)),
      boot_block_length_readings = first(block_length),
      kL_ci_lower   = quantile(kL, 0.025, na.rm = TRUE),
      kL_ci_upper   = quantile(kL, 0.975, na.rm = TRUE),
      Vmax_ci_lower = quantile(Vmax, 0.025, na.rm = TRUE),
      Vmax_ci_upper = quantile(Vmax, 0.975, na.rm = TRUE),
      DO0_ci_lower  = quantile(DO0, 0.025, na.rm = TRUE),
      DO0_ci_upper  = quantile(DO0, 0.975, na.rm = TRUE),
      r_kL_Vmax_boot = if (sd(kL, na.rm = TRUE) > 0)
        cor(kL, Vmax, use = "complete.obs") else NA_real_,
      .groups = "drop")
}

# ---- Base-Km per-vial file, with diagnostics joined -------------------------
base_out <- base_fits
if (run_kL_profile) base_out <- base_out %>% left_join(prof_range, by = "Sample_Name")
if (run_bootstrap)  base_out <- base_out %>% left_join(boot_ci, by = "Sample_Name")
write.csv(base_out,
          file.path(out_dir, paste0("YEP_CombinedModel_Km", Km_base, "_PerVial.csv")),
          row.names = FALSE)

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
