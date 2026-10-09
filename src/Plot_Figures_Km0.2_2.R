# =============================================================================
# YEP figures at Km = 0.2 mg/L (from Model_DO_biotic_abiotic_v3.R outputs)
# -----------------------------------------------------------------------------
#   Figure 1  : Vmax and kL, dry vs wet sediments (Wilcoxon rank-sum)
#   Figure 2  : Vmax and kL by treatment within each sediment type
#               (Kruskal-Wallis; Dunn tests with Holm adjustment if KW p < 0.05)
#   Message 1 : Wet sediments consumed O2 faster (measured rates) and had
#               higher abiotic rate constants (kL)
#   Message 2 : Dry sediments: kL ~ 0 (almost all biotic) and higher Vmax.
#               This is Figure 1, so no separate figure is made.
#   Message 3 : Treatment effects on measured rates, Vmax and kL
#
# Total O2 consumption uses the measured rates from the data package
# (segmented regression; YEP_INC_20261007.csv). Biotic vs abiotic is shown
# with the fitted Vmax (biotic) and kL (abiotic) directly.
#
# With normalize_by_mass = TRUE, rates, Vmax and kL are expressed per kg of
# dry sediment (multiplied by litres of water / kg dry sediment, where water
# includes pore water). Set it to FALSE for the per-litre versions (SI).
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
model_output_path <- file.path("modeling_outputs_v3",
                               "YEP_CombinedModel_Km0.2_PerVial.csv")
normalize_by_mass <- TRUE   # TRUE: per kg dry sediment; FALSE: per litre of water
Km      <- 0.2   # mg/L; must match the fitted file
rates_path <- "v2_data/YEP_INC_20261007.csv"   # measured rates (data package)
mass_path  <- file.path("v2_data", "v2_YEP_Sample_Data",
                        "v2_YEP_Sediment_Water_Mass_Volume.csv")
fig_dir <- if (normalize_by_mass) "Figures_Km0.2_per_kg" else "Figures_Km0.2_per_L"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

sediment_order  <- c("Dry Sediments", "Wet Sediments")
treatment_order <- c("Control", "Unburned + DOC", "High Burn + DOC")

# Colors from Figure 1 and Figure 2 scripts
sediment_colors <- c("Dry Sediments" = "#56B4E9", "Wet Sediments" = "#1B4F72")
treatment_colors <- c("Control" = "#0072B2", "Unburned + DOC" = "darkgreen",
                      "High Burn + DOC" = "#8B4513")

# Mass normalization: value x (litres of water / kg dry sediment).
# Water includes pore water (Water_Mass_g = added water + pore water).
# Vmax and measured rates: mg O2 L-1 h-1 -> mg O2 kg-1 h-1
# kL:                      h-1           -> L kg-1 h-1
variable_labels <- if (normalize_by_mass) c(
  "Vmax"     = "Biotic O\u2082 consumption at saturation,\nVmax (mg O\u2082 kg\u207b\u00b9 h\u207b\u00b9)",
  "kL"       = "Chemical reaction rate constant,\nkL (L kg\u207b\u00b9 h\u207b\u00b9)",
  "Measured" = "Measured O\u2082 consumption rate\n(mg O\u2082 kg\u207b\u00b9 h\u207b\u00b9)"
) else c(
  "Vmax"     = "Biotic O\u2082 consumption at saturation,\nVmax (mg O\u2082 L\u207b\u00b9 h\u207b\u00b9)",
  "kL"       = "Chemical reaction rate constant,\nkL (h\u207b\u00b9)",
  "Measured" = "Measured O\u2082 consumption rate\n(mg O\u2082 L\u207b\u00b9 h\u207b\u00b9)"
)
basis_note <- if (normalize_by_mass) "Values per kg dry sediment." else "Values per litre of water."

pub_theme <- theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black"),
    strip.background = element_rect(fill = "grey85", color = "grey20"),
    strip.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 11, hjust = 0),
    plot.margin = margin(10, 24, 10, 10)
  )

save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
  tryCatch(
    ggsave(file.path(fig_dir, paste0(name, ".pdf")), p,
           width = width, height = height, device = cairo_pdf),
    error = function(e) message("PDF not saved for ", name))
}

p_label <- function(p, prefix) {
  case_when(is.na(p) ~ paste(prefix, "p = NA"),
            p < 0.001 ~ paste(prefix, "p < 0.001"),
            TRUE ~ sprintf("%s p = %.3f", prefix, p))
}

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
measured <- read_csv(rates_path, show_col_types = FALSE) %>%
  filter(!str_detect(Sample_Name, "^YEP_INC-W")) %>%      # water-only controls
  transmute(Sample_Name,
            Measured_Rate = abs(Respiration_Rate_mg_DO_per_L_per_H))

mass <- read_csv(mass_path, skip = 2, show_col_types = FALSE) %>%
  filter(str_detect(Sample_Name, "_INC-")) %>%
  transmute(Sample_Name,
            Dry_Sediment_Mass_g = as.numeric(Dry_Sediment_Mass_g),
            Water_Mass_g = as.numeric(Water_Mass_g),
            L_water_per_kg = (Water_Mass_g / 1000) / (Dry_Sediment_Mass_g / 1000))

fits <- read_csv(model_output_path, show_col_types = FALSE) %>%
  filter(fit_ok) %>%
  left_join(measured, by = "Sample_Name") %>%
  left_join(mass, by = "Sample_Name") %>%
  mutate(
    Sediment_Type = factor(Sediment_Type, levels = sediment_order),
    Treatment = factor(Treatment, levels = treatment_order),
    norm = if (normalize_by_mass) L_water_per_kg else 1,
    Vmax_plot     = Vmax_mg_L_h * norm,
    kL_plot       = kL_per_h * norm,
    Measured_plot = Measured_Rate * norm
  )

if (any(abs(fits$Km_mg_L - Km) > 1e-9)) {
  stop("Km in the fitted file does not match the Km setting in this script.")
}
if (any(is.na(fits$Measured_Rate)) || nrow(fits) != nrow(measured)) {
  stop("Sample names in the rates file and the model output do not match.")
}
if (normalize_by_mass && any(is.na(fits$L_water_per_kg))) {
  stop("Some vials are missing from the mass/volume file.")
}
write_csv(fits %>% select(Sample_Name, Sediment_Type, Treatment,
                          Dry_Sediment_Mass_g, Water_Mass_g, L_water_per_kg,
                          Vmax_plot, kL_plot, Measured_plot),
          file.path(fig_dir, "Plotted_values.csv"))

to_long <- function(d, vars) {
  d %>%
    select(Sample_Name, Sediment_Type, Treatment, all_of(unname(vars))) %>%
    pivot_longer(all_of(unname(vars)), names_to = "Variable",
                 values_to = "Value") %>%
    mutate(Variable = factor(names(vars)[match(Variable, vars)],
                             levels = names(vars)))
}

# -----------------------------------------------------------------------------
# Dry vs wet comparison (Figure 1, Message 2)
# -----------------------------------------------------------------------------
dry_wet_stats <- function(long) {
  long %>%
    group_by(Variable) %>%
    summarise(
      mean_dry = mean(Value[Sediment_Type == "Dry Sediments"]),
      mean_wet = mean(Value[Sediment_Type == "Wet Sediments"]),
      p = wilcox.test(Value ~ Sediment_Type, exact = FALSE)$p.value,
      y_max = max(Value), y_min = min(Value),
      .groups = "drop") %>%
    mutate(y_rng = pmax(y_max - y_min, abs(y_max), 1e-6),
           label = p_label(p, "Wilcox"),
           x = 0.6, y = y_max + 0.2 * y_rng)
}

plot_dry_wet <- function(long, stats) {
  ggplot(long, aes(Sediment_Type, Value, fill = Sediment_Type)) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.12, seed = 1),
               alpha = 0.75, size = 2) +
    geom_blank(data = stats, aes(x = 1, y = y_max + 0.3 * y_rng),
               inherit.aes = FALSE) +
    geom_text(data = stats, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 4.2, fontface = "bold") +
    scale_fill_manual(values = sediment_colors) +
    scale_x_discrete(labels = c("Dry Sediments" = "Dry\nSediments",
                                "Wet Sediments" = "Wet\nSediments")) +
    labs(x = NULL, y = NULL, fill = "Sediment Type") +
    coord_cartesian(clip = "off") +
    facet_wrap(~Variable, nrow = 1, scales = "free_y",
               labeller = labeller(Variable = variable_labels)) +
    pub_theme
}

# -----------------------------------------------------------------------------
# Treatment comparison within sediment type (Figure 2, Message 3)
# -----------------------------------------------------------------------------
treatment_stats <- function(long) {
  kw <- long %>%
    group_by(Variable, Sediment_Type) %>%
    kruskal_test(Value ~ Treatment) %>%
    ungroup() %>%
    mutate(label = p_label(p, "KW"), kw_sig = p < 0.05)

  dunn_in <- long %>%
    semi_join(kw %>% filter(kw_sig), by = c("Variable", "Sediment_Type"))
  dunn <- if (nrow(dunn_in) > 0) {
    dunn_in %>%
      group_by(Variable, Sediment_Type) %>%
      dunn_test(Value ~ Treatment, p.adjust.method = "holm") %>%
      ungroup()
  } else {
    tibble(Variable = factor(), Sediment_Type = factor(), group1 = character(),
           group2 = character(), p = numeric(), p.adj = numeric())
  }

  # Panels in a row share the y-axis (facet_grid), so label positions use
  # the range of each variable across both sediment types
  ranges <- long %>%
    group_by(Variable) %>%
    summarise(y_max = max(Value), y_min = min(Value), .groups = "drop") %>%
    mutate(y_rng = pmax(y_max - y_min, abs(y_max), 1e-6)) %>%
    tidyr::crossing(Sediment_Type = factor(sediment_order, levels = sediment_order))

  kw_lab <- kw %>% left_join(ranges, by = c("Variable", "Sediment_Type")) %>%
    mutate(x = 0.55, y = y_max + 0.2 * y_rng)

  # Brackets only for pairwise comparisons with adjusted p < 0.05
  brackets <- dunn %>%
    filter(p.adj < 0.05) %>%
    mutate(x1 = match(group1, treatment_order),
           x2 = match(group2, treatment_order),
           xmid = (x1 + x2) / 2, width = abs(x2 - x1)) %>%
    left_join(ranges, by = c("Variable", "Sediment_Type")) %>%
    arrange(Variable, Sediment_Type, width) %>%
    group_by(Variable, Sediment_Type) %>%
    mutate(y_pos = y_max + (0.42 + 0.13 * (row_number() - 1)) * y_rng,
           y_tip = y_pos - 0.04 * y_rng,
           label = p_label(p.adj, "Dunn")) %>%
    ungroup()

  tops <- ranges %>% mutate(top = y_max + 0.35 * y_rng)
  if (nrow(brackets) > 0) {
    tops <- tops %>%
      left_join(brackets %>% group_by(Variable, Sediment_Type) %>%
                  summarise(top_br = max(y_pos + 0.12 * y_rng), .groups = "drop"),
                by = c("Variable", "Sediment_Type")) %>%
      mutate(top = pmax(top, top_br, na.rm = TRUE))
  }
  list(kw = kw, dunn = dunn, kw_lab = kw_lab, brackets = brackets, tops = tops)
}

plot_treatments <- function(long, st) {
  ggplot(long, aes(Treatment, Value, fill = Treatment)) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.12, seed = 1),
               alpha = 0.75, size = 2) +
    geom_blank(data = st$tops, aes(x = 1, y = top), inherit.aes = FALSE) +
    geom_text(data = st$kw_lab, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 4.2, fontface = "bold") +
    geom_segment(data = st$brackets, aes(x = x1, xend = x2, y = y_pos, yend = y_pos),
                 inherit.aes = FALSE, linewidth = 0.7) +
    geom_segment(data = st$brackets, aes(x = x1, xend = x1, y = y_tip, yend = y_pos),
                 inherit.aes = FALSE, linewidth = 0.7) +
    geom_segment(data = st$brackets, aes(x = x2, xend = x2, y = y_tip, yend = y_pos),
                 inherit.aes = FALSE, linewidth = 0.7) +
    geom_text(data = st$brackets, aes(x = xmid, y = y_pos, label = label),
              inherit.aes = FALSE, vjust = -0.35, size = 4.2, fontface = "bold") +
    scale_fill_manual(values = treatment_colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
    labs(x = NULL, y = NULL, fill = "Treatment") +
    coord_cartesian(clip = "off") +
    facet_grid(Variable ~ Sediment_Type, scales = "free_y",
               labeller = labeller(Variable = variable_labels)) +
    pub_theme +
    theme(strip.text.y = element_text(size = 11))
}

write_treatment_stats <- function(st, name) {
  write_csv(st$kw, file.path(fig_dir, paste0(name, "_KW_stats.csv")))
  write_csv(st$dunn, file.path(fig_dir, paste0(name, "_Dunn_stats.csv")))
}

# -----------------------------------------------------------------------------
# Figure 1: Vmax and kL, dry vs wet
# -----------------------------------------------------------------------------
fig1_long  <- to_long(fits, c(Vmax = "Vmax_plot", kL = "kL_plot"))
fig1_stats <- dry_wet_stats(fig1_long)
write_csv(fig1_stats %>% select(Variable, mean_dry, mean_wet, p),
          file.path(fig_dir, "Fig1_Vmax_kL_DryWet_Wilcox_stats.csv"))
save_fig(plot_dry_wet(fig1_long, fig1_stats) +
           labs(caption = paste0("Combined model fits with Km = ", Km, " mg/L. ", basis_note)),
         "Fig1_Vmax_kL_DryWet_Wilcox", 9, 5.5)

# -----------------------------------------------------------------------------
# Figure 2: Vmax and kL by treatment within sediment type
# -----------------------------------------------------------------------------
fig2_long  <- fig1_long
fig2_stats <- treatment_stats(fig2_long)
write_treatment_stats(fig2_stats, "Fig2_Vmax_kL_Treatment")
save_fig(plot_treatments(fig2_long, fig2_stats) +
           labs(caption = paste0("Kruskal-Wallis test within each sediment type; ",
                                 "Dunn tests (Holm-adjusted) shown where KW p < 0.05 ",
                                 "and pairwise p < 0.05. Km = ", Km, " mg/L. ", basis_note)),
         "Fig2_Vmax_kL_Treatment_KW_Dunn", 14, 10)

# -----------------------------------------------------------------------------
# Message 1: wet sediments consumed O2 faster and had higher abiotic rate
# constants. Measured rates and kL, dry vs wet.
# (Message 2 -- kL ~ 0 and higher Vmax in dry sediments -- is Figure 1.)
# -----------------------------------------------------------------------------
msg1_long  <- to_long(fits, c(Measured = "Measured_plot", kL = "kL_plot"))
msg1_stats <- dry_wet_stats(msg1_long)
write_csv(msg1_stats %>% select(Variable, mean_dry, mean_wet, p),
          file.path(fig_dir, "Msg1_MeasuredRate_kL_DryWet_stats.csv"))
save_fig(plot_dry_wet(msg1_long, msg1_stats) +
           labs(caption = paste0("Measured rates: segmented regression (data package). ",
                                 "kL: combined model, Km = ", Km, " mg/L.\n", basis_note)),
         "Msg1_Wet_Faster_Higher_Abiotic", 9, 5.5)

# -----------------------------------------------------------------------------
# Message 3: treatment effects on measured rates, Vmax and kL
# -----------------------------------------------------------------------------
msg3_long  <- to_long(fits, c(Measured = "Measured_plot", Vmax = "Vmax_plot",
                              kL = "kL_plot"))
msg3_stats <- treatment_stats(msg3_long)
write_treatment_stats(msg3_stats, "Msg3_Treatment")
save_fig(plot_treatments(msg3_long, msg3_stats) +
           labs(caption = paste0("Kruskal-Wallis test within each sediment type; ",
                                 "Dunn tests (Holm-adjusted) shown where KW p < 0.05 ",
                                 "and pairwise p < 0.05. Km = ", Km, " mg/L. ", basis_note)),
         "Msg3_Treatment_Effects", 14, 13)

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
message("\nFigure 1 (dry vs wet):")
print(fig1_stats %>% select(Variable, mean_dry, mean_wet, p))
message("\nMessage 1 (measured rates and kL, dry vs wet):")
print(msg1_stats %>% select(Variable, mean_dry, mean_wet, p))
message("\nTreatment tests (Kruskal-Wallis):")
print(msg3_stats$kw %>% select(Variable, Sediment_Type, statistic, p))
message("\nDunn tests (where KW p < 0.05):")
print(msg3_stats$dunn %>% select(Variable, Sediment_Type, group1, group2, p.adj))
message("\nFigures written to: ", fig_dir)

