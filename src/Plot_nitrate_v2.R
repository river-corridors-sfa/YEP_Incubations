rm(list = ls())

library(tidyverse)
library(rstatix)

# Run from repo root or from src/.
if (!dir.exists("v2_data") && dir.exists("../v2_data")) {
  setwd("..")
}

fig_dir <- file.path("Figures", "nitrate_v2")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

moisture_colors <- c(
  "Dry" = "#56B4E9",
  "Wet" = "#1B4F72"
)

treatment_order <- names(treatment_colors)
moisture_order <- names(moisture_colors)

ions_path <- file.path("v2_data", "v2_YEP_Sample_Data", "YEP_Sediment_Ions.csv")
sample_metadata_path <- file.path("v2_data", "v2_YEP_Sample_Name_Metadata.csv")

sample_metadata <- read_csv(
  sample_metadata_path,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    Sample_Material = str_squish(Sample_Material),
    Treatment = case_when(
      str_squish(Treatment) == "Columbia River synthetic river water" ~ "Control",
      str_squish(Treatment) == "Unburned Douglas fir leachate" ~ "Unburned + DOC",
      str_squish(Treatment) == "Burned (high severity) Douglas fir leachate" ~ "High Burn + DOC",
      TRUE ~ NA_character_
    ),
    Moisture = case_when(
      str_to_lower(str_squish(Sediment_Field_Moisture_Conditions)) == "dry" ~ "Dry",
      str_to_lower(str_squish(Sediment_Field_Moisture_Conditions)) == "wet" ~ "Wet",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Sample_Name, Sample_Material, Treatment, Moisture, Incubation_Day, Replicate_Identifier)

nitrate_data <- read_csv(
  ions_path,
  skip = 2,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE,
  name_repair = "minimal",
  trim_ws = TRUE
) %>%
  filter(str_detect(Sample_Name, "^YEP[12].*_SIN-[HSU]")) %>%
  select(Sample_Name, `00618_NO3_mg_per_L_as_N`) %>%
  mutate(
    Nitrate_mg_per_L_as_N = parse_number(as.character(`00618_NO3_mg_per_L_as_N`)),
    Nitrate_mg_per_L_as_N = na_if(Nitrate_mg_per_L_as_N, -9999)
  ) %>%
  select(Sample_Name, Nitrate_mg_per_L_as_N) %>%
  left_join(sample_metadata, by = "Sample_Name") %>%
  filter(
    Sample_Material == "Aqueous sample post-incubation (treatment solution and sediment)",
    !is.na(Treatment),
    !is.na(Moisture),
    !is.na(Nitrate_mg_per_L_as_N)
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Moisture = factor(Moisture, levels = moisture_order)
  )

write_csv(
  nitrate_data,
  file.path(fig_dir, "nitrate_plot_data.csv")
)

print(count(nitrate_data, Treatment, Moisture))

wilcox_stats <- nitrate_data %>%
  group_by(Treatment) %>%
  wilcox_test(Nitrate_mg_per_L_as_N ~ Moisture, exact = FALSE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  ungroup() %>%
  mutate(
    p_label = case_when(
      is.na(p.adj) ~ "Wilcox p = NA",
      p.adj < 0.001 ~ "Wilcox p < 0.001",
      TRUE ~ sprintf("Wilcox p = %.3f", p.adj)
    )
  )

write_csv(
  wilcox_stats,
  file.path(fig_dir, "Nitrate_DryWet_Wilcox_stats.csv")
)

kw_stats <- nitrate_data %>%
  group_by(Moisture) %>%
  kruskal_test(Nitrate_mg_per_L_as_N ~ Treatment) %>%
  ungroup() %>%
  mutate(
    kw_label = case_when(
      is.na(p) ~ "KW p = NA",
      p < 0.001 ~ "KW p < 0.001",
      TRUE ~ sprintf("KW p = %.3f", p)
    ),
    kw_sig = p < 0.05
  )

write_csv(
  kw_stats,
  file.path(fig_dir, "Nitrate_Treatment_KW_stats.csv")
)

dunn_input <- nitrate_data %>%
  semi_join(
    kw_stats %>%
      filter(kw_sig) %>%
      select(Moisture),
    by = "Moisture"
  )

if (nrow(dunn_input) > 0) {
  dunn_stats <- dunn_input %>%
    group_by(Moisture) %>%
    dunn_test(Nitrate_mg_per_L_as_N ~ Treatment, p.adjust.method = "holm") %>%
    add_significance("p.adj") %>%
    ungroup()
} else {
  dunn_stats <- tibble(
    Moisture = factor(levels = moisture_order),
    group1 = character(),
    group2 = character(),
    n1 = integer(),
    n2 = integer(),
    statistic = numeric(),
    p = numeric(),
    p.adj = numeric(),
    p.adj.signif = character()
  )
}

write_csv(
  dunn_stats,
  file.path(fig_dir, "Nitrate_Treatment_Dunn_stats.csv")
)

dry_wet_ranges <- nitrate_data %>%
  group_by(Treatment) %>%
  summarise(
    y_min = min(Nitrate_mg_per_L_as_N, na.rm = TRUE),
    y_max = max(Nitrate_mg_per_L_as_N, na.rm = TRUE),
    y_rng = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(y_rng = if_else(y_rng == 0, pmax(abs(y_max), 1), y_rng))

wilcox_labels <- wilcox_stats %>%
  left_join(dry_wet_ranges, by = "Treatment") %>%
  mutate(
    x = 0.55,
    y = y_max + 0.20 * y_rng,
    y_top = y_max + 0.30 * y_rng
  )

treatment_ranges <- nitrate_data %>%
  group_by(Moisture) %>%
  summarise(
    y_min = min(Nitrate_mg_per_L_as_N, na.rm = TRUE),
    y_max = max(Nitrate_mg_per_L_as_N, na.rm = TRUE),
    y_rng = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(y_rng = if_else(y_rng == 0, pmax(abs(y_max), 1), y_rng))

kw_labels <- kw_stats %>%
  left_join(treatment_ranges, by = "Moisture") %>%
  mutate(
    x = 0.55,
    y = y_max + 0.18 * y_rng,
    y_top = y_max + 0.28 * y_rng
  )

dunn_brackets <- tibble()

if (nrow(dunn_stats) > 0) {
  dunn_brackets <- dunn_stats %>%
    left_join(treatment_ranges, by = "Moisture") %>%
    group_by(Moisture) %>%
    mutate(
      x1 = match(group1, treatment_order),
      x2 = match(group2, treatment_order),
      xmid = (x1 + x2) / 2,
      y.position = y_max + (0.36 + 0.16 * row_number()) * y_rng,
      y.tip = y.position - 0.05 * y_rng
    ) %>%
    ungroup()
}

treatment_space <- treatment_ranges %>%
  transmute(Moisture, y_top = y_max + 0.32 * y_rng)

if (nrow(dunn_brackets) > 0) {
  dunn_space <- dunn_brackets %>%
    group_by(Moisture) %>%
    summarise(y_top = max(y.position + 0.08 * y_rng), .groups = "drop")

  treatment_space <- treatment_space %>%
    full_join(dunn_space, by = "Moisture", suffix = c("", "_dunn")) %>%
    mutate(y_top = pmax(y_top, y_top_dunn, na.rm = TRUE)) %>%
    select(Moisture, y_top)
}

pub_theme <- theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black"),
    strip.background = element_rect(fill = "grey85", color = "grey20"),
    strip.text = element_text(size = 16),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 24, 10, 10)
  )

fig_dry_wet <- ggplot(
  nitrate_data,
  aes(x = Moisture, y = Nitrate_mg_per_L_as_N, fill = Moisture)
) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.75,
    size = 2
  ) +
  scale_fill_manual(values = moisture_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = expression(NO[3]^"-" * "-N (mg L"^{-1} * ")"),
    fill = "Moisture"
  ) +
  geom_blank(
    data = wilcox_labels,
    aes(x = 1, y = y_top),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = wilcox_labels,
    aes(x = x, y = y, label = p_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.2,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(~Treatment, nrow = 1, scales = "free_y") +
  pub_theme

fig_treatment <- ggplot(
  nitrate_data,
  aes(x = Treatment, y = Nitrate_mg_per_L_as_N, fill = Treatment)
) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.75,
    size = 2
  ) +
  scale_fill_manual(values = treatment_colors, drop = FALSE) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12), drop = FALSE) +
  labs(
    x = NULL,
    y = expression(NO[3]^"-" * "-N (mg L"^{-1} * ")"),
    fill = "Treatment"
  ) +
  geom_blank(
    data = treatment_space,
    aes(x = 1, y = y_top),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = kw_labels,
    aes(x = x, y = y, label = kw_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.2,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(~Moisture, nrow = 1, scales = "free_y") +
  pub_theme

if (nrow(dunn_brackets) > 0) {
  fig_treatment <- fig_treatment +
    geom_segment(
      data = dunn_brackets,
      aes(x = x1, xend = x2, y = y.position, yend = y.position),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = dunn_brackets,
      aes(x = x1, xend = x1, y = y.tip, yend = y.position),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = dunn_brackets,
      aes(x = x2, xend = x2, y = y.tip, yend = y.position),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_text(
      data = dunn_brackets,
      aes(x = xmid, y = y.position, label = p.adj.signif),
      inherit.aes = FALSE,
      vjust = -0.3,
      size = 4.2,
      fontface = "bold"
    )
}

print(fig_dry_wet)
print(fig_treatment)

ggsave(
  file.path(fig_dir, "Nitrate_DryWet_by_Treatment.png"),
  fig_dry_wet,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "Nitrate_Treatment_by_DryWet.png"),
  fig_treatment,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "Nitrate_DryWet_by_Treatment.pdf"),
  fig_dry_wet,
  width = 14,
  height = 8
)

ggsave(
  file.path(fig_dir, "Nitrate_Treatment_by_DryWet.pdf"),
  fig_treatment,
  width = 12,
  height = 8
)

message("Done. Nitrate figures saved in: ", normalizePath(fig_dir))
