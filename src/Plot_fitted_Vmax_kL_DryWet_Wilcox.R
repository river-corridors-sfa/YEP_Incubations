rm(list = ls())

library(tidyverse)

# Run from repo root or from src/.
if (!dir.exists("modeling_outputs") && dir.exists("../modeling_outputs")) {
  setwd("..")
}

model_output_path <- "modeling_outputs/YEP_Complete_Analysis_Following_Paper_Enhanced.csv"
fig_dir <- "Figures"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

sediment_colors <- c(
  "Dry Sediments" = "#56B4E9",
  "Wet Sediments" = "#1B4F72"
)

param_order <- c("Vmax", "kL")
treatment_order <- c("Control", "Unburned + DOC", "High Burn + DOC")
sediment_order <- c("Dry Sediments", "Wet Sediments")
parameter_labels <- c(
  "Vmax" = "Biotic O₂ consumption\nat saturation (Vmax)",
  "kL" = "Chemical reaction\nrate constant (kL)"
)

model_fits <- read_csv(model_output_path, show_col_types = FALSE) %>%
  filter(Treatment_Type != "Unknown") %>%
  mutate(
    Sediment_Type = case_when(
      str_detect(Treatment_Type, "^Dry_") ~ "Dry Sediments",
      str_detect(Treatment_Type, "^Wet_") ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      str_detect(Treatment_Type, "Control") ~ "Control",
      str_detect(Treatment_Type, "Unburned") ~ "Unburned + DOC",
      str_detect(Treatment_Type, "HighBurn") ~ "High Burn + DOC",
      TRUE ~ NA_character_
    ),
    Sediment_Type = factor(Sediment_Type, levels = sediment_order),
    Treatment = factor(Treatment, levels = treatment_order)
  ) %>%
  select(
    Sample_Name,
    Treatment_Type,
    Treatment,
    Sediment_Type,
    Combined_Vmax_per_h,
    Combined_kL_per_h
  ) %>%
  filter(!is.na(Treatment), !is.na(Sediment_Type)) %>%
  pivot_longer(
    cols = c(Combined_Vmax_per_h, Combined_kL_per_h),
    names_to = "Parameter",
    values_to = "Fitted_Value"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      Combined_Vmax_per_h = "Vmax",
      Combined_kL_per_h = "kL"
    ),
    Parameter = factor(Parameter, levels = param_order)
  ) %>%
  filter(!is.na(Fitted_Value))

wilcox_stats <- model_fits %>%
  group_by(Parameter) %>%
  summarise(
    n_dry = sum(Sediment_Type == "Dry Sediments"),
    n_wet = sum(Sediment_Type == "Wet Sediments"),
    p = if (n_dry > 0 && n_wet > 0) {
      wilcox.test(Fitted_Value ~ Sediment_Type, exact = FALSE)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p) ~ "Wilcox p = NA",
      p < 0.001 ~ "Wilcox p < 0.001",
      TRUE ~ sprintf("Wilcox p = %.3f", p)
    )
  )

write_csv(
  wilcox_stats,
  file.path(fig_dir, "Fitted_Vmax_kL_DryWet_Wilcox_stats.csv")
)

label_positions <- model_fits %>%
  group_by(Parameter) %>%
  summarise(
    y_min = min(Fitted_Value, na.rm = TRUE),
    y_max = max(Fitted_Value, na.rm = TRUE),
    y_rng = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    y_rng = if_else(y_rng == 0, pmax(abs(y_max), 1), y_rng),
    x = 0.65,
    y = y_max + 0.28 * y_rng
  )

wilcox_labels <- wilcox_stats %>%
  left_join(label_positions, by = "Parameter")

pub_theme <- theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black"),
    strip.background = element_rect(fill = "grey85", color = "grey20"),
    strip.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 24, 10, 10)
  )

fig_combined <- ggplot(model_fits, aes(x = Sediment_Type, y = Fitted_Value, fill = Sediment_Type)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 2) +
  geom_blank(
    data = wilcox_labels,
    aes(x = 1, y = y),
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
  scale_fill_manual(values = sediment_colors) +
  scale_x_discrete(labels = c("Dry Sediments" = "Dry\nSediments", "Wet Sediments" = "Wet\nSediments")) +
  labs(
    x = NULL,
    y = expression("Fitted rate constant (" * h^{-1} * ")"),
    fill = "Sediment Type"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(
    ~Parameter,
    nrow = 1,
    scales = "free_y",
    labeller = labeller(Parameter = parameter_labels)
  ) +
  pub_theme

ggsave(
  file.path(fig_dir, "Fitted_Vmax_kL_DryWet_Wilcox.png"),
  fig_combined,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "Fitted_Vmax_kL_DryWet_Wilcox.pdf"),
  fig_combined,
  width = 8,
  height = 5
)
