rm(list = ls())

library(tidyverse)
library(rstatix)

# Run from repo root or from src/.
if (!dir.exists("modeling_outputs") && dir.exists("../modeling_outputs")) {
  setwd("..")
}

model_output_path <- "modeling_outputs/YEP_Complete_Analysis_Following_Paper_Enhanced.csv"
fig_dir <- "Figures"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

treatment_order <- c("Control", "Unburned + DOC", "High Burn + DOC")
sediment_order <- c("Dry Sediments", "Wet Sediments")
param_order <- c("Vmax", "kL")

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

kw_stats <- model_fits %>%
  group_by(Parameter, Sediment_Type) %>%
  kruskal_test(Fitted_Value ~ Treatment) %>%
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
  file.path(fig_dir, "Fitted_Vmax_kL_Treatment_KW_stats.csv")
)

dunn_input <- model_fits %>%
  semi_join(
    kw_stats %>%
      filter(kw_sig) %>%
      select(Parameter, Sediment_Type),
    by = c("Parameter", "Sediment_Type")
  )

if (nrow(dunn_input) > 0) {
  dunn_stats <- dunn_input %>%
    group_by(Parameter, Sediment_Type) %>%
    dunn_test(Fitted_Value ~ Treatment, p.adjust.method = "holm") %>%
    add_significance("p.adj") %>%
    ungroup()
} else {
  dunn_stats <- tibble(
    Parameter = factor(levels = param_order),
    Sediment_Type = factor(levels = sediment_order),
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
  file.path(fig_dir, "Fitted_Vmax_kL_Treatment_Dunn_stats.csv")
)

facet_ranges <- model_fits %>%
  group_by(Parameter, Sediment_Type) %>%
  summarise(
    y_min = min(Fitted_Value, na.rm = TRUE),
    y_max = max(Fitted_Value, na.rm = TRUE),
    y_rng = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(y_rng = if_else(y_rng == 0, pmax(abs(y_max), 1), y_rng))

kw_labels <- kw_stats %>%
  left_join(facet_ranges, by = c("Parameter", "Sediment_Type")) %>%
  mutate(
    x = 0.55,
    y = y_max + 0.20 * y_rng
  )

dunn_brackets <- dunn_stats %>%
  mutate(
    group1 = factor(group1, levels = treatment_order),
    group2 = factor(group2, levels = treatment_order),
    x1 = match(as.character(group1), treatment_order),
    x2 = match(as.character(group2), treatment_order),
    xmid = (x1 + x2) / 2,
    bracket_width = x2 - x1
  ) %>%
  left_join(facet_ranges, by = c("Parameter", "Sediment_Type")) %>%
  arrange(Parameter, Sediment_Type, bracket_width, x1, x2) %>%
  group_by(Parameter, Sediment_Type) %>%
  mutate(
    y.position = y_max + (0.15 + 0.13 * (row_number() - 1)) * y_rng,
    y.tip = y.position - 0.04 * y_rng,
    label = case_when(
      is.na(p.adj) ~ "p = NA",
      p.adj < 0.001 ~ "p < 0.001",
      TRUE ~ sprintf("p = %.3f", p.adj)
    )
  ) %>%
  ungroup()

annotation_tops <- facet_ranges %>%
  mutate(req_top = y_max + 0.35 * y_rng) %>%
  select(Parameter, Sediment_Type, req_top)

if (nrow(dunn_brackets) > 0) {
  bracket_tops <- dunn_brackets %>%
    group_by(Parameter, Sediment_Type) %>%
    summarise(req_top_br = max(y.position + 0.05 * y_rng, na.rm = TRUE), .groups = "drop")

  annotation_tops <- annotation_tops %>%
    left_join(bracket_tops, by = c("Parameter", "Sediment_Type")) %>%
    mutate(req_top = pmax(req_top, req_top_br, na.rm = TRUE)) %>%
    select(Parameter, Sediment_Type, req_top)
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

fig_vmax <- model_fits %>%
  filter(Parameter == "Vmax") %>%
  ggplot(aes(x = Treatment, y = Fitted_Value, fill = Treatment)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 2) +
  geom_blank(
    data = annotation_tops %>% filter(Parameter == "Vmax"),
    aes(x = 1, y = req_top),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = kw_labels %>% filter(Parameter == "Vmax"),
    aes(x = x, y = y, label = kw_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 5,
    fontface = "bold"
  ) +
  geom_segment(
    data = dunn_brackets %>% filter(Parameter == "Vmax"),
    aes(x = x1, xend = x2, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = dunn_brackets %>% filter(Parameter == "Vmax"),
    aes(x = x1, xend = x1, y = y.tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = dunn_brackets %>% filter(Parameter == "Vmax"),
    aes(x = x2, xend = x2, y = y.tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_text(
    data = dunn_brackets %>% filter(Parameter == "Vmax"),
    aes(x = xmid, y = y.position, label = label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
  labs(
    x = NULL,
    y = expression("Fitted " * V[max] * " (" * h^{-1} * ")"),
    fill = "Treatment"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(~Sediment_Type, nrow = 1, scales = "free_y") +
  pub_theme

fig_kl <- model_fits %>%
  filter(Parameter == "kL") %>%
  ggplot(aes(x = Treatment, y = Fitted_Value, fill = Treatment)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 2) +
  geom_blank(
    data = annotation_tops %>% filter(Parameter == "kL"),
    aes(x = 1, y = req_top),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = kw_labels %>% filter(Parameter == "kL"),
    aes(x = x, y = y, label = kw_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 5,
    fontface = "bold"
  ) +
  geom_segment(
    data = dunn_brackets %>% filter(Parameter == "kL"),
    aes(x = x1, xend = x2, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = dunn_brackets %>% filter(Parameter == "kL"),
    aes(x = x1, xend = x1, y = y.tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = dunn_brackets %>% filter(Parameter == "kL"),
    aes(x = x2, xend = x2, y = y.tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_text(
    data = dunn_brackets %>% filter(Parameter == "kL"),
    aes(x = xmid, y = y.position, label = label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
  labs(
    x = NULL,
    y = expression("Fitted " * k[L] * " (" * h^{-1} * ")"),
    fill = "Treatment"
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(~Sediment_Type, nrow = 1, scales = "free_y") +
  pub_theme

fig_combined <- ggplot(model_fits, aes(x = Treatment, y = Fitted_Value, fill = Treatment)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 2) +
  geom_blank(
    data = annotation_tops,
    aes(x = 1, y = req_top),
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
  geom_segment(
    data = dunn_brackets,
    aes(x = x1, xend = x2, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = dunn_brackets,
    aes(x = x1, xend = x1, y = y.tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = dunn_brackets,
    aes(x = x2, xend = x2, y = y.tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_text(
    data = dunn_brackets,
    aes(x = xmid, y = y.position, label = label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 4.2,
    fontface = "bold"
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
  labs(
    x = NULL,
    y = expression("Fitted value (" * h^{-1} * ")"),
    fill = "Treatment"
  ) +
  coord_cartesian(clip = "off") +
  facet_grid(Parameter ~ Sediment_Type, scales = "free_y") +
  pub_theme

print(fig_vmax)
print(fig_kl)
print(fig_combined)

ggsave(
  file.path(fig_dir, "Fitted_Vmax_kL_Treatment_KW_Dunn.png"),
  fig_combined,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "Fitted_Vmax_kL_Treatment_KW_Dunn.pdf"),
  fig_combined,
  width = 14,
  height = 10
)
