# ============================================================
# FULL SCRIPT (pub-ready + HOLM posthoc + NO/LESS OVERLAP STATS)
#
# Key behavior:
# - For 3-group treatment plots (Respiration, NPOC, CO2):
#     * Always show Kruskal-Wallis p-value per facet
#     * ONLY show Dunn posthoc (HOLM) if KW p < 0.05 (per facet)
# - For 2-group comparisons (Dry vs Wet within each treatment):
#     * Wilcoxon per treatment; Bonferroni across treatments
#     * Show a single label per facet (no brackets) to reduce clutter
# - Adds extra plotting space automatically for stats (geom_blank + free_y)
# - Disables clipping so labels won’t be cut off
# - Saves figures to ./Figures (creates folder if missing)
# ============================================================

rm(list = ls())

# ---------- Libraries ----------
library(tidyverse)
library(ggplot2)
library(rstatix)

# ---------- Colors ----------
treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

sediment_colors <- c(
  "Dry Sediments" = "#56B4E9",
  "Wet Sediments" = "#1B4F72"
)

# ---------- Publication theme ----------
pub_theme <- theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 16, color = "black"),
    strip.text   = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 30, 10, 10)
  )

# Wrap x labels (instead of angled)
wrap_x <- function(p, width = 14) {
  p + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = width))
}

# ---------- Output folder ----------
fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---------- Vectorized p formatter ----------
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ sprintf("p = %.3f", p)
  )
}

# ---------- Helper: per-facet y-range (handles negative values) ----------
facet_yrange <- function(df, y_col, facet_col) {
  df %>%
    group_by(.data[[facet_col]]) %>%
    summarise(
      y_min = min(.data[[y_col]], na.rm = TRUE),
      y_max = max(.data[[y_col]], na.rm = TRUE),
      y_rng = y_max - y_min,
      .groups = "drop"
    ) %>%
    mutate(y_rng = ifelse(y_rng == 0, abs(y_max), y_rng))
}

# ---------- Build Dunn bracket coordinates using y-range spacing ----------
make_brackets_range <- function(test_df, data_df, x_factor, y_col, facet_col,
                                step_increase = 0.14, tip_frac = 0.06, label_col = "p.adj.signif") {
  
  x_levels <- levels(data_df[[x_factor]])
  
  out <- test_df %>%
    mutate(
      x1 = match(.data$group1, x_levels),
      x2 = match(.data$group2, x_levels),
      xmid = (x1 + x2) / 2,
      label_txt = .data[[label_col]]
    )
  
  yr <- facet_yrange(data_df, y_col = y_col, facet_col = facet_col)
  
  out %>%
    left_join(yr, by = setNames(facet_col, facet_col)) %>%
    group_by(.data[[facet_col]]) %>%
    mutate(
      y.position = y_max + (step_increase * y_rng) * row_number(),
      y.tip      = y.position - tip_frac * y_rng
    ) %>%
    ungroup()
}

# ---------- Add manual Dunn brackets ----------
add_dunn_brackets <- function(p, br_df, label_size = 6, line_width = 0.7) {
  p +
    geom_segment(
      data = br_df,
      aes(x = x1, xend = x2, y = y.position, yend = y.position),
      inherit.aes = FALSE,
      linewidth = line_width
    ) +
    geom_segment(
      data = br_df,
      aes(x = x1, xend = x1, y = y.tip, yend = y.position),
      inherit.aes = FALSE,
      linewidth = line_width
    ) +
    geom_segment(
      data = br_df,
      aes(x = x2, xend = x2, y = y.tip, yend = y.position),
      inherit.aes = FALSE,
      linewidth = line_width
    ) +
    geom_text(
      data = br_df,
      aes(x = xmid, y = y.position, label = label_txt),
      inherit.aes = FALSE,
      vjust = -0.35,
      size = label_size,
      fontface = "bold"
    )
}

# ---------- Add KW label (top-left per facet) ----------
add_kw_labels <- function(p, kw_df, facet_col, y_col, data_df, label_size = 6) {
  yr <- facet_yrange(data_df, y_col = y_col, facet_col = facet_col)
  
  kw_anno <- kw_df %>%
    left_join(yr, by = setNames(facet_col, facet_col)) %>%
    mutate(
      x = 0.7,
      y = y_max + 0.30 * y_rng
    )
  
  p +
    geom_text(
      data = kw_anno,
      aes(x = x, y = y, label = kw_label),
      inherit.aes = FALSE,
      hjust = 0,
      size = label_size,
      fontface = "bold"
    )
}

# ---------- Reserve y-space for KW/Dunn annotations + enable free_y facets + no clipping ----------
add_space_for_annotations <- function(p, data_df, y_col, facet_col, br_df = NULL, kw_df = NULL,
                                      top_pad_frac = 0.25) {
  yr <- facet_yrange(data_df, y_col = y_col, facet_col = facet_col)
  
  req <- yr %>%
    transmute(!!facet_col := .data[[facet_col]],
              req_top = y_max + top_pad_frac * y_rng)
  
  if (!is.null(kw_df) && nrow(kw_df) > 0) {
    req <- req %>%
      left_join(yr, by = setNames(facet_col, facet_col)) %>%
      mutate(req_top = pmax(req_top, y_max + 0.40 * y_rng)) %>%
      select(all_of(facet_col), req_top)
  }
  
  if (!is.null(br_df) && nrow(br_df) > 0) {
    br_top <- br_df %>%
      group_by(.data[[facet_col]]) %>%
      summarise(req_top = max(y.position, na.rm = TRUE), .groups = "drop")
    req <- req %>%
      full_join(br_top, by = facet_col, suffix = c("", "_br")) %>%
      mutate(req_top = pmax(req_top, req_top_br, na.rm = TRUE)) %>%
      select(all_of(facet_col), req_top)
  }
  
  p +
    coord_cartesian(clip = "off") +
    facet_wrap(as.formula(paste0("~", facet_col)), nrow = 1, scales = "free_y") +
    geom_blank(
      data = req,
      aes(x = 1, y = req_top),
      inherit.aes = FALSE
    )
}

# ============================================================
# 1) LOAD DATA
# ============================================================

# Respiration (your new file)
respiration_data <- read.csv("YEP_INC_ReadyForBoye_20260213.csv")

# NPOC / CO2 from the YEP download folder
npoc_data <- read.csv(
  "YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_NPOC_TN.csv",
  skip = 2
) %>% filter(grepl("YEP", Sample_Name))

co2_data <- read.csv(
  "YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_CO2.csv",
  skip = 2
) %>% filter(grepl("YEP", Sample_Name))

# ============================================================
# 2) CLEAN DATA
# ============================================================

# Respiration clean (robust parsing from Sample_Name)
respiration_clean <- respiration_data %>%
  select(
    Sample_Name,
    Respiration_Rate_mg_DO_per_L_per_H,
    Respiration_Rate_mg_DO_per_kg_per_H
  ) %>%
  mutate(
    Respiration_Rate_mg_DO_per_L_per_H  = as.numeric(Respiration_Rate_mg_DO_per_L_per_H),
    Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)
  ) %>%
  mutate(
    Condition = case_when(
      stringr::str_detect(Sample_Name, "-H") ~ "H",
      stringr::str_detect(Sample_Name, "-U") ~ "U",
      stringr::str_detect(Sample_Name, "-S") ~ "S",
      TRUE ~ NA_character_
    ),
    Sediment_Type = case_when(
      stringr::str_detect(Sample_Name, "^YEP1") ~ "Dry Sediments",
      stringr::str_detect(Sample_Name, "^YEP2") ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    DOC_Treatment_short = case_when(
      Condition == "S" ~ "Control",
      Condition == "U" ~ "Unburned + DOC",
      Condition == "H" ~ "High Burn + DOC",
      TRUE ~ NA_character_
    ),
    DOC_Treatment_short = factor(
      DOC_Treatment_short,
      levels = c("Control", "Unburned + DOC", "High Burn + DOC")
    ),
    Resp_kg = -abs(Respiration_Rate_mg_DO_per_kg_per_H),
    Resp_L  = -abs(Respiration_Rate_mg_DO_per_L_per_H)
  ) %>%
  filter(!is.na(Sediment_Type), !is.na(DOC_Treatment_short)) %>%
  filter(!is.na(Resp_kg), Resp_kg != -9999)

# NPOC clean
npoc_final <- npoc_data %>%
  filter(!is.na(Sample_Name), Sample_Name != "", !stringr::str_detect(Sample_Name, "^#")) %>%
  select(Sample_Name, Extractable_NPOC_mg_per_kg) %>%
  filter(!is.na(Extractable_NPOC_mg_per_kg), Extractable_NPOC_mg_per_kg != "-9999") %>%
  mutate(Extractable_NPOC_mg_per_kg = as.numeric(Extractable_NPOC_mg_per_kg)) %>%
  mutate(
    Site = case_when(
      stringr::str_detect(Sample_Name, "^YEP1") ~ "YEP1",
      stringr::str_detect(Sample_Name, "^YEP2") ~ "YEP2",
      TRUE ~ NA_character_
    ),
    Condition = case_when(
      stringr::str_detect(Sample_Name, "-H") ~ "H",
      stringr::str_detect(Sample_Name, "-U") ~ "U",
      stringr::str_detect(Sample_Name, "-S") ~ "S",
      TRUE ~ NA_character_
    ),
    Sediment_Type = case_when(
      Site == "YEP1" ~ "Dry Sediments",
      Site == "YEP2" ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    DOC_Treatment_short = case_when(
      Condition == "S" ~ "Control",
      Condition == "U" ~ "Unburned + DOC",
      Condition == "H" ~ "High Burn + DOC",
      TRUE ~ NA_character_
    ),
    DOC_Treatment_short = factor(
      DOC_Treatment_short,
      levels = c("Control", "Unburned + DOC", "High Burn + DOC")
    )
  ) %>%
  filter(!is.na(Sediment_Type), !is.na(DOC_Treatment_short))

# CO2 clean
co2_final <- co2_data %>%
  filter(!is.na(Sample_Name), Sample_Name != "", !stringr::str_detect(Sample_Name, "^#")) %>%
  select(Sample_Name, Partial_Pressure_CO2_moles_per_L) %>%
  filter(!is.na(Partial_Pressure_CO2_moles_per_L), Partial_Pressure_CO2_moles_per_L != "-9999") %>%
  mutate(Partial_Pressure_CO2_moles_per_L = as.numeric(Partial_Pressure_CO2_moles_per_L)) %>%
  mutate(
    Site = case_when(
      stringr::str_detect(Sample_Name, "^YEP1") ~ "YEP1",
      stringr::str_detect(Sample_Name, "^YEP2") ~ "YEP2",
      TRUE ~ NA_character_
    ),
    Condition = case_when(
      stringr::str_detect(Sample_Name, "-H") ~ "H",
      stringr::str_detect(Sample_Name, "-U") ~ "U",
      stringr::str_detect(Sample_Name, "-S") ~ "S",
      TRUE ~ NA_character_
    ),
    Sediment_Type = case_when(
      Site == "YEP1" ~ "Dry Sediments",
      Site == "YEP2" ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    DOC_Treatment_short = case_when(
      Condition == "S" ~ "Control",
      Condition == "U" ~ "Unburned + DOC",
      Condition == "H" ~ "High Burn + DOC",
      TRUE ~ NA_character_
    ),
    DOC_Treatment_short = factor(
      DOC_Treatment_short,
      levels = c("Control", "Unburned + DOC", "High Burn + DOC")
    )
  ) %>%
  filter(!is.na(Sediment_Type), !is.na(DOC_Treatment_short))

# ============================================================
# 3) PLOTS + STATS
# ============================================================

# ---------- A) Respiration per kg: KW always; Dunn(HOLM) only if KW sig ----------
kw_resp <- respiration_clean %>%
  group_by(Sediment_Type) %>%
  kruskal_test(Resp_kg ~ DOC_Treatment_short) %>%
  ungroup() %>%
  mutate(
    kw_p = p,
    kw_label = paste0("KW ", fmt_p(kw_p)),
    kw_sig = kw_p < 0.05
  )

dunn_resp <- respiration_clean %>%
  semi_join(kw_resp %>% filter(kw_sig) %>% select(Sediment_Type), by = "Sediment_Type") %>%
  group_by(Sediment_Type) %>%
  dunn_test(Resp_kg ~ DOC_Treatment_short, p.adjust.method = "holm") %>%
  add_significance("p.adj") %>%
  ungroup()

dunn_resp_br <- NULL
if (nrow(dunn_resp) > 0) {
  dunn_resp_br <- make_brackets_range(
    test_df = dunn_resp,
    data_df = respiration_clean,
    x_factor = "DOC_Treatment_short",
    y_col = "Resp_kg",
    facet_col = "Sediment_Type",
    step_increase = 0.16,
    tip_frac = 0.06,
    label_col = "p.adj.signif"
  )
}

fig_resp_treat <- ggplot(respiration_clean,
                         aes(x = DOC_Treatment_short, y = Resp_kg, fill = DOC_Treatment_short)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 2) +
  scale_fill_manual(values = treatment_colors) +
  labs(x = NULL,
       y = "Respiration Rate (mg O\u2082 kg\u207b\u00b9 h\u207b\u00b9)",
       fill = "Treatment") +
  pub_theme

fig_resp_treat <- wrap_x(fig_resp_treat, width = 12)
fig_resp_treat <- add_space_for_annotations(
  fig_resp_treat, respiration_clean, "Resp_kg", "Sediment_Type",
  br_df = dunn_resp_br, kw_df = kw_resp, top_pad_frac = 0.35
)
fig_resp_treat <- add_kw_labels(fig_resp_treat, kw_resp, "Sediment_Type", "Resp_kg", respiration_clean, label_size = 6)
if (!is.null(dunn_resp_br) && nrow(dunn_resp_br) > 0) {
  fig_resp_treat <- add_dunn_brackets(fig_resp_treat, dunn_resp_br, label_size = 6, line_width = 0.7)
}
print(fig_resp_treat)

# ---------- B) Dry vs Wet within each treatment: Wilcoxon + Bonferroni; label only ----------
wilcox_resp <- respiration_clean %>%
  group_by(DOC_Treatment_short) %>%
  wilcox_test(Resp_kg ~ Sediment_Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  ungroup() %>%
  mutate(w_label = paste0("Wilcox ", fmt_p(p.adj)))

yr_w <- facet_yrange(respiration_clean, "Resp_kg", "DOC_Treatment_short")

wilcox_lab <- wilcox_resp %>%
  left_join(yr_w, by = "DOC_Treatment_short") %>%
  mutate(x = 0.7, y = y_max + 0.30 * y_rng)

fig_resp_sed <- ggplot(respiration_clean,
                       aes(x = Sediment_Type, y = Resp_kg, fill = Sediment_Type)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 2) +
  scale_fill_manual(values = sediment_colors) +
  labs(x = NULL,
       y = "Respiration Rate (mg O\u2082 kg\u207b\u00b9 h\u207b\u00b9)",
       fill = "Sediment Type") +
  pub_theme

fig_resp_sed <- wrap_x(fig_resp_sed, width = 12) +
  coord_cartesian(clip = "off") +
  facet_wrap(~DOC_Treatment_short, nrow = 1, scales = "free_y") +
  geom_blank(data = wilcox_lab, aes(x = 1, y = y), inherit.aes = FALSE) +
  geom_text(
    data = wilcox_lab,
    aes(x = x, y = y, label = w_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 6,
    fontface = "bold"
  )

print(fig_resp_sed)

# ---------- D) NPOC: KW always; Dunn(HOLM) only if KW sig ----------
kw_npoc <- npoc_final %>%
  group_by(Sediment_Type) %>%
  kruskal_test(Extractable_NPOC_mg_per_kg ~ DOC_Treatment_short) %>%
  ungroup() %>%
  mutate(
    kw_p = p,
    kw_label = paste0("KW ", fmt_p(kw_p)),
    kw_sig = kw_p < 0.05
  )

dunn_npoc <- npoc_final %>%
  semi_join(kw_npoc %>% filter(kw_sig) %>% select(Sediment_Type), by = "Sediment_Type") %>%
  group_by(Sediment_Type) %>%
  dunn_test(Extractable_NPOC_mg_per_kg ~ DOC_Treatment_short, p.adjust.method = "holm") %>%
  add_significance("p.adj") %>%
  ungroup()

dunn_npoc_br <- NULL
if (nrow(dunn_npoc) > 0) {
  dunn_npoc_br <- make_brackets_range(
    test_df = dunn_npoc,
    data_df = npoc_final,
    x_factor = "DOC_Treatment_short",
    y_col = "Extractable_NPOC_mg_per_kg",
    facet_col = "Sediment_Type",
    step_increase = 0.16,
    tip_frac = 0.06,
    label_col = "p.adj.signif"
  )
}

fig_npoc <- ggplot(npoc_final,
                   aes(x = DOC_Treatment_short, y = Extractable_NPOC_mg_per_kg, fill = DOC_Treatment_short)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 2) +
  scale_fill_manual(values = treatment_colors) +
  labs(x = NULL,
       y = "NPOC (mg kg\u207b\u00b9 dry sediment)",
       fill = "Treatment") +
  pub_theme

fig_npoc <- wrap_x(fig_npoc, width = 12)
fig_npoc <- add_space_for_annotations(
  fig_npoc, npoc_final, "Extractable_NPOC_mg_per_kg", "Sediment_Type",
  br_df = dunn_npoc_br, kw_df = kw_npoc, top_pad_frac = 0.35
)
fig_npoc <- add_kw_labels(fig_npoc, kw_npoc, "Sediment_Type", "Extractable_NPOC_mg_per_kg", npoc_final, label_size = 6)
if (!is.null(dunn_npoc_br) && nrow(dunn_npoc_br) > 0) {
  fig_npoc <- add_dunn_brackets(fig_npoc, dunn_npoc_br, label_size = 6, line_width = 0.7)
}
print(fig_npoc)

# ---------- E) CO2: KW always; Dunn(HOLM) only if KW sig ----------
kw_co2 <- co2_final %>%
  group_by(Sediment_Type) %>%
  kruskal_test(Partial_Pressure_CO2_moles_per_L ~ DOC_Treatment_short) %>%
  ungroup() %>%
  mutate(
    kw_p = p,
    kw_label = paste0("KW ", fmt_p(kw_p)),
    kw_sig = kw_p < 0.05
  )

dunn_co2 <- co2_final %>%
  semi_join(kw_co2 %>% filter(kw_sig) %>% select(Sediment_Type), by = "Sediment_Type") %>%
  group_by(Sediment_Type) %>%
  dunn_test(Partial_Pressure_CO2_moles_per_L ~ DOC_Treatment_short, p.adjust.method = "holm") %>%
  add_significance("p.adj") %>%
  ungroup()

dunn_co2_br <- NULL
if (nrow(dunn_co2) > 0) {
  dunn_co2_br <- make_brackets_range(
    test_df = dunn_co2,
    data_df = co2_final,
    x_factor = "DOC_Treatment_short",
    y_col = "Partial_Pressure_CO2_moles_per_L",
    facet_col = "Sediment_Type",
    step_increase = 0.18,
    tip_frac = 0.06,
    label_col = "p.adj.signif"
  )
}

fig_co2 <- ggplot(co2_final,
                  aes(x = DOC_Treatment_short, y = Partial_Pressure_CO2_moles_per_L, fill = DOC_Treatment_short)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 2) +
  scale_fill_manual(values = treatment_colors) +
  labs(x = NULL,
       y = expression(CO[2]~"(mol L"^{-1}*")"),
       fill = "Treatment") +
  pub_theme

fig_co2 <- wrap_x(fig_co2, width = 12)
fig_co2 <- add_space_for_annotations(
  fig_co2, co2_final, "Partial_Pressure_CO2_moles_per_L", "Sediment_Type",
  br_df = dunn_co2_br, kw_df = kw_co2, top_pad_frac = 0.40
)
fig_co2 <- add_kw_labels(fig_co2, kw_co2, "Sediment_Type", "Partial_Pressure_CO2_moles_per_L", co2_final, label_size = 6)
if (!is.null(dunn_co2_br) && nrow(dunn_co2_br) > 0) {
  fig_co2 <- add_dunn_brackets(fig_co2, dunn_co2_br, label_size = 6, line_width = 0.7)
}
# scientific formatting after add_space_for_annotations (keeps expansion)
fig_co2 <- fig_co2 + scale_y_continuous(labels = scales::scientific)

print(fig_co2)

# ============================================================
# 4) SAVE FIGURES (Figures folder)
# ============================================================
ggsave(file.path(fig_dir, "Respiration_Treatment_PubReady_Holm.png"), fig_resp_treat, width = 14, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Respiration_DryVsWet_Wilcox_PubReady.png"), fig_resp_sed, width = 14, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "NPOC_Treatment_PubReady_Holm.png"), fig_npoc, width = 14, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "CO2_Treatment_PubReady_Holm.png"), fig_co2, width = 14, height = 8, dpi = 300)

message("Done. Figures saved in: ", normalizePath(fig_dir))
