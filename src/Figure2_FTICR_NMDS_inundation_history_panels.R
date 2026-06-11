# ================================
# Figure 2 FTICR-MS presence/absence NMDS by inundation history
# ================================

library(tidyverse)
library(vegan)

# Run from repo root or from src/.
if (!dir.exists("Figures") && dir.exists("../Figures")) {
  setwd("..")
}

set.seed(123)

out_dir <- file.path("Figures", "fticr_nmds_v2")

site_scores_path <- file.path(out_dir, "fticr_pa_nmds_site_scores.csv")
dry_permanova_path <- file.path(out_dir, "fticr_pa_permanova_dry_only.csv")
wet_permanova_path <- file.path(out_dir, "fticr_pa_permanova_wet_only.csv")
fticr_data_path <- file.path(
  "v2_data", "v2_YEP_Sample_Data", "FTICR",
  "YEP_Sediment_CoreMS_Processed_ICR_Data.csv"
)

needed_files <- c(
  site_scores_path,
  dry_permanova_path,
  wet_permanova_path,
  fticr_data_path
)
missing_files <- needed_files[!file.exists(needed_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing NMDS/PERMANOVA output files. Run ",
    "src/FTICR_NMDS_presence_absence_envfit.R first. Missing: ",
    paste(missing_files, collapse = ", ")
  )
}

treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

moisture_shapes <- c(
  "Dry" = 16,
  "Wet" = 17
)

site_scores <- read_csv(site_scores_path, show_col_types = FALSE) %>%
  mutate(
    Treatment = factor(Treatment, levels = names(treatment_colors)),
    Moisture = factor(Moisture, levels = names(moisture_shapes))
  ) %>%
  filter(!is.na(Treatment), !is.na(Moisture))

fticr_raw <- read_csv(
  fticr_data_path,
  show_col_types = FALSE,
  name_repair = "minimal",
  trim_ws = TRUE
)

fticr_matrix <- fticr_raw %>%
  select(all_of(site_scores$Sample_Name)) %>%
  mutate(
    across(
      everything(),
      ~ {
        x <- as.character(.x)
        x <- str_replace(x, "^<", "")
        x[x %in% c("", "NA", "N/A", "-9999")] <- NA_character_
        suppressWarnings(as.numeric(x))
      }
    )
  )

pa <- t(ifelse(as.matrix(fticr_matrix) > 0, 1, 0))
pa[is.na(pa)] <- 0
storage.mode(pa) <- "numeric"
rownames(pa) <- site_scores$Sample_Name
colnames(pa) <- paste0("mol_", seq_len(nrow(fticr_raw)))
pa <- pa[, colSums(pa) > 0, drop = FALSE]

nmds <- metaMDS(
  pa,
  distance = "jaccard",
  binary = TRUE,
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

dry_permanova <- read_csv(dry_permanova_path, show_col_types = FALSE) %>%
  filter(Term == "Treatment") %>%
  mutate(Moisture = "Dry")

wet_permanova <- read_csv(wet_permanova_path, show_col_types = FALSE) %>%
  filter(Term == "Treatment") %>%
  mutate(Moisture = "Wet")

permanova_labels <- bind_rows(dry_permanova, wet_permanova) %>%
  mutate(
    Moisture = factor(Moisture, levels = levels(site_scores$Moisture)),
    p_text = case_when(
      is.na(`Pr(>F)`) ~ "NA",
      `Pr(>F)` < 0.001 ~ "< 0.001",
      TRUE ~ paste0("= ", sprintf("%.3f", `Pr(>F)`))
    ),
    label = paste0(Moisture, " treatment PERMANOVA")
  )

x_min <- min(site_scores$NMDS1, na.rm = TRUE)
x_max <- max(site_scores$NMDS1, na.rm = TRUE)
y_min <- min(site_scores$NMDS2, na.rm = TRUE)
y_max <- max(site_scores$NMDS2, na.rm = TRUE)

plot_labels <- tibble(
  x = x_min + 0.03 * (x_max - x_min),
  y = y_max - c(0.06, 0.12, 0.18) * (y_max - y_min),
  label = c(
    paste0("NMDS stress = ", sprintf("%.2f", nmds$stress)),
    paste0(
      "Dry PERMANOVA",
      ": R² = ", sprintf("%.2f", permanova_labels %>% filter(Moisture == "Dry") %>% pull(R2)),
      ", p = ", sprintf("%.3f", permanova_labels %>% filter(Moisture == "Dry") %>% pull(`Pr(>F)`))
    ),
    paste0(
      "Wet PERMANOVA",
      ": R² = ", sprintf("%.2f", permanova_labels %>% filter(Moisture == "Wet") %>% pull(R2)),
      ", p = ", sprintf("%.3f", permanova_labels %>% filter(Moisture == "Wet") %>% pull(`Pr(>F)`))
    )
  )
)

write_csv(
  permanova_labels %>%
    select(Moisture, Term, Df, SumOfSqs, R2, F, `Pr(>F)`, label),
  file.path(out_dir, "figure2_fticr_pa_nmds_inundation_history_permanova_labels.csv")
)

nmds_inundation_plot <- ggplot(
  site_scores,
  aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Moisture)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text(
    data = plot_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.4,
    color = "black"
  ) +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  scale_shape_manual(values = moisture_shapes, drop = FALSE) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Treatment",
    shape = "Sediment moisture"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.6),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    aspect.ratio = 1,
    plot.margin = margin(8, 12, 8, 8)
  )

print(nmds_inundation_plot)

ggsave(
  file.path(out_dir, "figure2_fticr_pa_nmds_by_inundation_history.png"),
  nmds_inundation_plot,
  width = 6.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "figure2_fticr_pa_nmds_by_inundation_history.pdf"),
  nmds_inundation_plot,
  width = 6.5,
  height = 6.5
)
