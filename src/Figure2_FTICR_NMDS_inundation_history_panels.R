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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

site_scores_path <- file.path(out_dir, "fticr_pa_nmds_site_scores.csv")
dry_permanova_path <- file.path(out_dir, "fticr_pa_permanova_dry_only.csv")
wet_permanova_path <- file.path(out_dir, "fticr_pa_permanova_wet_only.csv")
sample_metadata_path <- file.path("v2_data", "v2_YEP_Sample_Name_Metadata.csv")
fticr_data_path <- file.path(
  "v2_data", "v2_YEP_Sample_Data", "FTICR",
  "YEP_Sediment_CoreMS_Processed_ICR_Data.csv"
)

needed_files <- c(sample_metadata_path, fticr_data_path)


treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

moisture_shapes <- c(
  "Dry" = 16,
  "Wet" = 17
)

target_sample_material <- "Aqueous sample post-incubation (treatment solution and sediment)"

to_number <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^<", "")
  x[x %in% c("", "NA", "N/A", "-9999")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

sample_design <- read_csv(
  sample_metadata_path,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    Sample_Material = str_squish(Sample_Material),
    Incubation_Day = str_squish(Incubation_Day),
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
    ),
    Treatment = factor(Treatment, levels = names(treatment_colors)),
    Moisture = factor(Moisture, levels = names(moisture_shapes)),
    Incubation_Day = factor(Incubation_Day, levels = c("Day_1", "Day_2", "Day_3"))
  )

fticr_raw <- read_csv(
  fticr_data_path,
  show_col_types = FALSE,
  name_repair = "minimal",
  trim_ws = TRUE
)

sample_cols <- setdiff(names(fticr_raw), "Calibrated_Mass")

fticr_numeric <- fticr_raw %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), to_number))

pa_all <- t(ifelse(as.matrix(fticr_numeric) > 0, 1, 0))
pa_all[is.na(pa_all)] <- 0
storage.mode(pa_all) <- "numeric"
rownames(pa_all) <- sample_cols
colnames(pa_all) <- paste0("mol_", seq_len(nrow(fticr_raw)))

fticr_sample_info <- tibble(Sample_Name = rownames(pa_all)) %>%
  left_join(sample_design, by = "Sample_Name") %>%
  filter(
    Sample_Material == target_sample_material,
    !is.na(Treatment),
    !is.na(Moisture),
    !is.na(Incubation_Day)
  ) %>%
  arrange(Sample_Name)

pa <- pa_all[fticr_sample_info$Sample_Name, , drop = FALSE]
pa <- pa[, colSums(pa) > 0, drop = FALSE]
pa <- pa[rowSums(pa) > 0, , drop = FALSE]

fticr_sample_info <- fticr_sample_info %>%
  filter(Sample_Name %in% rownames(pa)) %>%
  arrange(match(Sample_Name, rownames(pa))) %>%
  mutate(
    Treatment = droplevels(Treatment),
    Moisture = droplevels(Moisture),
    Incubation_Day = droplevels(Incubation_Day)
  )

stopifnot(identical(fticr_sample_info$Sample_Name, rownames(pa)))

dist_jaccard <- vegdist(pa, method = "jaccard", binary = TRUE)

nmds <- metaMDS(
  pa,
  distance = "jaccard",
  binary = TRUE,
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

site_scores <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(fticr_sample_info, by = "Sample_Name") %>%
  mutate(
    Treatment = factor(Treatment, levels = names(treatment_colors)),
    Moisture = factor(Moisture, levels = names(moisture_shapes))
  ) %>%
  filter(!is.na(Treatment), !is.na(Moisture))

write_csv(site_scores, site_scores_path)

run_within_moisture_permanova <- function(moisture_level) {
  sample_info_sub <- fticr_sample_info %>%
    filter(Moisture == moisture_level) %>%
    droplevels()

  sample_names_sub <- sample_info_sub$Sample_Name
  dist_sub <- as.dist(
    as.matrix(dist_jaccard)[sample_names_sub, sample_names_sub]
  )

  set.seed(123)
  adonis2(
    dist_sub ~ Treatment + Incubation_Day,
    data = sample_info_sub,
    permutations = 9999,
    by = "margin"
  )
}

dry_permanova_all <- run_within_moisture_permanova("Dry")
wet_permanova_all <- run_within_moisture_permanova("Wet")

write_csv(
  as.data.frame(dry_permanova_all) %>% rownames_to_column("Term"),
  dry_permanova_path
)

write_csv(
  as.data.frame(wet_permanova_all) %>% rownames_to_column("Term"),
  wet_permanova_path
)

dry_permanova <- as.data.frame(dry_permanova_all) %>%
  rownames_to_column("Term") %>%
  filter(Term == "Treatment") %>%
  mutate(Moisture = "Dry")

wet_permanova <- as.data.frame(wet_permanova_all) %>%
  rownames_to_column("Term") %>%
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

dry_label <- permanova_labels %>% filter(Moisture == "Dry")
wet_label <- permanova_labels %>% filter(Moisture == "Wet")

plot_labels <- tibble(
  x = x_min + 0.03 * (x_max - x_min),
  y = y_max - c(0.06, 0.12, 0.18) * (y_max - y_min),
  label = c(
    paste0("NMDS stress = ", sprintf("%.2f", nmds$stress)),
    paste0(
      "Dry PERMANOVA",
      ": R2 = ", sprintf("%.2f", dry_label %>% pull(R2)),
      ", p ", dry_label %>% pull(p_text)
    ),
    paste0(
      "Wet PERMANOVA",
      ": R2 = ", sprintf("%.2f", wet_label %>% pull(R2)),
      ", p ", wet_label %>% pull(p_text)
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
  file.path(out_dir, "Figure3_fticr_pa_nmds_by_inundation_history.png"),
  nmds_inundation_plot,
  width = 6.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "Figure3_fticr_pa_nmds_by_inundation_history.pdf"),
  nmds_inundation_plot,
  width = 6.5,
  height = 6.5
)
