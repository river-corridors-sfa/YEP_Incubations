# =========================
# NMDS + PERMANOVA for incubated OM data
# =========================

# Packages
library(tidyverse)
library(vegan)
library(ggplot2)

# -------------------------
# 1. Read files
# -------------------------
om_data <- read_csv("C:/Users/gara009/OneDrive - PNNL/RC-SFA - Documents/Study_YEP/FTICR/05_PublishReadyData/Processed Data/YEP_Updated_SN12-Processed_Data.csv")

# -------------------------
# 2. Clean sample names
# -------------------------
sample_cols <- colnames(om_data)[-1]

clean_sample_names <- sample_cols %>%
  str_remove("_p075\\.corems$") %>%
  str_remove("_p2\\.corems$") %>%
  str_remove("\\.corems$")

colnames(om_data) <- c("Calibrated_mz", clean_sample_names)

# -------------------------
# 3. Keep only sediment-incubated samples
# -------------------------
sediment_samples <- colnames(om_data)[-1] %>%
  .[str_detect(., "^YEP[0-9]+[A-Z]_SIR-[HSU][0-9]+$")]

om_sub <- om_data %>%
  select(Calibrated_mz, all_of(sediment_samples))

# -------------------------
# 4. Convert intensities to presence/absence
# -------------------------
om_pa <- om_sub %>%
  column_to_rownames("Calibrated_mz") %>%
  t() %>%
  as.data.frame()

om_pa[] <- lapply(om_pa, as.numeric)
om_pa[is.na(om_pa)] <- 0
om_pa_bin <- ifelse(om_pa > 0, 1, 0) %>% as.data.frame()

# -------------------------
# 5. Create sample info from sample names
# -------------------------
sample_info <- tibble(Sample_Name = rownames(om_pa_bin)) %>%
  mutate(
    Treatment_code = str_extract(Sample_Name, "(?<=_SIR-)[HSU]"),
    Treatment = case_when(
      Treatment_code == "S" ~ "Control",
      Treatment_code == "U" ~ "Unburned + DOC",
      Treatment_code == "H" ~ "High Burn + DOC"
    ),
    Moisture = case_when(
      str_detect(Sample_Name, "^YEP1") ~ "Dry",
      str_detect(Sample_Name, "^YEP2") ~ "Wet"
    ),
    Treatment = factor(
      Treatment,
      levels = c("Control", "Unburned + DOC", "High Burn + DOC")
    ),
    Moisture = factor(Moisture, levels = c("Dry", "Wet"))
  )

# -------------------------
# 6. Distance matrix
# -------------------------
dist_jaccard <- vegdist(om_pa_bin, method = "jaccard", binary = TRUE)

# -------------------------
# 7. Helper function for p-value formatting
# -------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  paste0("= ", sprintf("%.3f", p))
}

# -------------------------
# 8. NMDS for all samples
# -------------------------
set.seed(123)

nmds_all <- metaMDS(
  dist_jaccard,
  k = 2,
  trymax = 200,
  trace = FALSE
)

scores_all <- scores(nmds_all, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_info, by = "Sample_Name")

# -------------------------
# 9. PERMANOVA for all samples
# -------------------------
set.seed(123)

permanova_all <- adonis2(
  dist_jaccard ~ Treatment + Moisture,
  data = sample_info,
  permutations = 999
)

print("PERMANOVA - all samples:")
print(as.data.frame(permanova_all))

all_tab <- as.data.frame(permanova_all)

r2_treat_all <- all_tab$R2[1]
p_treat_all  <- all_tab$`Pr(>F)`[1]

r2_moist_all <- all_tab$R2[2]
p_moist_all  <- all_tab$`Pr(>F)`[2]

label_all <- paste0(
  "Stress = ", round(nmds_all$stress, 3),
  "\nTreatment: R² = ", round(r2_treat_all, 3), ", p ", fmt_p(p_treat_all),
  "\nMoisture: R² = ", round(r2_moist_all, 3), ", p ", fmt_p(p_moist_all)
)

# -------------------------
# 10. Subset dry and wet
# -------------------------
dry_samples <- sample_info %>%
  filter(Moisture == "Dry") %>%
  pull(Sample_Name)

wet_samples <- sample_info %>%
  filter(Moisture == "Wet") %>%
  pull(Sample_Name)

om_pa_dry <- om_pa_bin[dry_samples, , drop = FALSE]
om_pa_wet <- om_pa_bin[wet_samples, , drop = FALSE]

sample_info_dry <- sample_info %>%
  filter(Moisture == "Dry")

sample_info_wet <- sample_info %>%
  filter(Moisture == "Wet")

# -------------------------
# 11. Dry-only NMDS + PERMANOVA
# -------------------------
dist_dry <- vegdist(om_pa_dry, method = "jaccard", binary = TRUE)

set.seed(123)
nmds_dry <- metaMDS(
  dist_dry,
  k = 2,
  trymax = 200,
  trace = FALSE
)

scores_dry <- scores(nmds_dry, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_info_dry, by = "Sample_Name")

set.seed(123)
permanova_dry <- adonis2(
  dist_dry ~ Treatment,
  data = sample_info_dry,
  permutations = 999
)

print("PERMANOVA - dry samples only:")
print(as.data.frame(permanova_dry))

dry_tab <- as.data.frame(permanova_dry)
r2_dry <- dry_tab$R2[1]
p_dry  <- dry_tab$`Pr(>F)`[1]

label_dry <- paste0(
  "Stress = ", round(nmds_dry$stress, 3),
  "\nTreatment: R² = ", round(r2_dry, 3), ", p ", fmt_p(p_dry)
)

# -------------------------
# 12. Wet-only NMDS + PERMANOVA
# -------------------------
dist_wet <- vegdist(om_pa_wet, method = "jaccard", binary = TRUE)

set.seed(123)
nmds_wet <- metaMDS(
  dist_wet,
  k = 2,
  trymax = 200,
  trace = FALSE
)

scores_wet <- scores(nmds_wet, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_info_wet, by = "Sample_Name")

set.seed(123)
permanova_wet <- adonis2(
  dist_wet ~ Treatment,
  data = sample_info_wet,
  permutations = 999
)

print("PERMANOVA - wet samples only:")
print(as.data.frame(permanova_wet))

wet_tab <- as.data.frame(permanova_wet)
r2_wet <- wet_tab$R2[1]
p_wet  <- wet_tab$`Pr(>F)`[1]

label_wet <- paste0(
  "Stress = ", round(nmds_wet$stress, 3),
  "\nTreatment: R² = ", round(r2_wet, 3), ", p ", fmt_p(p_wet)
)

# -------------------------
# 13. Plot settings
# -------------------------
treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

moisture_shapes <- c(
  "Dry" = 16,
  "Wet" = 17
)

# -------------------------
# 14. Plot all samples (no ellipses)
# -------------------------
plot_all <- ggplot(
  scores_all,
  aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Moisture)
) +
  geom_point(size = 3.2) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = moisture_shapes) +
  theme_bw(base_size = 13) +
  labs(
    title = "NMDS of sediment-incubated OM samples",
    subtitle = label_all,
    x = "NMDS1",
    y = "NMDS2",
    color = "Treatment",
    shape = "Sediment moisture"
  )

print(plot_all)

# -------------------------
# 15. Plot dry only (with ellipses)
# -------------------------
plot_dry <- ggplot(
  scores_dry,
  aes(x = NMDS1, y = NMDS2, color = Treatment)
) +
  stat_ellipse(
    aes(fill = Treatment),
    geom = "polygon",
    alpha = 0.18,
    color = NA,
    level = 0.95
  ) +
  geom_point(size = 3.2, shape = 16) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  theme_bw(base_size = 13) +
  labs(
    title = "NMDS of dry sediment-incubated OM samples",
    subtitle = label_dry,
    x = "NMDS1",
    y = "NMDS2",
    color = "Treatment",
    fill = "Treatment"
  )

print(plot_dry)

# -------------------------
# 16. Plot wet only (with ellipses)
# -------------------------
plot_wet <- ggplot(
  scores_wet,
  aes(x = NMDS1, y = NMDS2, color = Treatment)
) +
  stat_ellipse(
    aes(fill = Treatment),
    geom = "polygon",
    alpha = 0.18,
    color = NA,
    level = 0.95
  ) +
  geom_point(size = 3.2, shape = 17) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  theme_bw(base_size = 13) +
  labs(
    title = "NMDS of wet sediment-incubated OM samples",
    subtitle = label_wet,
    x = "NMDS1",
    y = "NMDS2",
    color = "Treatment",
    fill = "Treatment"
  )

print(plot_wet)

# -------------------------
# 17. Save figures
# -------------------------
ggsave("NMDS_all_samples_no_ellipses.png", plot_all, width = 8, height = 6, dpi = 300)
ggsave("NMDS_dry_only_with_ellipses.png", plot_dry, width = 8, height = 6, dpi = 300)
ggsave("NMDS_wet_only_with_ellipses.png", plot_wet, width = 8, height = 6, dpi = 300)