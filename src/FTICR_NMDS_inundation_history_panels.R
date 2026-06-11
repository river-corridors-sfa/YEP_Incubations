# ================================
# FTICR-MS NMDS Analysis
# ================================

library(tidyverse)
library(vegan)

# Setup ----
set.seed(123)

# Set working directory
if (!dir.exists("Figures") && dir.exists("../Figures")) setwd("..")

# Create output folder
out_dir <- file.path("Figures", "fticr_nmds_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define file paths
data_dir <- "v2_data"
metadata_file <- file.path(data_dir, "v2_YEP_Sample_Name_Metadata.csv")
fticr_file <- file.path(data_dir, "v2_YEP_Sample_Data", "FTICR", 
                        "YEP_Sediment_CoreMS_Processed_ICR_Data.csv")

# Define visual mappings
treatment_colors <- c("Control" = "#0072B2", 
                      "Unburned + DOC" = "darkgreen",
                      "High Burn + DOC" = "#8B4513")
moisture_shapes <- c("Dry" = 16, "Wet" = 17)

# Target samples: post-incubation aqueous samples
target_sample_type <- "Aqueous sample post-incubation (treatment solution and sediment)"


# Load and Clean Data ----

## Read metadata and standardize treatment/moisture names
metadata <- read_csv(metadata_file, show_col_types = FALSE, trim_ws = TRUE) %>%
  mutate(
    # Simplify treatment names
    Treatment = case_when(
      str_detect(Treatment, "Columbia River") ~ "Control",
      str_detect(Treatment, "Unburned") ~ "Unburned + DOC",
      str_detect(Treatment, "Burned") ~ "High Burn + DOC"
    ),
    # Standardize moisture categories
    Moisture = str_to_title(str_squish(Sediment_Field_Moisture_Conditions)),
    # Convert to factors with defined order
    Treatment = factor(Treatment, levels = names(treatment_colors)),
    Moisture = factor(Moisture, levels = names(moisture_shapes)),
    Incubation_Day = factor(Incubation_Day, levels = c("Day_1", "Day_2", "Day_3"))
  )

## Read FTICR data (rows = molecules, columns = samples)
fticr_raw <- read_csv(fticr_file, show_col_types = FALSE, trim_ws = TRUE)


# Create Presence/Absence Matrix ----

## Extract sample columns (exclude mass column)
sample_cols <- setdiff(names(fticr_raw), "Calibrated_Mass")

## Convert to presence/absence (1 if detected, 0 if not)
pa_matrix <- fticr_raw %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), ~as.numeric(. > 0 & !is.na(.)))) %>%
  as.matrix() %>%
  t()  # Transpose: rows = samples, columns = molecules

rownames(pa_matrix) <- sample_cols
colnames(pa_matrix) <- paste0("mol_", seq_len(nrow(fticr_raw)))


# Filter Samples ----

## Keep only target samples with complete metadata
sample_info <- tibble(Sample_Name = rownames(pa_matrix)) %>%
  left_join(metadata, by = "Sample_Name") %>%
  filter(
    Sample_Material == target_sample_type,
    !is.na(Treatment), !is.na(Moisture), !is.na(Incubation_Day)
  )

## Filter presence/absence matrix to match
pa_matrix <- pa_matrix[sample_info$Sample_Name, , drop = FALSE]

## Remove molecules never detected & samples with no detections
pa_matrix <- pa_matrix[, colSums(pa_matrix) > 0, drop = FALSE]
pa_matrix <- pa_matrix[rowSums(pa_matrix) > 0, , drop = FALSE]

## Update sample info to match filtered matrix
sample_info <- sample_info %>%
  filter(Sample_Name %in% rownames(pa_matrix)) %>%
  arrange(match(Sample_Name, rownames(pa_matrix)))

stopifnot(identical(sample_info$Sample_Name, rownames(pa_matrix)))


# Run NMDS Ordination ----

## Calculate Jaccard distance (presence/absence similarity)
dist_jaccard <- vegdist(pa_matrix, method = "jaccard", binary = TRUE)

## Run NMDS to create 2D representation
nmds <- metaMDS(pa_matrix, distance = "jaccard", binary = TRUE, k = 2, 
                trymax = 200, autotransform = FALSE, trace = FALSE)

## Extract sample coordinates and add metadata
site_scores <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_info, by = "Sample_Name")

write_csv(site_scores, file.path(out_dir, "fticr_pa_nmds_site_scores.csv"))


# Statistical Testing: PERMANOVA ----
# Test if treatments differ significantly, separately for dry and wet samples

run_permanova <- function(moisture_level) {
  # Subset data for this moisture level
  samples_subset <- sample_info %>% filter(Moisture == moisture_level)
  dist_subset <- as.dist(as.matrix(dist_jaccard)[samples_subset$Sample_Name, 
                                                 samples_subset$Sample_Name])
  
  # Run PERMANOVA
  set.seed(123)
  adonis2(dist_subset ~ Treatment + Incubation_Day, 
          data = samples_subset, permutations = 9999, by = "margin")
}

## Run for dry and wet separately
dry_permanova <- run_permanova("Dry")
wet_permanova <- run_permanova("Wet")

## Save full results
write_csv(as.data.frame(dry_permanova) %>% rownames_to_column("Term"),
          file.path(out_dir, "fticr_pa_permanova_dry_only.csv"))
write_csv(as.data.frame(wet_permanova) %>% rownames_to_column("Term"),
          file.path(out_dir, "fticr_pa_permanova_wet_only.csv"))


# Create Plot Labels ----

## Extract treatment effect statistics
extract_treatment_stats <- function(permanova_result, moisture_name) {
  as.data.frame(permanova_result) %>%
    rownames_to_column("Term") %>%
    filter(Term == "Treatment") %>%
    mutate(Moisture = moisture_name)
}

permanova_stats <- bind_rows(
  extract_treatment_stats(dry_permanova, "Dry"),
  extract_treatment_stats(wet_permanova, "Wet")
) %>%
  mutate(
    p_text = case_when(
      `Pr(>F)` < 0.001 ~ "< 0.001",
      TRUE ~ paste0("= ", sprintf("%.3f", `Pr(>F)`))
    )
  )

## Create text labels for plot
plot_labels <- tibble(
  x = min(site_scores$NMDS1) + 0.03 * diff(range(site_scores$NMDS1)),
  y = max(site_scores$NMDS2) - c(0.06, 0.12, 0.18) * diff(range(site_scores$NMDS2)),
  label = c(
    sprintf("NMDS stress = %.2f", nmds$stress),
    sprintf("Dry PERMANOVA: R² = %.2f, p %s", 
            permanova_stats$R2[1], permanova_stats$p_text[1]),
    sprintf("Wet PERMANOVA: R² = %.2f, p %s", 
            permanova_stats$R2[2], permanova_stats$p_text[2])
  )
)


# Create NMDS Plot ----

nmds_plot <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2, 
                                     color = Treatment, shape = Moisture)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text(data = plot_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4, color = "black") +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  scale_shape_manual(values = moisture_shapes, drop = FALSE) +
  labs(x = "NMDS1", y = "NMDS2", color = "Treatment", shape = "Sediment moisture") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right", aspect.ratio = 1)

print(nmds_plot)

# Save plot
ggsave(file.path(out_dir, "Figure3_fticr_pa_nmds_by_inundation_history.png"),
       nmds_plot, width = 6.5, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "Figure3_fticr_pa_nmds_by_inundation_history.pdf"),
       nmds_plot, width = 6.5, height = 6.5)