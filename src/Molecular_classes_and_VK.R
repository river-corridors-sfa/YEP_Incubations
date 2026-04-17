# ==========================================================
# Clear environment
# ==========================================================
rm(list = ls())
gc()

# ==========================================================
# Packages
# ==========================================================
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ==========================================================
# 1. Read both files
# ==========================================================
om_data <- read_csv(
  "C:/Users/gara009/OneDrive - PNNL/RC-SFA - Documents/Study_YEP/FTICR/05_PublishReadyData/Processed Data/YEP_Updated_SN12-Processed_Data.csv"
)

mol_data <- read_csv(
  "C:/Users/gara009/OneDrive - PNNL/RC-SFA - Documents/Study_YEP/FTICR/05_PublishReadyData/Processed Data/YEP_Updated_SN12-Processed_Mol.csv"
)

# ==========================================================
# 2. Standardise merge key
# ==========================================================
om_data  <- om_data  %>% rename(Calibrated_mz = `Calibrated m/z`)
mol_data <- mol_data %>% rename(Calibrated_mz = `Calibrated m/z`)

# ==========================================================
# 3. Clean sample column names in intensity file
# ==========================================================
property_col  <- "Calibrated_mz"
sample_cols   <- setdiff(colnames(om_data), property_col)

clean_names <- sample_cols %>%
  str_remove("_p075\\.corems$") %>%
  str_remove("_p2\\.corems$") %>%
  str_remove("\\.corems$")

colnames(om_data) <- c(property_col, clean_names)

# Keep only sediment-incubated samples
sediment_samples <- clean_names[
  str_detect(clean_names, "^YEP[0-9]+[A-Z]_SIR-[HSU][0-9]+$")
]

om_data <- om_data %>% select(Calibrated_mz, all_of(sediment_samples))

# ==========================================================
# 4. Rename key molecular property columns
# ==========================================================
mol_data <- mol_data %>%
  rename(
    AImod = AI_Mod,
    DBE   = DBE_1,
    HC    = HtoC_ratio,
    OC    = OtoC_ratio
  )

cat("Molecular properties file:", nrow(mol_data), "formulas\n")
cat("Intensity file:", nrow(om_data), "masses x", length(sediment_samples), "samples\n")

# ==========================================================
# 5. Merge by rounded Calibrated_mz
# ==========================================================
mol_data <- mol_data %>% mutate(Calibrated_mz = round(Calibrated_mz, 6))
om_data  <- om_data  %>% mutate(Calibrated_mz = round(Calibrated_mz, 6))

merged <- mol_data %>% inner_join(om_data, by = "Calibrated_mz")

cat("Matched peaks after merge:", nrow(merged), "\n")

# ==========================================================
# 6. Assign Van Krevelen classes (Alan's classification)
# ==========================================================
merged <- merged %>%
  mutate(
    Class = case_when(
      AImod >  0.66 & HC <  1.5                  ~ "Condensed Aromatic",
      AImod <= 0.66 & AImod > 0.5 & HC < 1.5     ~ "Aromatic",
      AImod <= 0.5  & HC <  1.5                   ~ "Unsaturated/Lignin",
      HC >= 2.0     & OC >= 0.9                   ~ "Carbohydrate",
      HC >= 2.0     & OC <  0.9                   ~ "Lipid",
      HC <  2.0     & HC >= 1.5 & N == 0          ~ "Aliphatic",
      HC <  2.0     & HC >= 1.5 & N >  0          ~ "Aliphatic_N",
      TRUE                                         ~ "Unassigned"
    )
  )

# ==========================================================
# 7. Pivot to long: one row per peak x sample (present only)
# ==========================================================
long_df <- merged %>%
  pivot_longer(
    cols      = all_of(sediment_samples),
    names_to  = "Sample_Name",
    values_to = "Intensity"
  ) %>%
  filter(!is.na(Intensity), Intensity > 0)

cat("Total peak-sample observations (present only):", nrow(long_df), "\n")

# ==========================================================
# 8. Add treatment / moisture metadata
# ==========================================================
long_df <- long_df %>%
  mutate(
    Treatment_code = str_extract(Sample_Name, "(?<=_SIR-)[HSU]"),
    Treatment = factor(case_when(
      Treatment_code == "S" ~ "Control",
      Treatment_code == "U" ~ "Unburned + DOC",
      Treatment_code == "H" ~ "High Burn + DOC"
    ), levels = c("Control", "Unburned + DOC", "High Burn + DOC")),
    Moisture = factor(case_when(
      str_detect(Sample_Name, "^YEP1") ~ "Dry",
      str_detect(Sample_Name, "^YEP2") ~ "Wet"
    ), levels = c("Dry", "Wet"))
  )

# ==========================================================
# 9. Per-sample mean & median of AImod, DBE, NOSC
# ==========================================================
sample_summaries <- long_df %>%
  group_by(Sample_Name, Treatment, Moisture) %>%
  summarise(
    mean_AImod   = mean(AImod, na.rm = TRUE),
    median_AImod = median(AImod, na.rm = TRUE),
    mean_DBE     = mean(DBE, na.rm = TRUE),
    median_DBE   = median(DBE, na.rm = TRUE),
    mean_NOSC    = mean(NOSC, na.rm = TRUE),
    median_NOSC  = median(NOSC, na.rm = TRUE),
    n_peaks      = n(),
    .groups      = "drop"
  )

cat("\nPer-sample summary (first rows):\n")
print(head(sample_summaries, 10))

# ==========================================================
# 10. Pivot summaries for faceted plotting
# ==========================================================
summaries_long <- sample_summaries %>%
  pivot_longer(
    cols      = mean_AImod:median_NOSC,
    names_to  = c("Statistic", "Property"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    Statistic = str_to_title(Statistic),
    Property  = toupper(Property)
  )

# ==========================================================
# 11. VK class counts and % relative abundance per sample
# ==========================================================
class_counts <- long_df %>%
  count(Sample_Name, Treatment, Moisture, Class, name = "n_peaks") %>%
  group_by(Sample_Name) %>%
  mutate(rel_abund = (n_peaks / sum(n_peaks)) * 100) %>%
  ungroup() %>%
  mutate(Class = factor(Class, levels = c(
    "Condensed Aromatic", "Aromatic", "Unsaturated/Lignin",
    "Aliphatic", "Aliphatic_N", "Lipid", "Carbohydrate", "Unassigned"
  )))

class_data <- class_counts %>% filter(Class != "Unassigned")

# ==========================================================
# 12. Palettes and comparison lists
# ==========================================================
treatment_colors <- c(
  "Control"          = "#0072B2",
  "Unburned + DOC"   = "darkgreen",
  "High Burn + DOC"  = "#8B4513"
)

moisture_colors <- c(
  "Dry" = "#E69F00",
  "Wet" = "#56B4E9"
)

class_palette <- c(
  "Condensed Aromatic"  = "#1b0e3d",
  "Aromatic"            = "#e41a1c",
  "Unsaturated/Lignin"  = "#ff7f00",
  "Aliphatic"           = "#4daf4a",
  "Aliphatic_N"         = "#a6d854",
  "Lipid"               = "#984ea3",
  "Carbohydrate"        = "#377eb8",
  "Unassigned"          = "grey70"
)

treatment_pairs <- list(
  c("Control", "Unburned + DOC"),
  c("Control", "High Burn + DOC"),
  c("Unburned + DOC", "High Burn + DOC")
)

# ==========================================================
# PART A — Molecular property box plots
# ==========================================================

# ----------------------------------------------------------
# A1. Dry vs Wet WITHIN each treatment
#     x = Treatment, dodged by Moisture
#     Rows = Mean/Median, Cols = AImod/DBE/NOSC
# ----------------------------------------------------------
p_A1 <- ggplot(
  summaries_long,
  aes(x = Treatment, y = Value, fill = Moisture)
) +
  geom_boxplot(
    position = position_dodge(0.75), alpha = 0.7, outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1),
    size = 2.2, alpha = 0.7, show.legend = FALSE
  ) +
  stat_compare_means(
    aes(group = Moisture), method = "wilcox.test",
    label = "p.format", vjust = -0.5, size = 3.2, bracket.size = 0.4
  ) +
  facet_wrap(Statistic ~ Property, scales = "free_y", ncol = 3)+
  scale_fill_manual(values = moisture_colors) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL, y = "Value", fill = "Moisture"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_A1)

# ----------------------------------------------------------
# A2. Treatment comparisons WITHIN each moisture level
#     Separate plots for Mean and Median
#     x = Treatment with pairwise brackets
#     Rows = Moisture, Cols = Property
# ----------------------------------------------------------
for (stat_name in c("Mean", "Median")) {
  
  plot_data <- summaries_long %>% filter(Statistic == stat_name)
  
  p <- ggplot(plot_data, aes(x = Treatment, y = Value, fill = Treatment)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2.2, alpha = 0.7, show.legend = FALSE) +
    stat_compare_means(
      comparisons = treatment_pairs, method = "wilcox.test",
      label = "p.format", size = 3, step.increase = 0.12
    ) +
    facet_grid(Moisture ~ Property, scales = "free_y") +
    scale_fill_manual(values = treatment_colors) +
    theme_bw(base_size = 13) +
    labs(
      title    = paste0(stat_name, " molecular properties: treatment comparisons"),
      subtitle = "Pairwise Wilcoxon rank-sum p-values shown",
      x = NULL, y = paste(stat_name, "value")
    ) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  assign(paste0("p_A2_", tolower(stat_name)), p)
  print(p)
}

# ==========================================================
# PART B — VK class % relative abundance
# ==========================================================

# ----------------------------------------------------------
# B1. Dry vs Wet WITHIN each treatment — faceted by Class
#     free_y so small classes are visible
# ----------------------------------------------------------
p_B1 <- ggplot(
  class_data,
  aes(x = Treatment, y = rel_abund, fill = Moisture)
) +
  geom_boxplot(
    position = position_dodge(0.75), alpha = 0.7, outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1),
    size = 2, alpha = 0.7, show.legend = FALSE
  ) +
  stat_compare_means(
    aes(group = Moisture), method = "wilcox.test",
    label = "p.format", vjust = -0.5, size = 2.8
  ) +
  facet_wrap(~ Class, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = moisture_colors) +
  theme_bw(base_size = 12) +
  labs(
    title    = "VK class relative abundance: Dry vs Wet by treatment",
    subtitle = "Wilcoxon rank-sum p-values shown (note: y-axes differ by class)",
    x = NULL, y = "Relative abundance (%)", fill = "Moisture"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

print(p_B1)

# ----------------------------------------------------------
# B2. Treatment comparisons WITHIN each moisture level
#     faceted by Moisture (rows) x Class (cols)
#     free_y so small classes are visible
# ----------------------------------------------------------
p_B2 <- ggplot(
  class_data,
  aes(x = Treatment, y = rel_abund, fill = Treatment)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, show.legend = FALSE) +
  stat_compare_means(
    comparisons = treatment_pairs, method = "wilcox.test",
    label = "p.signif", size = 3, step.increase = 0.14
  ) +
  facet_grid(Moisture ~ Class, scales = "free_y") +
  scale_fill_manual(values = treatment_colors) +
  theme_bw(base_size = 11) +
  labs(
    title    = "VK class relative abundance: treatment comparisons by moisture",
    subtitle = "Pairwise Wilcoxon: ns p>0.05, * p<0.05, ** p<0.01, *** p<0.001 (y-axes differ by class)",
    x = NULL, y = "Relative abundance (%)"
  ) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))

print(p_B2)

# ----------------------------------------------------------
# B3. Stacked bar — mean RA, Dry vs Wet side by side
#     fixed y (all bars sum to ~100%)
# ----------------------------------------------------------
class_means <- class_counts %>%
  group_by(Treatment, Moisture, Class) %>%
  summarise(
    mean_rel_abund = mean(rel_abund, na.rm = TRUE),
    .groups = "drop"
  )

p_B3 <- ggplot(
  class_means,
  aes(x = Moisture, y = mean_rel_abund, fill = Class)
) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~ Treatment, nrow = 1) +
  scale_fill_manual(values = class_palette) +
  theme_bw(base_size = 13) +
  labs(
    title = "Mean VK class composition: Dry vs Wet within each treatment",
    x = NULL, y = "Mean relative abundance (%)", fill = "VK Class"
  )

print(p_B3)

# ==========================================================
# PART C — Wilcoxon test tables
# ==========================================================

# ----------------------------------------------------------
# C1. Molecular properties: Dry vs Wet within each treatment
# ----------------------------------------------------------
wilcox_moisture_props <- summaries_long %>%
  group_by(Statistic, Property, Treatment) %>%
  wilcox_test(Value ~ Moisture) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

cat("\n===== Dry vs Wet within treatment (molecular properties) =====\n")
print(wilcox_moisture_props, n = Inf)

# ----------------------------------------------------------
# C2. Molecular properties: pairwise treatments within moisture
# ----------------------------------------------------------
wilcox_treat_props <- summaries_long %>%
  group_by(Statistic, Property, Moisture) %>%
  wilcox_test(Value ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

cat("\n===== Pairwise treatment comparisons within moisture (molecular properties) =====\n")
print(wilcox_treat_props, n = Inf)

# ----------------------------------------------------------
# C3. VK classes: Dry vs Wet within each treatment
# ----------------------------------------------------------
wilcox_moisture_class <- class_data %>%
  group_by(Class, Treatment) %>%
  wilcox_test(rel_abund ~ Moisture) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

cat("\n===== Dry vs Wet within treatment (VK classes) =====\n")
print(wilcox_moisture_class, n = Inf)

# ----------------------------------------------------------
# C4. VK classes: pairwise treatments within moisture
# ----------------------------------------------------------
wilcox_treat_class <- class_data %>%
  group_by(Class, Moisture) %>%
  wilcox_test(rel_abund ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

cat("\n===== Pairwise treatment comparisons within moisture (VK classes) =====\n")
print(wilcox_treat_class, n = Inf)

# ==========================================================
# SAVE FIGURES
# ==========================================================
ggsave("A1_mol_props_dry_vs_wet.png",            p_A1,            width = 12, height = 8,  dpi = 300)
ggsave("A2_mol_props_treatment_comp_mean.png",   p_A2_mean,       width = 12, height = 8,  dpi = 300)
ggsave("A2_mol_props_treatment_comp_median.png", p_A2_median,     width = 12, height = 8,  dpi = 300)
ggsave("B1_VK_class_dry_vs_wet.png",             p_B1,            width = 14, height = 10, dpi = 300)
ggsave("B2_VK_class_treatment_comp.png",         p_B2,            width = 16, height = 8,  dpi = 300)
ggsave("B3_VK_class_stacked_bar.png",            p_B3,            width = 10, height = 6,  dpi = 300)

# ==========================================================
# SAVE TABLES
# ==========================================================
write_csv(sample_summaries,      "per_sample_mol_property_summaries.csv")
write_csv(class_counts,          "per_sample_VK_class_counts.csv")
write_csv(wilcox_moisture_props, "wilcox_dry_vs_wet_mol_properties.csv")
write_csv(wilcox_treat_props,    "wilcox_treatment_pairwise_mol_properties.csv")
write_csv(wilcox_moisture_class, "wilcox_dry_vs_wet_VK_classes.csv")
write_csv(wilcox_treat_class,    "wilcox_treatment_pairwise_VK_classes.csv")

cat("\n===== Done — all plots and tables saved =====\n")

# ==========================================================
# PART B1b — Standard VK classes from mol file
# ==========================================================

# ----------------------------------------------------------
# Check what classes exist in each scheme
# ----------------------------------------------------------
cat("\n--- bs1_class levels ---\n")
long_df %>% count(bs1_class, sort = TRUE) %>% print(n = Inf)

cat("\n--- bs2_class levels ---\n")
long_df %>% count(bs2_class, sort = TRUE) %>% print(n = Inf)

cat("\n--- bs3_class levels ---\n")
long_df %>% count(bs3_class, sort = TRUE) %>% print(n = Inf)

# ----------------------------------------------------------
# Function to build the B1-style plot for any class column
# ----------------------------------------------------------
make_vk_class_plot <- function(df, class_col, plot_title) {
  
  # Per-sample class counts and % RA
  class_ct <- df %>%
    rename(VK_Class = !!sym(class_col)) %>%
    filter(!is.na(VK_Class), VK_Class != "") %>%
    count(Sample_Name, Treatment, Moisture, VK_Class, name = "n_peaks") %>%
    group_by(Sample_Name) %>%
    mutate(rel_abund = (n_peaks / sum(n_peaks)) * 100) %>%
    ungroup()
  
  # Ordered factor by median abundance (largest first)
  class_order <- class_ct %>%
    group_by(VK_Class) %>%
    summarise(med = median(rel_abund), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(VK_Class)
  
  class_ct$VK_Class <- factor(class_ct$VK_Class, levels = class_order)
  
  # Plot
  p <- ggplot(
    class_ct,
    aes(x = Treatment, y = rel_abund, fill = Moisture)
  ) +
    geom_boxplot(
      position = position_dodge(0.75), alpha = 0.7, outlier.shape = NA
    ) +
    geom_point(
      position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1),
      size = 2, alpha = 0.7, show.legend = FALSE
    ) +
    stat_compare_means(
      aes(group = Moisture), method = "wilcox.test",
      label = "p.format", vjust = -0.5, size = 2.8
    ) +
    facet_wrap(~ VK_Class, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = moisture_colors) +
    theme_bw(base_size = 12) +
    labs(
      title    = plot_title,
      subtitle = "Wilcoxon rank-sum p-values shown (y-axes differ by class)",
      x = NULL, y = "Relative abundance (%)", fill = "Moisture"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  
  # Wilcoxon table
  wilcox_tbl <- class_ct %>%
    group_by(VK_Class, Treatment) %>%
    wilcox_test(rel_abund ~ Moisture) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  return(list(plot = p, class_counts = class_ct, wilcox = wilcox_tbl))
}

# ----------------------------------------------------------
# Build plots for each classification scheme
# ----------------------------------------------------------
res_bs1 <- make_vk_class_plot(
  long_df, "bs1_class",
  "Standard VK classes (bs1): Dry vs Wet by treatment"
)

res_bs2 <- make_vk_class_plot(
  long_df, "bs2_class",
  "Standard VK classes (bs2): Dry vs Wet by treatment"
)

res_bs3 <- make_vk_class_plot(
  long_df, "bs3_class",
  "Standard VK classes (bs3): Dry vs Wet by treatment"
)

# ----------------------------------------------------------
# Print plots
# ----------------------------------------------------------
print(res_bs1$plot)
print(res_bs2$plot)
print(res_bs3$plot)

# ----------------------------------------------------------
# Print Wilcoxon tables
# ----------------------------------------------------------
cat("\n===== Dry vs Wet within treatment (bs1_class) =====\n")
print(res_bs1$wilcox, n = Inf)

cat("\n===== Dry vs Wet within treatment (bs2_class) =====\n")
print(res_bs2$wilcox, n = Inf)

cat("\n===== Dry vs Wet within treatment (bs3_class) =====\n")
print(res_bs3$wilcox, n = Inf)

# ----------------------------------------------------------
# Save plots
# ----------------------------------------------------------
ggsave("B1b_VK_bs1_class_dry_vs_wet.png", res_bs1$plot, width = 14, height = 10, dpi = 300)
ggsave("B1b_VK_bs2_class_dry_vs_wet.png", res_bs2$plot, width = 14, height = 10, dpi = 300)
ggsave("B1b_VK_bs3_class_dry_vs_wet.png", res_bs3$plot, width = 14, height = 10, dpi = 300)

# ----------------------------------------------------------
# Save tables
# ----------------------------------------------------------
write_csv(res_bs1$wilcox,        "wilcox_dry_vs_wet_bs1_class.csv")
write_csv(res_bs2$wilcox,        "wilcox_dry_vs_wet_bs2_class.csv")
write_csv(res_bs3$wilcox,        "wilcox_dry_vs_wet_bs3_class.csv")
write_csv(res_bs1$class_counts,  "per_sample_bs1_class_counts.csv")
write_csv(res_bs2$class_counts,  "per_sample_bs2_class_counts.csv")
write_csv(res_bs3$class_counts,  "per_sample_bs3_class_counts.csv")