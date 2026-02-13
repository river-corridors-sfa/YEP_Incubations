rm(list=ls(all=T))
# Load required libraries
library(ggplot2)
library(dplyr)
library(FSA)  # for Dunn's test
library(rcompanion)  # alternative for Dunn's test
library(ggpubr)  # for adding stats to plots

# Read your data
data <- read.csv("modeling_outputs/YEP_Complete_Analysis_Following_Paper_Enhanced.csv")

# Filter out the "Unknown" treatment types and create a sediment condition variable
data_filtered <- data %>%
  filter(Treatment_Type != "Unknown") %>%
  mutate(
    Sediment_Condition = ifelse(grepl("^Wet_", Treatment_Type), "Wet", "Dry"),
    Treatment = case_when(
      grepl("HighBurn", Treatment_Type) ~ "High Burn + DOC",
      grepl("Unburned", Treatment_Type) ~ "Unburned + DOC",
      grepl("Control", Treatment_Type) ~ "Control",
      TRUE ~ "Other"
    )
  )

# Split data by sediment condition
data_wet <- data_filtered %>% filter(Sediment_Condition == "Wet")
data_dry <- data_filtered %>% filter(Sediment_Condition == "Dry")

# === VMAX ANALYSIS ===
cat("\n========== VMAX ANALYSIS ==========\n")

# Wet sediments - Vmax
cat("\n--- WET SEDIMENTS - Vmax ---\n")
kw_vmax_wet <- kruskal.test(Combined_Vmax_per_h ~ Treatment, data = data_wet)
print(kw_vmax_wet)

if(kw_vmax_wet$p.value < 0.05) {
  dunn_vmax_wet <- dunnTest(Combined_Vmax_per_h ~ Treatment, 
                            data = data_wet, 
                            method = "bonferroni")
  print(dunn_vmax_wet)
} else {
  cat("No significant differences found, skipping post-hoc test\n")
}

# Dry sediments - Vmax
cat("\n--- DRY SEDIMENTS - Vmax ---\n")
kw_vmax_dry <- kruskal.test(Combined_Vmax_per_h ~ Treatment, data = data_dry)
print(kw_vmax_dry)

if(kw_vmax_dry$p.value < 0.05) {
  dunn_vmax_dry <- dunnTest(Combined_Vmax_per_h ~ Treatment, 
                            data = data_dry, 
                            method = "bonferroni")
  print(dunn_vmax_dry)
} else {
  cat("No significant differences found, skipping post-hoc test\n")
}

# === kL ANALYSIS ===
cat("\n========== kL ANALYSIS ==========\n")

# Wet sediments - kL
cat("\n--- WET SEDIMENTS - kL ---\n")
kw_kl_wet <- kruskal.test(Combined_kL_per_h ~ Treatment, data = data_wet)
print(kw_kl_wet)

if(kw_kl_wet$p.value < 0.05) {
  dunn_kl_wet <- dunnTest(Combined_kL_per_h ~ Treatment, 
                          data = data_wet, 
                          method = "bonferroni")
  print(dunn_kl_wet)
} else {
  cat("No significant differences found, skipping post-hoc test\n")
}

# Dry sediments - kL
cat("\n--- DRY SEDIMENTS - kL ---\n")
kw_kl_dry <- kruskal.test(Combined_kL_per_h ~ Treatment, data = data_dry)
print(kw_kl_dry)

if(kw_kl_dry$p.value < 0.05) {
  dunn_kl_dry <- dunnTest(Combined_kL_per_h ~ Treatment, 
                          data = data_dry, 
                          method = "bonferroni")
  print(dunn_kl_dry)
} else {
  cat("No significant differences found, skipping post-hoc test\n")
}

# === CREATE PLOTS WITH STATS ===

# Vmax plot with separate stats for wet and dry
plot_vmax <- ggplot(data_filtered, aes(x = Sediment_Condition, y = Combined_Vmax_per_h, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) +
  stat_compare_means(aes(group = Treatment), 
                     method = "kruskal.test",
                     label = "p.format",
                     label.y = max(data_filtered$Combined_Vmax_per_h, na.rm = TRUE) * 1.05) +
  labs(
    title = 'A',
    x = "Sediment Condition",
    y = "Combined Vmax (per h)",
    fill = "Treatment") +
  theme_bw() 

# kL plot with separate stats for wet and dry
plot_kl <- ggplot(data_filtered, aes(x = Sediment_Condition, y = Combined_kL_per_h, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) +
  stat_compare_means(aes(group = Treatment), 
                     method = "kruskal.test",
                     label = "p.format",
                     label.y = max(data_filtered$Combined_kL_per_h, na.rm = TRUE) * 1.05) +
  labs(
    title ="B",
    x = "Sediment Condition",
    y = "Combined kL (per h)",
    fill = "Treatment") +
  theme_bw() 

# Display plots
print(plot_vmax)
print(plot_kl)

# Save plots
ggsave("Combined_Vmax_Plot_with_Stats.png", plot_vmax, width = 10, height = 7, dpi = 300)
ggsave("Combined_kL_Plot_with_Stats.png", plot_kl, width = 10, height = 7, dpi = 300)

library(ggplot2)
library(ggpubr)
library(patchwork)

# Define custom colors for treatments
treatment_colors <- c(
  "Control" = "#0072B2",              # Blue
  "Unburned + DOC" = "darkgreen",       # Green
  "High Burn + DOC" = "#8B4513"       # Brown (saddle brown)
)

# Alternative if you want a darker/blacker brown:
# "High Burn + DOC" = "#3E2723"  # Very dark brown (almost black)
# "High Burn + DOC" = "#1A1A1A"  # Nearly black

# Create Vmax plot with custom colors
plot_vmax <- ggplot(data_filtered, aes(x = Sediment_Condition, y = Combined_Vmax_per_h, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) +
  stat_compare_means(aes(group = Treatment),
                     method = "kruskal.test",
                     label = "p.format",
                     label.y = max(data_filtered$Combined_Vmax_per_h, na.rm = TRUE) * 1.05) +
  scale_fill_manual(values = treatment_colors) +  # Apply custom colors
  labs(
    title = 'A',
    x = "Sediment Condition",
    y = "Combined Vmax (per h)",
    fill = "Treatment") +
  theme_bw() +
  theme(legend.position = "none")

# Create kL plot with custom colors
plot_kl <- ggplot(data_filtered, aes(x = Sediment_Condition, y = Combined_kL_per_h, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) +
  stat_compare_means(aes(group = Treatment),
                     method = "kruskal.test",
                     label = "p.format",
                     label.y = max(data_filtered$Combined_kL_per_h, na.rm = TRUE) * 1.05) +
  scale_fill_manual(values = treatment_colors) +  # Apply custom colors
  labs(
    title = "B",
    x = "Sediment Condition",
    y = "Combined kL (per h)",
    fill = "Treatment") +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots with shared legend at bottom
combined_plot <- plot_vmax + plot_kl + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

# Display combined plot
print(combined_plot)

# Save combined plot
ggsave("Combined_Vmax_kL_Plot.png", combined_plot, width = 12, height = 6, dpi = 300)
