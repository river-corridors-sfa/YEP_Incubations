rm(list=ls())

# Load required libraries
library(tidyverse)
library(ggplot2)
library(car)
library(lme4)
library(effectsize)
library(broom)
library(gridExtra)
library(emmeans)

# ====== Read the data files =====
respiration_data <- read.csv("Data/segmented_respiration_analysis.csv")
npoc_data <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_NPOC_TN.csv", skip = 2) %>%
  filter(grepl('YEP',Sample_Name))
co2_data <- read.csv('data/co2_production_rates.csv')
mass_data <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Water_Mass_Volume.csv", skip = 2)%>%
  filter(grepl('YEP',Sample_Name))

# ====== Data cleaning =======
# Function to extract treatment information from sample names
extract_treatments <- function(sample_name) {
  if (str_detect(sample_name, "YEP[12][ABC]_")) {
    site <- str_extract(sample_name, "YEP[12]")
    replicate <- str_extract(sample_name, "(?<=YEP[12])[ABC]")
    condition <- str_extract(sample_name, "(?<=-)[HSU]")
  } else if (str_detect(sample_name, "-W")) {
    site <- "Water"
    replicate <- str_extract(sample_name, "\\d+$")
    condition <- str_extract(sample_name, "(?<=W)[HSU]")
  } else {
    site <- NA
    replicate <- NA
    condition <- NA
  }
  
  return(data.frame(
    Sample_Name = sample_name,
    Site = site,
    Replicate = replicate,
    Condition = condition
  ))
}

# Create treatment lookup table for respiration data
respiration_treatments <- map_dfr(respiration_data$Sample_Name, extract_treatments) %>%
  filter(!is.na(Site), !is.na(Condition), Site != "Water") %>%
  mutate(
    Sediment_Type = case_when(
      Site == "YEP1" ~ "Dry Sediments",
      Site == "YEP2" ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    DOC_Treatment = case_when(
      Condition == "H" ~ "PyOM Added (High Burn)",
      Condition == "U" ~ "DOM Added (Unburned)", 
      Condition == "S" ~ "Synthetic Water (Control)",
      TRUE ~ NA_character_
    ),
    DOC_Added = case_when(
      Condition %in% c("H", "U") ~ "DOC Added",
      Condition == "S" ~ "No DOC (Control)",
      TRUE ~ NA_character_
    )
  )

# Clean respiration data using left_join
respiration_clean <- respiration_data %>%
  left_join(respiration_treatments, by = "Sample_Name") %>%
  filter(!is.na(Site), !is.na(Condition)) %>%
  mutate(
    rate_negative = Rate_mg_L_h_break_1,
    Treatment_Combo = paste(Sediment_Type, DOC_Treatment, sep = " + ")
  )

# Clean NPOC data
npoc_clean <- npoc_data %>%
  filter(!is.na(Sample_Name), Sample_Name != "", !str_detect(Sample_Name, "^#")) %>%
  dplyr::select(Sample_Name, Extractable_NPOC_mg_per_kg) %>%
  filter(!is.na(Extractable_NPOC_mg_per_kg), 
         Extractable_NPOC_mg_per_kg != "-9999") %>%
  mutate(Extractable_NPOC_mg_per_kg = as.numeric(Extractable_NPOC_mg_per_kg))

# Create treatment lookup for NPOC data
npoc_treatments <- map_dfr(npoc_clean$Sample_Name, extract_treatments) %>%
  filter(!is.na(Site), Site != "Water") %>%
  mutate(
    Sediment_Type = case_when(
      Site == "YEP1" ~ "Dry Sediments",
      Site == "YEP2" ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    DOC_Treatment = case_when(
      Condition == "H" ~ "PyOM Added (High Burn)",
      Condition == "U" ~ "DOM Added (Unburned)", 
      Condition == "S" ~ "Synthetic Water (Control)",
      TRUE ~ NA_character_
    ),
    DOC_Added = case_when(
      Condition %in% c("H", "U") ~ "DOC Added",
      Condition == "S" ~ "No DOC (Control)",
      TRUE ~ NA_character_
    )
  )

# Join NPOC data with treatments
npoc_final <- npoc_clean %>%
  left_join(npoc_treatments, by = "Sample_Name") %>%
  filter(!is.na(Condition))

# Clean CO2 data
co2_clean <- co2_data %>%
  filter(!is.na(Sample_Name), Sample_Name != "", !str_detect(Sample_Name, "^#")) %>%
  dplyr::select(Sample_Name, CO2_Production_Rate_mol_per_L_per_H) %>%
  mutate(CO2_Production_Rate_mol_per_L_per_H= as.numeric(CO2_Production_Rate_mol_per_L_per_H))

# Create treatment lookup for CO2 data
co2_treatments <- map_dfr(co2_clean$Sample_Name, extract_treatments) %>%
  filter(!is.na(Site), Site != "Water") %>%
  mutate(
    Sediment_Type = case_when(
      Site == "YEP1" ~ "Dry Sediments",
      Site == "YEP2" ~ "Wet Sediments",
      TRUE ~ NA_character_
    ),
    DOC_Treatment = case_when(
      Condition == "H" ~ "PyOM Added (High Burn)",
      Condition == "U" ~ "DOM Added (Unburned)", 
      Condition == "S" ~ "Synthetic Water (Control)",
      TRUE ~ NA_character_
    ),
    DOC_Added = case_when(
      Condition %in% c("H", "U") ~ "DOC Added",
      Condition == "S" ~ "No DOC (Control)",
      TRUE ~ NA_character_
    )
  )

# Join CO2 data with treatments
co2_final <- co2_clean %>%
  left_join(co2_treatments, by = "Sample_Name") %>%
  filter(!is.na(Site))

# Clean mass data for mixed effects model
mass_clean <- mass_data %>%
  filter(!is.na(Sample_Name), Sample_Name != "", !str_detect(Sample_Name, "^#")) %>%
  dplyr::select(Sample_Name, Water_Mass_g, Dry_Sediment_Mass_g) %>%
  filter(!is.na(Water_Mass_g), !is.na(Dry_Sediment_Mass_g),
         Water_Mass_g != "-9999", Dry_Sediment_Mass_g != "-9999") %>%
  mutate(
    Water_Mass_g = as.numeric(Water_Mass_g),
    Dry_Sediment_Mass_g = as.numeric(Dry_Sediment_Mass_g)
  )
# ======== Plots ======
#Boxplots of respiration rates per dry mass
# Normalize per dry mass
respiration_clean = respiration_clean %>%
  left_join(mass_clean ,
            by = 'Sample_Name') %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(is.na(Dry_Sediment_Mass_g), -9999,
                                                      ((rate_negative*Water_Mass_g*0.001)/(Dry_Sediment_Mass_g*0.001)))) %>%
  filter(!is.na(rate_negative)) 

fig2 <- ggplot(respiration_clean, aes(x = DOC_Treatment, y = Respiration_Rate_mg_DO_per_kg_per_H, fill = Sediment_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             alpha = 0.6) +
  labs(
    x = " ",
    y = "Respiration Rate (mg O₂ kg⁻¹ h⁻¹)",
    fill = "Sediment Type"
  ) +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Dry Sediments" = "lightblue", "Wet Sediments" = "darkblue"))

print(fig2)

fig3 <- ggplot(respiration_clean, aes(fill = DOC_Treatment, y = Respiration_Rate_mg_DO_per_kg_per_H, x = Sediment_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             alpha = 0.6) +
  labs(
    x = " ",
    y = "Respiration Rate (mg O₂ kg⁻¹ h⁻¹)",
    fill = "Sediment Type"
  ) +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot(respiration_clean, aes(x = DOC_Treatment, y = rate_negative, fill = Sediment_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             alpha = 0.6) +
  labs(
    x = " ",
    y = "Respiration Rate (mg O₂ L⁻¹ h⁻¹)",
    fill = "Sediment Type"
  ) +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Dry Sediments" = "lightblue", "Wet Sediments" = "darkblue"))


# ANOVA for respiration rates
respiration_anova <- aov(Respiration_Rate_mg_DO_per_kg_per_H ~ Sediment_Type * DOC_Treatment, data = respiration_clean)
print("ANOVA Results for Respiration Rates:")
print(summary(respiration_anova))
print("Type III ANOVA:")
print(Anova(respiration_anova, type = "III"))

# Post-hoc tests
print("Post-hoc comparisons:")
respiration_emm <- emmeans(respiration_anova, ~ Sediment_Type * DOC_Treatment)
print(pairs(respiration_emm))

# ADDITIONAL PLOTS with corrected treatments

# NPOC boxplots
fig_npoc <- ggplot(npoc_final, aes(x = DOC_Treatment, y = Extractable_NPOC_mg_per_kg, fill = Sediment_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             alpha = 0.6) +
  labs(x = "DOC Treatment",
    y = "NPOC (mg kg⁻¹ dry sediment)",
    fill = "Sediment Type"
  ) +
  theme_bw() +
  theme(legend.position = 'bottom',  axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Dry Sediments" = "lightgreen", "Wet Sediments" = "darkgreen"))

print(fig_npoc)

# CO2 boxplots
fig_co2 <- ggplot(co2_final, aes(x = DOC_Treatment, y = CO2_Production_Rate_mol_per_L_per_H, fill = Sediment_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             alpha = 0.6) +
  labs( x = "DOC Treatment",
    y = "CO₂Production (mol L⁻¹)",
    fill = "Sediment Type"
  ) +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Dry Sediments" = "coral", "Wet Sediments" = "brown")) +
  scale_y_continuous(labels = scales::scientific)

print(fig_co2)

# EFFECT SIZE CALCULATIONS - Separate by Sediment Type
print("=== EFFECT SIZE ANALYSIS: DOC TREATMENTS vs CONTROL (BY SEDIMENT TYPE) ===")

# Separate data by sediment type first
dry_sediment_data <- respiration_clean %>% filter(Sediment_Type == "Dry Sediments")
wet_sediment_data <- respiration_clean %>% filter(Sediment_Type == "Wet Sediments")

# 1. PyOM (H) vs Synthetic Water Control (S) - BY SEDIMENT TYPE
print("1. PyOM (H) vs Synthetic Water Control (S):")

# Dry sediments: H vs S
dry_h_vs_s_data <- dry_sediment_data %>% filter(Condition %in% c("H", "S"))
if(nrow(dry_h_vs_s_data) > 0) {
  dry_h_vs_s_effect <- cohens_d(Respiration_Rate_mg_DO_per_kg_per_H ~ Condition, data = dry_h_vs_s_data)
  print("   DRY SEDIMENTS - PyOM vs Control - Cohen's d:")
  print(dry_h_vs_s_effect)
} else {
  dry_h_vs_s_effect <- NULL
}

# Wet sediments: H vs S  
wet_h_vs_s_data <- wet_sediment_data %>% filter(Condition %in% c("H", "S"))
if(nrow(wet_h_vs_s_data) > 0) {
  wet_h_vs_s_effect <- cohens_d(Respiration_Rate_mg_DO_per_kg_per_H ~ Condition, data = wet_h_vs_s_data)
  print("   WET SEDIMENTS - PyOM vs Control - Cohen's d:")
  print(wet_h_vs_s_effect)
} else {
  wet_h_vs_s_effect <- NULL
}

# 2. Unburned DOM (U) vs Synthetic Water Control (S) - BY SEDIMENT TYPE
print("2. Unburned DOM (U) vs Synthetic Water Control (S):")

# Dry sediments: U vs S
dry_u_vs_s_data <- dry_sediment_data %>% filter(Condition %in% c("U", "S"))
if(nrow(dry_u_vs_s_data) > 0) {
  dry_u_vs_s_effect <- cohens_d(Respiration_Rate_mg_DO_per_kg_per_H ~ Condition, data = dry_u_vs_s_data)
  print("   DRY SEDIMENTS - Unburned DOM vs Control - Cohen's d:")
  print(dry_u_vs_s_effect)
} else {
  dry_u_vs_s_effect <- NULL
}

# Wet sediments: U vs S
wet_u_vs_s_data <- wet_sediment_data %>% filter(Condition %in% c("U", "S"))
if(nrow(wet_u_vs_s_data) > 0) {
  wet_u_vs_s_effect <- cohens_d(Respiration_Rate_mg_DO_per_kg_per_H ~ Condition, data = wet_u_vs_s_data)
  print("   WET SEDIMENTS - Unburned DOM vs Control - Cohen's d:")
  print(wet_u_vs_s_effect)
} else {
  wet_u_vs_s_effect <- NULL
}

# 3. Combined DOC treatments (H+U) vs Control (S) - BY SEDIMENT TYPE
print("3. Combined DOC Added vs Control:")

if(nrow(dry_sediment_data) > 0) {
  dry_doc_effect <- cohens_d(Respiration_Rate_mg_DO_per_kg_per_H ~ DOC_Added, data = dry_sediment_data)
  print("   DRY SEDIMENTS - Combined DOC vs Control - Cohen's d:")
  print(dry_doc_effect)
}

if(nrow(wet_sediment_data) > 0) {
  wet_doc_effect <- cohens_d(Respiration_Rate_mg_DO_per_kg_per_H ~ DOC_Added, data = wet_sediment_data)
  print("   WET SEDIMENTS - Combined DOC vs Control - Cohen's d:")
  print(wet_doc_effect)
}

# 4. Overall model effect sizes
overall_effects <- eta_squared(respiration_anova)
print("Overall Model Effect Sizes (Eta squared):")
print(overall_effects)

# 5. Create comprehensive effect size summary table
effect_size_summary <- data.frame(
  Comparison = c(
    "PyOM vs Control (Dry Sediments)",
    "PyOM vs Control (Wet Sediments)", 
    "Unburned DOM vs Control (Dry Sediments)",
    "Unburned DOM vs Control (Wet Sediments)",
    "Combined DOC vs Control (Dry Sediments)",
    "Combined DOC vs Control (Wet Sediments)",
    "Sediment Type (Main Effect)",
    "DOC Treatment (Main Effect)",
    "Sediment × DOC (Interaction)"
  ),
  Effect_Size = c(
    if(!is.null(dry_h_vs_s_effect)) abs(dry_h_vs_s_effect$Cohens_d) else NA,
    if(!is.null(wet_h_vs_s_effect)) abs(wet_h_vs_s_effect$Cohens_d) else NA,
    if(!is.null(dry_u_vs_s_effect)) abs(dry_u_vs_s_effect$Cohens_d) else NA,
    if(!is.null(wet_u_vs_s_effect)) abs(wet_u_vs_s_effect$Cohens_d) else NA,
    if(exists("dry_doc_effect")) abs(dry_doc_effect$Cohens_d) else NA,
    if(exists("wet_doc_effect")) abs(wet_doc_effect$Cohens_d) else NA,
    overall_effects$Eta2[1],  # Sediment Type
    overall_effects$Eta2[2],  # DOC Treatment  
    overall_effects$Eta2[3]   # Interaction
  ),
  Effect_Type = c(
    rep("Cohen's d", 6),
    rep("Eta squared", 3)
  )
) %>%
  mutate(
    Interpretation = case_when(
      Effect_Type == "Cohen's d" & Effect_Size < 0.2 ~ "Negligible",
      Effect_Type == "Cohen's d" & Effect_Size < 0.5 ~ "Small", 
      Effect_Type == "Cohen's d" & Effect_Size < 0.8 ~ "Medium",
      Effect_Type == "Cohen's d" & Effect_Size >= 0.8 ~ "Large",
      Effect_Type == "Eta squared" & Effect_Size < 0.01 ~ "Small",
      Effect_Type == "Eta squared" & Effect_Size < 0.06 ~ "Medium", 
      Effect_Type == "Eta squared" & Effect_Size >= 0.06 ~ "Large",
      is.na(Effect_Size) ~ "Unable to calculate",
      TRUE ~ "Unable to calculate"
    )
  )

print("=== COMPREHENSIVE EFFECT SIZE SUMMARY TABLE ===")
print(effect_size_summary)

# 6. Specific contrasts within each sediment type
print("=== PLANNED CONTRASTS BY SEDIMENT TYPE ===")

# Contrasts for dry sediments
print("DRY SEDIMENTS:")
dry_contrast_results <- emmeans(respiration_anova, ~ DOC_Treatment | Sediment_Type)
dry_contrast_pairs <- pairs(dry_contrast_results)
print(dry_contrast_pairs)

# Contrasts for wet sediments  
print("WET SEDIMENTS:")
wet_contrast_results <- emmeans(respiration_anova, ~ DOC_Treatment | Sediment_Type)
wet_contrast_pairs <- pairs(wet_contrast_results)
print(wet_contrast_pairs)

# Effect sizes for specific contrasts by sediment type
contrast_effects <- eff_size(emmeans(respiration_anova, ~ DOC_Treatment | Sediment_Type), 
                             sigma = sigma(respiration_anova), 
                             edf = df.residual(respiration_anova))
print("Effect sizes for contrasts by sediment type:")
print(contrast_effects)

# ===== ABSOLUTE DIFFERENCES ANALYSIS (Original Units) =====
print("=== ABSOLUTE DIFFERENCES IN ORIGINAL UNITS (mg O₂ kg⁻¹ h⁻¹) ===")

# Calculate group means for each sediment type and treatment combination
group_means <- respiration_clean %>%
  group_by(Sediment_Type, DOC_Treatment) %>%
  summarise(
    mean_rate = mean(Respiration_Rate_mg_DO_per_kg_per_H, na.rm = TRUE),
    sd_rate = sd(Respiration_Rate_mg_DO_per_kg_per_H, na.rm = TRUE),
    n = n(),
    se_rate = sd_rate / sqrt(n),
    .groups = 'drop'
  )

print("Group Means:")
print(group_means)

# Calculate absolute differences by sediment type
absolute_differences <- data.frame(
  Sediment_Type = character(),
  Comparison = character(),
  Group1_Mean = numeric(),
  Group2_Mean = numeric(),
  Absolute_Difference = numeric(),
  Percent_Difference = numeric(),
  Combined_SE = numeric(),
  stringsAsFactors = FALSE
)

sediment_types <- unique(respiration_clean$Sediment_Type)

for(sed_type in sediment_types) {
  # Get means for this sediment type
  means_subset <- group_means %>% filter(Sediment_Type == sed_type)
  
  # PyOM vs Control
  pyom_mean <- means_subset$mean_rate[means_subset$DOC_Treatment == "PyOM Added (High Burn)"]
  control_mean <- means_subset$mean_rate[means_subset$DOC_Treatment == "Synthetic Water (Control)"]
  pyom_se <- means_subset$se_rate[means_subset$DOC_Treatment == "PyOM Added (High Burn)"]
  control_se <- means_subset$se_rate[means_subset$DOC_Treatment == "Synthetic Water (Control)"]
  
  if(length(pyom_mean) > 0 && length(control_mean) > 0) {
    abs_diff_pyom <- abs(pyom_mean - control_mean)
    percent_diff_pyom <- abs_diff_pyom / abs(control_mean) * 100
    combined_se_pyom <- sqrt(pyom_se^2 + control_se^2)
    
    absolute_differences <- rbind(absolute_differences, data.frame(
      Sediment_Type = sed_type,
      Comparison = "PyOM vs Control",
      Group1_Mean = pyom_mean,
      Group2_Mean = control_mean,
      Absolute_Difference = abs_diff_pyom,
      Percent_Difference = percent_diff_pyom,
      Combined_SE = combined_se_pyom
    ))
  }
  
  # DOM vs Control
  dom_mean <- means_subset$mean_rate[means_subset$DOC_Treatment == "DOM Added (Unburned)"]
  dom_se <- means_subset$se_rate[means_subset$DOC_Treatment == "DOM Added (Unburned)"]
  
  if(length(dom_mean) > 0 && length(control_mean) > 0) {
    abs_diff_dom <- abs(dom_mean - control_mean)
    percent_diff_dom <- abs_diff_dom / abs(control_mean) * 100
    combined_se_dom <- sqrt(dom_se^2 + control_se^2)
    
    absolute_differences <- rbind(absolute_differences, data.frame(
      Sediment_Type = sed_type,
      Comparison = "DOM vs Control",
      Group1_Mean = dom_mean,
      Group2_Mean = control_mean,
      Absolute_Difference = abs_diff_dom,
      Percent_Difference = percent_diff_dom,
      Combined_SE = combined_se_dom
    ))
  }
  
  # PyOM vs DOM
  if(length(pyom_mean) > 0 && length(dom_mean) > 0) {
    abs_diff_treatments <- abs(pyom_mean - dom_mean)
    percent_diff_treatments <- abs_diff_treatments / abs(dom_mean) * 100
    combined_se_treatments <- sqrt(pyom_se^2 + dom_se^2)
    
    absolute_differences <- rbind(absolute_differences, data.frame(
      Sediment_Type = sed_type,
      Comparison = "PyOM vs DOM",
      Group1_Mean = pyom_mean,
      Group2_Mean = dom_mean,
      Absolute_Difference = abs_diff_treatments,
      Percent_Difference = percent_diff_treatments,
      Combined_SE = combined_se_treatments
    ))
  }
}

# Add practical significance interpretation
absolute_differences <- absolute_differences %>%
  mutate(
    Practical_Significance = case_when(
      Absolute_Difference > 10 ~ "Large (>10 mg O₂ kg⁻¹ h⁻¹)",
      Absolute_Difference > 5 ~ "Medium (5-10 mg O₂ kg⁻¹ h⁻¹)",
      Absolute_Difference > 2 ~ "Small (2-5 mg O₂ kg⁻¹ h⁻¹)",
      TRUE ~ "Negligible (<2 mg O₂ kg⁻¹ h⁻¹)"
    )
  )

print("=== ABSOLUTE DIFFERENCES TABLE ===")
print(absolute_differences)

# Create visualization
library(ggplot2)

p_abs_diff <- ggplot(absolute_differences, aes(x = reorder(paste(Sediment_Type, Comparison), Absolute_Difference), 
                                               y = Absolute_Difference, 
                                               fill = Sediment_Type)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = Absolute_Difference - Combined_SE, 
                    ymax = Absolute_Difference + Combined_SE), 
                width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(title = "Absolute Differences Between Treatments",
       subtitle = "Error bars show combined standard error",
       x = "Comparison",
       y = "Absolute Difference (mg O₂ kg⁻¹ h⁻¹)",
       fill = "Sediment Type") +
  geom_text(aes(label = paste0(round(Absolute_Difference, 1), "\n(", round(Percent_Difference, 1), "%)")),
            hjust = -0.1, size = 3)

print(p_abs_diff)

# Summary of largest differences
cat("\n=== LARGEST ABSOLUTE DIFFERENCES ===\n")
largest_diffs <- absolute_differences %>%
  arrange(desc(Absolute_Difference)) %>%
  head(3)

for(i in 1:nrow(largest_diffs)) {
  cat(paste0(i, ". ", largest_diffs$Sediment_Type[i], " - ", largest_diffs$Comparison[i], 
             ": ", round(largest_diffs$Absolute_Difference[i], 2), " mg O₂ kg⁻¹ h⁻¹ ",
             "(", round(largest_diffs$Percent_Difference[i], 1), "% difference)\n"))
}
# ===== Export results =======
# Save plots
ggsave("Figures/Figure2_Respiration_by_DOC_Treatment.png", fig2, width = 12, height = 8, dpi = 300)
ggsave("Figures/Respiration_by_DOC_Treatment.png", fig3, width = 12, height = 8, dpi = 300)
ggsave("Figures/NPOC_by_DOC_Treatment.png", fig_npoc, width = 12, height = 8, dpi = 300)
ggsave("Figures/CO2_Production_by_DOC_Treatment.png", fig_co2, width = 12, height = 8, dpi = 300)
ggsave("Figures/Absolute_differences.png", p_abs_diff, width = 12, height = 8, dpi = 300)

# Save effect size table
write.csv(effect_size_summary, "Data/Effect_Size_Summary_DOC_vs_Control.csv", row.names = FALSE)

print("Analysis complete! The effect size analysis specifically compares:")
print("1. PyOM treatment (H) vs Synthetic water control (S)")
print("2. Unburned DOM treatment (U) vs Synthetic water control (S)")
print("3. Combined DOC treatments (H+U) vs Control (S)")
print("4. Effects within each sediment type (dry vs wet)")