# Load required libraries
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)

# Read in the data files

# Read in the data files
do_timeseries <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv", skip = 14, header = FALSE, stringsAsFactors = FALSE)
colnames(do_timeseries) <- c("Field_Name","Sample_Name","IGSN","Material","DateTime",
                  "Elapsed_Seconds","Temperature_degreesC","DO_mg_per_L",
                  "Firesting_Serial_Number","Methods_Deviation")

npoc <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_NPOC_TN.csv", skip = 2) %>%
  filter(grepl('YEP',Sample_Name))
co2 <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_CO2.csv", skip = 2) %>%
  filter(grepl('YEP',Sample_Name))

# Function to extract sample ID without suffix (INC, SOC, GAS)
extract_sample_id <- function(sample_name) {
  sample_name <- gsub("_INC-", "_", sample_name)
  sample_name <- gsub("_SOC-", "_", sample_name)
  sample_name <- gsub("_GAS-", "_", sample_name)
  return(sample_name)
}

# ============================================================
# PROCESS DO TIME SERIES DATA
# ============================================================

# Calculate total O2 consumed for each sample
# Trim first 2 minutes and last 2 minutes (120 seconds each)
do_summary <- do_timeseries %>%
  mutate(
    Sample_ID = extract_sample_id(Sample_Name),
    Elapsed_Minutes = Elapsed_Seconds / 60,
    DO_mg_per_L = as.numeric(DO_mg_per_L)
  ) %>%
  group_by(Sample_ID) %>%
  mutate(
    Max_Time = max(Elapsed_Minutes, na.rm = TRUE),
    Min_Time = min(Elapsed_Minutes, na.rm = TRUE)
  ) %>%
  # Filter to remove first 2 min and last 2 min
  filter(Elapsed_Minutes >= (Min_Time + 2),
         Elapsed_Minutes <= (Max_Time - 2)) %>%
  summarise(
    # Get initial and final DO (after trimming)
    Initial_DO_mg_L = first(DO_mg_per_L),
    Final_DO_mg_L = last(DO_mg_per_L),
    Initial_Time_min = first(Elapsed_Minutes),
    Final_Time_min = last(Elapsed_Minutes),
    Total_Incubation_Time_min = Final_Time_min - Initial_Time_min,
    
    # Calculate total O2 consumed
    Total_O2_consumed_mg_L = Initial_DO_mg_L - Final_DO_mg_L,
    
    # Calculate average rate over this period
    Avg_O2_Rate_mg_L_h = (Total_O2_consumed_mg_L / Total_Incubation_Time_min) * 60,
    
    .groups = 'drop'
  )

# Check the results
print("DO Summary after trimming first/last 2 minutes:")
print(head(do_summary))

# Prepare NPOC data - CONVERT TO NUMERIC
npoc_clean <- npoc %>%
  filter(Material == "Sediment") %>%
  mutate(
    Sample_ID = extract_sample_id(Sample_Name),
    NPOC_mg_per_L = as.numeric(Extractable_NPOC_mg_per_L)
  ) %>%
  select(Sample_ID, NPOC_mg_per_L)

# Prepare CO2 data - CONVERT TO NUMERIC and REMOVE -9999 values
co2_clean <- co2 %>%
  mutate(
    Sample_ID = extract_sample_id(Sample_Name),
    CO2_moles_per_L = as.numeric(Partial_Pressure_CO2_moles_per_L)
  ) %>%
  # Remove -9999 values (missing data)
  filter(CO2_moles_per_L != -9999, !is.na(CO2_moles_per_L)) %>%
  select(Sample_ID, CO2_moles_per_L)

# Merge all data
merged_data <- do_summary %>%
  left_join(npoc_clean, by = "Sample_ID") %>%
  left_join(co2_clean, by = "Sample_ID")

# REMOVE samples with missing CO2 data
merged_data <- merged_data %>%
  filter(!is.na(CO2_moles_per_L))

# Add grouping variables
merged_data <- merged_data %>%
  mutate(
    Site = ifelse(grepl("^YEP1", Sample_ID), "YEP1", "YEP2"),
    Sediment_Type = ifelse(Site == "YEP1", "Dry", "Wet"),
    Treatment_Letter = sub(".*_([HUS])\\d+$", "\\1", Sample_ID),
    Treatment = case_when(
      Treatment_Letter == "S" ~ "Control (No DOM)",
      Treatment_Letter == "U" ~ "Unburned DOM",
      Treatment_Letter == "H" ~ "High Severity DOM",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment))

# ============================================================
# STOICHIOMETRIC CALCULATIONS WITH TOTAL O2 CONSUMED
# ============================================================

# Constants
water_volume_L <- 0.030  # 30 mL = 0.030 L
C_molar_mass <- 12  # g/mol
O2_molar_mass <- 32  # g/mol
CO2_molar_mass <- 44  # g/mol

# Calculate molar concentrations
stoich_calculations <- merged_data %>%
  mutate(
    # DOC added (convert from mg C/L to mol C/L)
    DOC_added_mg_C_per_L = case_when(
      Treatment == "Control (No DOM)" ~ 0,
      Treatment == "Unburned DOM" ~ 9,
      Treatment == "High Severity DOM" ~ 9,
      TRUE ~ NA_real_
    ),
    DOC_added_mol_per_L = (DOC_added_mg_C_per_L / 1000) / C_molar_mass,
    
    # NPOC measured (convert from mg C/L to mol C/L)
    NPOC_mol_per_L = (NPOC_mg_per_L / 1000) / C_molar_mass,
    
    # O2 consumed - USING TOTAL CONSUMPTION
    O2_consumed_mg_per_L = Total_O2_consumed_mg_L,
    
    # Convert O2 to moles
    O2_consumed_mol_per_L = (O2_consumed_mg_per_L / 1000) / O2_molar_mass,
    
    # CO2 is already in mol/L
    CO2_produced_mol_per_L = CO2_moles_per_L,
    
    # ============================================================
    # CALCULATE RATIOS
    # ============================================================
    
    # Respiratory Quotient: CO2 produced / O2 consumed
    RQ = ifelse(O2_consumed_mol_per_L > 0, 
                CO2_produced_mol_per_L / O2_consumed_mol_per_L, 
                NA_real_),
    
    # Carbon mineralized based on O2 consumption
    # Assuming C + O2 -> CO2, moles C mineralized = moles O2 consumed
    C_mineralized_from_O2_mol_per_L = O2_consumed_mol_per_L,
    
    # Carbon mineralized based on CO2 production
    C_mineralized_from_CO2_mol_per_L = CO2_produced_mol_per_L,
    
    # Ratio of C mineralized to DOC added
    C_mineralized_to_DOC_ratio_O2 = ifelse(DOC_added_mol_per_L > 0,
                                           C_mineralized_from_O2_mol_per_L / DOC_added_mol_per_L,
                                           NA_real_),
    
    C_mineralized_to_DOC_ratio_CO2 = ifelse(DOC_added_mol_per_L > 0,
                                            C_mineralized_from_CO2_mol_per_L / DOC_added_mol_per_L,
                                            NA_real_),
    
    # Percent of DOC mineralized
    Percent_DOC_Mineralized_O2 = C_mineralized_to_DOC_ratio_O2 * 100,
    Percent_DOC_Mineralized_CO2 = C_mineralized_to_DOC_ratio_CO2 * 100
  )

# ============================================================
# SUMMARY STATISTICS WITH MOLAR RATIOS
# ============================================================

summary_molar <- stoich_calculations %>%
  group_by(Sediment_Type, Treatment) %>%
  summarise(
    # Sample size
    n = n(),
    
    # Incubation time (trimmed)
    `Mean Incubation Time (min)` = mean(Total_Incubation_Time_min, na.rm = TRUE),
    `SD Incubation Time (min)` = sd(Total_Incubation_Time_min, na.rm = TRUE),
    
    # DOC added (mol/L)
    `DOC Added (mol C/L)` = mean(DOC_added_mol_per_L, na.rm = TRUE),
    
    # Total O2 consumed (mg/L)
    `Total O2 Consumed Mean (mg/L)` = mean(O2_consumed_mg_per_L, na.rm = TRUE),
    `Total O2 Consumed SD (mg/L)` = sd(O2_consumed_mg_per_L, na.rm = TRUE),
    
    # O2 consumed (mol/L)
    `O2 Consumed Mean (mol/L)` = mean(O2_consumed_mol_per_L, na.rm = TRUE),
    `O2 Consumed SD (mol/L)` = sd(O2_consumed_mol_per_L, na.rm = TRUE),
    
    # CO2 produced (mol/L)
    `CO2 Produced Mean (mol/L)` = mean(CO2_produced_mol_per_L, na.rm = TRUE),
    `CO2 Produced SD (mol/L)` = sd(CO2_produced_mol_per_L, na.rm = TRUE),
    
    # C mineralized from O2 (mol/L)
    `C Mineralized (from O2) Mean (mol/L)` = mean(C_mineralized_from_O2_mol_per_L, na.rm = TRUE),
    `C Mineralized (from O2) SD (mol/L)` = sd(C_mineralized_from_O2_mol_per_L, na.rm = TRUE),
    
    # C mineralized from CO2 (mol/L)
    `C Mineralized (from CO2) Mean (mol/L)` = mean(C_mineralized_from_CO2_mol_per_L, na.rm = TRUE),
    `C Mineralized (from CO2) SD (mol/L)` = sd(C_mineralized_from_CO2_mol_per_L, na.rm = TRUE),
    
    # Respiratory Quotient
    `RQ Mean (CO2/O2)` = mean(RQ, na.rm = TRUE),
    `RQ SD` = sd(RQ, na.rm = TRUE),
    
    # Ratio of C mineralized to DOC added
    `C_mineralized:DOC_added (O2 basis) Mean` = mean(C_mineralized_to_DOC_ratio_O2, na.rm = TRUE),
    `C_mineralized:DOC_added (O2 basis) SD` = sd(C_mineralized_to_DOC_ratio_O2, na.rm = TRUE),
    
    `C_mineralized:DOC_added (CO2 basis) Mean` = mean(C_mineralized_to_DOC_ratio_CO2, na.rm = TRUE),
    `C_mineralized:DOC_added (CO2 basis) SD` = sd(C_mineralized_to_DOC_ratio_CO2, na.rm = TRUE),
    
    # Percent DOC mineralized
    `Percent DOC Mineralized (O2 basis) Mean` = mean(Percent_DOC_Mineralized_O2, na.rm = TRUE),
    `Percent DOC Mineralized (O2 basis) SD` = sd(Percent_DOC_Mineralized_O2, na.rm = TRUE),
    
    `Percent DOC Mineralized (CO2 basis) Mean` = mean(Percent_DOC_Mineralized_CO2, na.rm = TRUE),
    `Percent DOC Mineralized (CO2 basis) SD` = sd(Percent_DOC_Mineralized_CO2, na.rm = TRUE),
    
    .groups = 'drop'
  ) %>%
  arrange(Sediment_Type, Treatment)

# Print summary
print(summary_molar)

# Save comprehensive summary
write.csv(summary_molar, "summary_molar_ratios_TOTAL_O2_cleaned.csv", row.names = FALSE)

# ============================================================
# VISUALIZATION: MOLAR RATIOS
# ============================================================

# Plot 1: RQ by treatment
p1 <- ggplot(stoich_calculations, aes(x = Treatment, y = RQ, fill = Sediment_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_bw() +
  labs(title = "Respiratory Quotient (RQ = CO2/O2)",
       subtitle = "Based on total O2 consumed (trimmed first/last 2 min). Red line = RQ of 1",
       y = "RQ (mol CO2 / mol O2)",
       x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("RQ_by_treatment_TOTAL_O2_cleaned.png", p1, width = 10, height = 6)

# Plot 2: Scatter plot of O2 consumed vs CO2 produced
p2 <- ggplot(stoich_calculations, 
             aes(x = O2_consumed_mol_per_L, 
                 y = CO2_produced_mol_per_L, 
                 color = Treatment,
                 shape = Sediment_Type)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "CO2 Produced vs O2 Consumed (Molar Basis)",
       subtitle = "Based on total O2 consumed. Red dashed line shows 1:1 ratio",
       x = "O2 Consumed (mol/L)",
       y = "CO2 Produced (mol/L)") +
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific) +
  theme(legend.position = "bottom")

ggsave("O2_vs_CO2_scatter_TOTAL_O2_cleaned.png", p2, width = 10, height = 6)

# Plot 3: Ratio of C mineralized to DOC added (O2 basis)
p3 <- stoich_calculations %>%
  filter(Treatment != "Control (No DOM)") %>%
  ggplot(aes(x = Treatment, y = C_mineralized_to_DOC_ratio_O2, fill = Sediment_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Ratio of C Mineralized to DOC Added (based on O2)",
       subtitle = "Red line = 1.0 (all added DOC mineralized)",
       y = "C Mineralized : DOC Added",
       x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("C_mineralized_to_DOC_ratio_O2_cleaned.png", p3, width = 10, height = 6)

# Plot 4: Ratio of C mineralized to DOC added (CO2 basis)
p4 <- stoich_calculations %>%
  filter(Treatment != "Control (No DOM)") %>%
  ggplot(aes(x = Treatment, y = C_mineralized_to_DOC_ratio_CO2, fill = Sediment_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Ratio of C Mineralized to DOC Added (based on CO2)",
       subtitle = "Red line = 1.0 (all added DOC mineralized)",
       y = "C Mineralized : DOC Added",
       x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("C_mineralized_to_DOC_ratio_CO2_cleaned.png", p4, width = 10, height = 6)

# Plot 5: Percent DOM respired (O2 basis)
p5 <- stoich_calculations %>%
  filter(Treatment != "Control (No DOM)") %>%
  ggplot(aes(x = Treatment, y = Percent_DOC_Mineralized_O2, fill = Sediment_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Percent of Added DOM Respired (Based on O2 Consumption)",
       subtitle = "Red dashed line indicates 100% of added DOM respired",
       y = "Percent of Added DOM Respired (%)",
       x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("Percent_DOM_respired_O2_cleaned.png", p5, width = 10, height = 6)

# Plot 6: Percent DOM respired (CO2 basis)
p6 <- stoich_calculations %>%
  filter(Treatment != "Control (No DOM)") %>%
  ggplot(aes(x = Treatment, y = Percent_DOC_Mineralized_CO2, fill = Sediment_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Percent of Added DOM Respired (Based on CO2 Production)",
       subtitle = "Red dashed line indicates 100% of added DOM respired",
       y = "Percent of Added DOM Respired (%)",
       x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("Percent_DOM_respired_CO2_cleaned.png", p6, width = 10, height = 6)

# Plot 7: Compare O2-based vs CO2-based C mineralization
comparison_data <- stoich_calculations %>%
  select(Sample_ID, Sediment_Type, Treatment,
         `From O2` = C_mineralized_from_O2_mol_per_L,
         `From CO2` = C_mineralized_from_CO2_mol_per_L) %>%
  pivot_longer(cols = c(`From O2`, `From CO2`),
               names_to = "Method",
               values_to = "C_mineralized_mol_per_L")

p7 <- ggplot(comparison_data, aes(x = Treatment, y = C_mineralized_mol_per_L, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~Sediment_Type) +
  theme_bw() +
  labs(title = "Carbon Mineralized: O2-based vs CO2-based Estimates",
       y = "C Mineralized (mol C/L)",
       x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::scientific)

ggsave("C_mineralized_comparison_TOTAL_O2_cleaned.png", p7, width = 10, height = 6)

# ============================================================
# INDIVIDUAL SAMPLE CALCULATIONS (for supplementary materials)
# ============================================================

individual_samples_molar <- stoich_calculations %>%
  select(Sample_ID, Sediment_Type, Treatment,
         Total_Incubation_Time_min,
         Initial_DO_mg_L,
         Final_DO_mg_L,
         DOC_added_mol_per_L,
         O2_consumed_mol_per_L,
         CO2_produced_mol_per_L,
         C_mineralized_from_O2_mol_per_L,
         C_mineralized_from_CO2_mol_per_L,
         RQ,
         C_mineralized_to_DOC_ratio_O2,
         C_mineralized_to_DOC_ratio_CO2,
         Percent_DOC_Mineralized_O2,
         Percent_DOC_Mineralized_CO2) %>%
  arrange(Sediment_Type, Treatment, Sample_ID)

write.csv(individual_samples_molar, "individual_sample_molar_ratios_TOTAL_O2_cleaned.csv", row.names = FALSE)

# ============================================================
# INTERPRETATION GUIDE
# ============================================================

cat("\n=== DATA CLEANING NOTE ===\n")
cat("Sample YEP2A_S2 was removed due to missing CO2 data (-9999 flag)\n\n")

cat("=== INTERPRETATION GUIDE ===\n\n")
cat("All values are in MOLAR UNITS (mol/L):\n\n")
cat("1. DOC Added: Moles of carbon added as dissolved organic matter\n")
cat("   - Control: 0 mol C/L\n")
cat("   - Unburned & High Severity: 9 mg C/L = 7.5 × 10^-4 mol C/L\n\n")
cat("2. O2 Consumed: Moles of O2 consumed during ENTIRE incubation (trimmed first/last 2 min)\n\n")
cat("3. CO2 Produced: Moles of CO2 measured at end of incubation\n\n")
cat("4. C Mineralized: Moles of carbon mineralized\n")
cat("   - From O2: Assumes C + O2 -> CO2 (1:1 stoichiometry)\n")
cat("   - From CO2: Direct measurement\n\n")
cat("5. RQ (Respiratory Quotient): mol CO2 / mol O2\n")
cat("   - RQ = 1: Complete aerobic carbohydrate oxidation\n")
cat("   - RQ < 1: Protein or lipid oxidation, or CO2 dissolution\n")
cat("   - RQ > 1: Fermentation, denitrification, or native carbon respiration\n\n")
cat("6. C_mineralized:DOC_added Ratio:\n")
cat("   - Ratio = 1: All added DOC was mineralized\n")
cat("   - Ratio < 1: Only partial mineralization of added DOC\n")
cat("   - Ratio > 1: Native sediment carbon also mineralized\n\n")
cat("NOTE: The non-linear O2 consumption in YEP2 (wet) samples suggests a shift\n")
cat("      from aerobic to anaerobic metabolism. High RQ values (>1) support this.\n\n")

print("Analysis complete! Files saved:")
print("1. summary_molar_ratios_TOTAL_O2_cleaned.csv")
print("2. individual_sample_molar_ratios_TOTAL_O2_cleaned.csv")
print("3. RQ_by_treatment_TOTAL_O2_cleaned.png")
print("4. O2_vs_CO2_scatter_TOTAL_O2_cleaned.png")
print("5. C_mineralized_to_DOC_ratio_O2_cleaned.png")
print("6. C_mineralized_to_DOC_ratio_CO2_cleaned.png")
print("7. Percent_DOM_respired_O2_cleaned.png")
print("8. Percent_DOM_respired_CO2_cleaned.png")
print("9. C_mineralized_comparison_TOTAL_O2_cleaned.png")