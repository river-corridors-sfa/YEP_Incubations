# Load required libraries
library(dplyr)
library(tidyr)

# Read in the data files
respiration <- read.csv("data/segmented_respiration_analysis.csv")
npoc <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_NPOC_TN.csv", skip = 2) %>%
  filter(grepl('YEP',Sample_Name))
co2 <- read.csv("YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_CO2.csv", skip = 2) %>%
  filter(grepl('YEP',Sample_Name))

# Function to extract sample ID without suffix (INC, SOC, GAS)
extract_sample_id <- function(sample_name) {
  # Remove the suffixes _INC-, _SOC-, _GAS-
  sample_name <- gsub("_INC-", "_", sample_name)
  sample_name <- gsub("_SOC-", "_", sample_name)
  sample_name <- gsub("_GAS-", "_", sample_name)
  return(sample_name)
}

# Prepare respiration data
resp_clean <- respiration %>%
  mutate(Sample_ID = extract_sample_id(Sample_Name)) %>%
  select(Sample_ID, Rate_break_1 = Rate_mg_L_h_break_1)

# Prepare NPOC data - CONVERT TO NUMERIC
npoc_clean <- npoc %>%
  filter(Material == "Sediment") %>%
  mutate(
    Sample_ID = extract_sample_id(Sample_Name),
    # Convert to numeric, which will turn non-numeric values to NA
    NPOC_mg_per_kg = as.numeric(Extractable_NPOC_mg_per_kg)
  ) %>%
  select(Sample_ID, NPOC_mg_per_kg)

# Prepare CO2 data - CONVERT TO NUMERIC
co2_clean <- co2 %>%
  mutate(
    Sample_ID = extract_sample_id(Sample_Name),
    # Convert to numeric
    CO2_moles_per_L = as.numeric(Partial_Pressure_CO2_moles_per_L)
  ) %>%
  select(Sample_ID, CO2_moles_per_L)

# Merge all data
merged_data <- resp_clean %>%
  left_join(npoc_clean, by = "Sample_ID") %>%
  left_join(co2_clean, by = "Sample_ID")

# Add grouping variables
merged_data <- merged_data %>%
  mutate(
    # Extract site (YEP1 or YEP2)
    Site = ifelse(grepl("^YEP1", Sample_ID), "YEP1", "YEP2"),
    
    # Define sediment type
    Sediment_Type = ifelse(Site == "YEP1", "Dry", "Wet"),
    
    # Extract treatment letter
    Treatment_Letter = sub(".*_([HUS])\\d+$", "\\1", Sample_ID),
    
    # Define treatment name
    Treatment = case_when(
      Treatment_Letter == "S" ~ "Control (No DOM)",
      Treatment_Letter == "U" ~ "Unburned DOM",
      Treatment_Letter == "H" ~ "High Severity DOM",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment))  # Remove any samples that don't match treatment pattern

# Create summary table with proper column names
summary_table <- merged_data %>%
  group_by(Sediment_Type, Treatment) %>%
  summarise(
    # Max O2 consumption rate statistics (mg/L/h)
    `Max O2 Consumption Rate Mean (mg/L/h)` = mean(Rate_break_1, na.rm = TRUE),
    `Max O2 Consumption Rate Median (mg/L/h)` = median(Rate_break_1, na.rm = TRUE),
    `Max O2 Consumption Rate SD (mg/L/h)` = sd(Rate_break_1, na.rm = TRUE),
    
    # NPOC statistics (mg/kg dry sediment)
    `NPOC Mean (mg/kg dry sediment)` = mean(NPOC_mg_per_kg, na.rm = TRUE),
    `NPOC Median (mg/kg dry sediment)` = median(NPOC_mg_per_kg, na.rm = TRUE),
    `NPOC SD (mg/kg dry sediment)` = sd(NPOC_mg_per_kg, na.rm = TRUE),
    
    # CO2 statistics (mol/L)
    `CO2 Mean (mol/L)` = mean(CO2_moles_per_L, na.rm = TRUE),
    `CO2 Median (mol/L)` = median(CO2_moles_per_L, na.rm = TRUE),
    `CO2 SD (mol/L)` = sd(CO2_moles_per_L, na.rm = TRUE),
    
    n = n(),
    .groups = 'drop'  # This removes the grouping warning
  ) %>%
  arrange(Sediment_Type, Treatment)

# Print the table
print(summary_table)

# Save as CSV for easy pasting into manuscript
#write.csv(summary_table, "summary_statistics_table.csv", row.names = FALSE)

# Create a formatted version for manuscript (rounded values)
manuscript_table <- summary_table %>%
  mutate(
    across(contains("Max O2"), ~round(.x, 3)),
    across(contains("NPOC"), ~round(.x, 2)),
    across(contains("CO2"), ~format(.x, scientific = TRUE, digits = 3))
  )

# Print formatted table
print(manuscript_table)

# Save formatted version
write.csv(manuscript_table, "Data/summary_statistics_formatted.csv", row.names = FALSE)
