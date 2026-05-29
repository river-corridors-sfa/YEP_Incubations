rm(list = ls())

library(tidyverse)

# Run from repo root or from src/.
if (!dir.exists("v2_data") && dir.exists("../v2_data")) {
  setwd("..")
}

out_dir <- "Tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

respiration_path <- file.path(
  "v2_data", "v2_YEP_Sample_Data",
  "YEP_Sediment_Incubations_Respiration_Rates.csv"
)
npoc_path <- file.path("v2_data", "v2_YEP_Sample_Data", "YEP_Sediment_NPOC_TN.csv")
gas_path <- file.path("v2_data", "v2_YEP_Sample_Data", "v2_YEP_Sediment_CO2_CH4_N2O.csv")
ions_path <- file.path("v2_data", "v2_YEP_Sample_Data", "YEP_Sediment_Ions.csv")

read_yep_table <- function(path) {
  read_csv(
    path,
    skip = 2,
    na = c("", "NA", "N/A"),
    show_col_types = FALSE,
    name_repair = "minimal",
    trim_ws = TRUE
  ) %>%
    filter(str_detect(Sample_Name, "^YEP"))
}

to_number <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^<", "")
  x[x %in% c("", "NA", "N/A", "-9999")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

make_sample_id <- function(sample_name) {
  sample_name %>%
    str_replace("_INC-", "_") %>%
    str_replace("_SOC-", "_") %>%
    str_replace("_GAS-", "_") %>%
    str_replace("_SIN-", "_")
}

respiration <- read_yep_table(respiration_path) %>%
  filter(str_detect(Sample_Name, "^YEP[12].*_INC-[HSU]")) %>%
  transmute(
    Sample_ID = make_sample_id(Sample_Name),
    Total_O2_Consumption_Rate_mg_h_kg = to_number(
      Normalized_Respiration_Rate_mg_DO_per_H_per_kg_dry_sediment
    )
  )

npoc <- read_yep_table(npoc_path) %>%
  filter(Material == "Sediment", str_detect(Sample_Name, "^YEP[12].*_SOC-[HSU]")) %>%
  transmute(
    Sample_ID = make_sample_id(Sample_Name),
    DOC_mg_kg = to_number(Extractable_NPOC_mg_per_kg)
  )

co2 <- read_yep_table(gas_path) %>%
  filter(Material == "Sediment", str_detect(Sample_Name, "^YEP[12].*_GAS-[HSU]")) %>%
  transmute(
    Sample_ID = make_sample_id(Sample_Name),
    CO2_Production_Rate_mol_L_h = to_number(Rate_CO2_moles_per_L_per_hr)
  )

nitrate <- read_yep_table(ions_path) %>%
  filter(str_detect(Sample_Name, "^YEP[12].*_SIN-[HSU]")) %>%
  transmute(
    Sample_ID = make_sample_id(Sample_Name),
    NO3_N_mg_L = to_number(`00618_NO3_mg_per_L_as_N`)
  )

table_data <- respiration %>%
  full_join(co2, by = "Sample_ID") %>%
  full_join(npoc, by = "Sample_ID") %>%
  full_join(nitrate, by = "Sample_ID") %>%
  mutate(
    Sediment_Type = case_when(
      str_detect(Sample_ID, "^YEP1") ~ "Dry",
      str_detect(Sample_ID, "^YEP2") ~ "Wet",
      TRUE ~ NA_character_
    ),
    Treatment_Letter = str_extract(Sample_ID, "[HSU](?=\\d+$)"),
    Treatment = case_when(
      Treatment_Letter == "S" ~ "Control",
      Treatment_Letter == "U" ~ "Unburned DOM",
      Treatment_Letter == "H" ~ "High Severity DOM",
      TRUE ~ NA_character_
    ),
    Sediment_Type = factor(Sediment_Type, levels = c("Dry", "Wet")),
    Treatment = factor(Treatment, levels = c("Control", "Unburned DOM", "High Severity DOM"))
  ) %>%
  filter(!is.na(Sediment_Type), !is.na(Treatment))

summary_table <- table_data %>%
  group_by(Sediment_Type, Treatment) %>%
  summarise(
    n_respiration = sum(!is.na(Total_O2_Consumption_Rate_mg_h_kg)),
    `Mean Total O2 Consumption Rate (mg h-1 kg-1)` =
      mean(Total_O2_Consumption_Rate_mg_h_kg, na.rm = TRUE),
    `Median Total O2 Consumption Rate (mg h-1 kg-1)` =
      median(Total_O2_Consumption_Rate_mg_h_kg, na.rm = TRUE),
    `SD Total O2 Consumption Rate (mg h-1 kg-1)` =
      sd(Total_O2_Consumption_Rate_mg_h_kg, na.rm = TRUE),
    n_co2 = sum(!is.na(CO2_Production_Rate_mol_L_h)),
    `Mean CO2 production rate (mol L-1 h-1)` =
      mean(CO2_Production_Rate_mol_L_h, na.rm = TRUE),
    `Median CO2 production rate (mol L-1 h-1)` =
      median(CO2_Production_Rate_mol_L_h, na.rm = TRUE),
    `SD CO2 production rate (mol L-1 h-1)` =
      sd(CO2_Production_Rate_mol_L_h, na.rm = TRUE),
    n_doc = sum(!is.na(DOC_mg_kg)),
    `Mean DOC (mg kg-1)` = mean(DOC_mg_kg, na.rm = TRUE),
    `Median DOC (mg kg-1)` = median(DOC_mg_kg, na.rm = TRUE),
    `SD DOC (mg kg-1)` = sd(DOC_mg_kg, na.rm = TRUE),
    n_no3 = sum(!is.na(NO3_N_mg_L)),
    `Mean NO3-N (mg L-1)` = mean(NO3_N_mg_L, na.rm = TRUE),
    `Median NO3-N (mg L-1)` = median(NO3_N_mg_L, na.rm = TRUE),
    `SD NO3-N (mg L-1)` = sd(NO3_N_mg_L, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Sediment_Type, Treatment)

formatted_table <- summary_table %>%
  mutate(
    across(matches("^(Mean|Median|SD).*O2"), ~ round(.x, 2)),
    across(matches("^(Mean|Median|SD).*CO2"), ~ format(.x, scientific = TRUE, digits = 3)),
    across(matches("^(Mean|Median|SD).*DOC"), ~ round(.x, 2)),
    across(matches("^(Mean|Median|SD).*NO3"), ~ round(.x, 2))
  )

write_csv(summary_table, file.path(out_dir, "Table1_summary_statistics.csv"))
write_csv(formatted_table, file.path(out_dir, "Table1_summary_statistics_formatted.csv"))

print(formatted_table)
