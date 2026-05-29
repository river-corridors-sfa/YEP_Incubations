# ================================
#
# Calculate respiration for YEP experiment
#
# Author: Vanessa Garayburu-Caruso (vanessa.garayburu-caruso@pnnl.gov)
#
# ==== Loading libraries =========
rm(list = ls())
library(dplyr)
library(tidyr)
# ==== set wd ==========

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("./..")

# ==== Read in the data files ====
respiration <- read.csv("./YEP_Sediment_Incubations_Respiration_Rates.csv", skip = 2) %>%
  filter(grepl("YEP", Sample_Name))

mass <- read.csv("./YEP_Sediment_Water_Mass_Volume.csv", skip = 2) %>%
  filter(grepl('YEP',Sample_Name)) %>%
  dplyr::mutate(dplyr::across(
    c("Field_Moist_Sediment_Mass_g","Dry_Sediment_Mass_g","Water_Mass_g","Wet_Sediment_Volume_mL"),
    as.numeric
  )) %>%
  dplyr::select(-Methods_Deviation)



data = merge(respiration, mass, by.x = "Sample_Name", by.y = "Sample_Name", all.x = TRUE) 

data_final = data %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(is.na(Dry_Sediment_Mass_g), -9999,
                                                      ((Respiration_Rate_mg_DO_per_L_per_H*Water_Mass_g*0.001)/(Dry_Sediment_Mass_g*0.001)))) %>%
  mutate(across(where(is.numeric) & !Respiration_p_value, ~ round(.x, 2))) %>%
  filter(!is.na(Respiration_Rate_mg_DO_per_L_per_H)) %>%
  dplyr::select(Field_Name, Sample_Name, IGSN,Material, "Respiration_Rate_mg_DO_per_L_per_H",          
                "Respiration_R_Squared",                       
                "Respiration_R_Squared_Adj",                   
                "Respiration_p_value" ,
                "Respiration_Rate_mg_DO_per_kg_per_H",
                "Total_Incubation_Time_Min" ,                  
                "Number_Points_In_Respiration_Regression",     
                "Number_Points_Removed_Respiration_Regression",
                "DO_Concentration_At_Incubation_Time_Zero",    
                "Methods_Deviation")
# ==== Save results ====
write.csv(data_final, "./YEP_Respiration.csv", row.names = FALSE)
