# Clear workspace and load libraries
rm(list = ls())
library(dplyr)
library(ggplot2)

# ===== PARAMETERS =====
# Your experimental conditions
temp_c <- 21  # temperature of the incubation chamber
temp_k <- temp_c + 273.15
salinity_psu <- 0.1  # estimate of the salinity of columbia river synthetic water

# Ambient CO2 concentration (atmospheric equilibrium experimental conditions)
ambient_co2_mol_per_L <- 1.6e-05  # ~400 ppm atmospheric CO2 at 21°C, 0.1 PSU
plot_toggle <- TRUE

# ===== READ CO2 DATA =====
co2_data <- read.csv('YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_CO2.csv', skip = 2) %>%
  filter(grepl("YEP", Sample_Name)) %>%
  filter(Partial_Pressure_CO2_moles_per_L != -9999) %>%
  dplyr::select(Sample_Name, Partial_Pressure_CO2_ppm, Partial_Pressure_CO2_moles_per_L) %>%
  mutate(Partial_Pressure_CO2_moles_per_L = as.numeric(Partial_Pressure_CO2_moles_per_L))

# ===== READ INCUBATION TIME DATA =====
# Read the previous results to get incubation times
incubation_data <-  read.csv('YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv', skip = 2) %>%
  filter(grepl("YEP", Sample_Name)) %>%
  filter(DO_mg_per_L != -9999) %>%
  dplyr::select(Sample_Name, Elapsed_Seconds, DO_mg_per_L, Temperature_degreesC) %>%
  mutate(across(c(Elapsed_Seconds, DO_mg_per_L), as.numeric))

# Create plots directory
if(plot_toggle) {
  if(!dir.exists("plots_co2")) {
    dir.create("plots_co2")
    cat("Created 'plots_co2' directory\n")
  }
}

# ===== RESULTS DATAFRAME =====
results <- data.frame(
  Sample_Name = character(),
  CO2_Sample_Name = character(),
  Initial_CO2_mol_per_L = numeric(),
  Final_CO2_mol_per_L = numeric(),
  CO2_Change_mol_per_L = numeric(),
  Net_CO2_Production_mol_per_L = numeric(),
  Total_Incubation_Time_Min = numeric(),
  Total_Incubation_Time_H = numeric(),
  CO2_Production_Rate_mol_per_L_per_Min = numeric(),
  CO2_Production_Rate_mol_per_L_per_H = numeric(),
  Final_CO2_ppm = numeric(),
  Temperature_C = numeric(),
  Salinity_PSU = numeric(),
  stringsAsFactors = FALSE
)

# ===== MAIN ANALYSIS LOOP =====
for(i in 1:nrow(co2_data)) {
  co2_sample_name <- co2_data$Sample_Name[i]
  
  cat("\n=== Processing CO2 sample", i, "of", nrow(co2_data), ":", co2_sample_name, "===\n")
  
  # Convert CO2 sample name to match incubation data (GAS -> INC)
  inc_sample_name <- gsub("_GAS-", "_INC-", co2_sample_name)
  
  # Get the corresponding incubation data from test file
  sample_inc_data <- incubation_data %>%
    filter(Sample_Name == inc_sample_name) %>%
    arrange(Elapsed_Seconds)
  
  if(nrow(sample_inc_data) == 0) {
    cat("  No matching incubation data found for:", inc_sample_name, "\n")
    next
  }
  
  # Calculate incubation time directly from test file - USE ALL TIME POINTS
  start_time <- min(sample_inc_data$Elapsed_Seconds)
  end_time <- max(sample_inc_data$Elapsed_Seconds)
  total_time_sec <- end_time - start_time
  total_time_min <- total_time_sec / 60
  total_time_h <- total_time_min / 60
  
  # Get CO2 data
  final_co2_mol_per_L <- as.numeric(co2_data$Partial_Pressure_CO2_moles_per_L[i])
  final_co2_ppm <- as.numeric(co2_data$Partial_Pressure_CO2_ppm[i])
  
  # Calculate CO2 production using formula: (Final CO2 - Ambient CO2) / Incubation Time
  initial_co2_mol_per_L <- ambient_co2_mol_per_L  # Assumed ambient
  co2_change <- final_co2_mol_per_L - initial_co2_mol_per_L
  net_co2_production <- co2_change  # Same as change since we're measuring vs ambient
  
  # Calculate production rates
  rate_per_min <- net_co2_production / total_time_min
  rate_per_h <- net_co2_production / total_time_h
  
  cat("  Original sample:", inc_sample_name, "\n")
  cat("  Incubation time:", round(total_time_min, 2), "minutes (", round(total_time_h, 3), "hours)\n")
  cat("  Initial CO2 (assumed ambient):", format(initial_co2_mol_per_L, scientific = TRUE), "mol/L\n")
 # cat("  Final CO2:", format(final_co2_mol_per_L, scientific = TRUE), "mol/L (", round(final_co2_ppm, 1), "ppm)\n")
  cat("  Net CO2 production:", format(net_co2_production, scientific = TRUE), "mol/L\n")
  cat("  Production rate:", format(rate_per_h, scientific = TRUE), "mol/L/h\n")
  
  # Store results
  result_row <- data.frame(
    Sample_Name = inc_sample_name,
    CO2_Sample_Name = co2_sample_name,
    Initial_CO2_mol_per_L = initial_co2_mol_per_L,
    Final_CO2_mol_per_L = final_co2_mol_per_L,
    CO2_Change_mol_per_L = co2_change,
    Net_CO2_Production_mol_per_L = net_co2_production,
    Total_Incubation_Time_Min = total_time_min,
    Total_Incubation_Time_H = total_time_h,
    CO2_Production_Rate_mol_per_L_per_Min = rate_per_min,
    CO2_Production_Rate_mol_per_L_per_H = rate_per_h,
    Final_CO2_ppm = final_co2_ppm,
    Temperature_C = temp_c,
    Salinity_PSU = salinity_psu,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(results, result_row)
}

# ===== CREATE VISUALIZATIONS =====
if(plot_toggle && nrow(results) > 0) {
  # Plot CO2 production rates
  p1 <- ggplot(results, aes(x = reorder(Sample_Name, CO2_Production_Rate_mol_per_L_per_H), 
                            y = CO2_Production_Rate_mol_per_L_per_H * 1e6)) +  # Convert to µmol/L/h for better readability
    geom_col(fill = "darkgreen", alpha = 0.7) +
    coord_flip() +
    theme_bw() +
    labs(title = "CO2 Production Rates by Sample",
         x = "Sample",
         y = "CO2 Production Rate (µmol/L/h)",
         caption = paste("Conditions:", temp_c, "°C,", salinity_psu, "PSU",
                         "\nAmbient CO2:", format(ambient_co2_mol_per_L, scientific = TRUE), "mol/L")) +
    theme(axis.text.y = element_text(size = 8))
  
  print(p1)
  ggsave("plots_co2/CO2_production_rates_by_sample.pdf", p1, width = 12, height = 8)
  
  # Plot final CO2 concentrations
  p2 <- ggplot(results, aes(x = reorder(Sample_Name, Final_CO2_mol_per_L), 
                            y = Final_CO2_ppm)) +
    geom_col(fill = "brown", alpha = 0.7) +
    geom_hline(yintercept = 400, linetype = "dashed", color = "red") +  # atmospheric CO2 in ppm
    coord_flip() +
    theme_bw() +
    labs(title = "Final CO2 Concentrations by Sample",
         x = "Sample",
         y = "Final CO2 Concentration (ppm)",
         caption = paste("Red line = Atmospheric CO2 (~400 ppm)",
                         "\nConditions:", temp_c, "°C,", salinity_psu, "PSU")) +
    theme(axis.text.y = element_text(size = 8))
  
  print(p2)
  ggsave("plots_co2/Final_CO2_concentrations_by_sample.pdf", p2, width = 12, height = 8)
  
  # Create a comparison plot: CO2 production vs DO consumption
  # Read DO respiration data for comparison
  if(file.exists("sediment_respiration_comparison_results.csv")) {
    do_rates <- read.csv("sediment_respiration_comparison_results.csv") %>%
      select(Sample_Name, original_rate_mg_L_per_h)
    
    # Merge with CO2 data
    comparison_data <- results %>%
      left_join(do_rates, by = "Sample_Name") %>%
      filter(!is.na(original_rate_mg_L_per_h))
    
    if(nrow(comparison_data) > 0) {
      p3 <- ggplot(comparison_data, aes(x = original_rate_mg_L_per_h, 
                                        y = CO2_Production_Rate_mol_per_L_per_H * 1e6)) +
        geom_point(size = 3, alpha = 0.7, color = "blue") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        theme_bw() +
        labs(title = "CO2 Production vs DO Consumption",
             x = "DO Consumption Rate (mg/L/h)",
             y = "CO2 Production Rate (µmol/L/h)",
             caption = paste("Each point represents one sample",
                             "\nConditions:", temp_c, "°C,", salinity_psu, "PSU")) +
        geom_text(aes(label = gsub("YEP._INC-", "", Sample_Name)), 
                  vjust = -0.5, hjust = 0.5, size = 2)
      
      print(p3)
      ggsave("plots_co2/CO2_vs_DO_comparison.pdf", p3, width = 10, height = 8)
      
      # Calculate correlation
      correlation <- cor(comparison_data$original_rate_mg_L_per_h, 
                         comparison_data$CO2_Production_Rate_mol_per_L_per_H, 
                         use = "complete.obs")
      cat("Correlation between DO consumption and CO2 production:", round(correlation, 3), "\n")
    }
  }
  
  # Create incubation time comparison plot
  p4 <- ggplot(results, aes(x = reorder(Sample_Name, Total_Incubation_Time_Min), 
                            y = Total_Incubation_Time_Min)) +
    geom_col(fill = "purple", alpha = 0.7) +
    coord_flip() +
    theme_bw() +
    labs(title = "Incubation Times by Sample",
         x = "Sample",
         y = "Total Incubation Time (minutes)") +
    theme(axis.text.y = element_text(size = 8))
  
  print(p4)
  ggsave("plots_co2/Incubation_times_by_sample.pdf", p4, width = 12, height = 8)
}

# ===== RESULTS SUMMARY =====
cat("\n=== CO2 ANALYSIS COMPLETE ===\n")
cat("Successfully processed:", nrow(results), "samples\n")
cat("Experimental conditions:", temp_c, "°C,", salinity_psu, "PSU\n")
cat("Assumed ambient CO2:", format(ambient_co2_mol_per_L, scientific = TRUE), "mol/L\n")

# Summary statistics
cat("\nCO2 Production Rate Summary:\n")
cat("Mean rate (mol/L/h):", format(mean(results$CO2_Production_Rate_mol_per_L_per_H, na.rm = TRUE), scientific = TRUE), "\n")
cat("Mean rate (µmol/L/h):", round(mean(results$CO2_Production_Rate_mol_per_L_per_H, na.rm = TRUE) * 1e6, 2), "\n")
cat("Rate range (µmol/L/h):", round(range(results$CO2_Production_Rate_mol_per_L_per_H, na.rm = TRUE) * 1e6, 2), "\n")

cat("\nIncubation Time Summary:\n")
cat("Mean time (min):", round(mean(results$Total_Incubation_Time_Min, na.rm = TRUE), 2), "\n")
cat("Time range (min):", round(range(results$Total_Incubation_Time_Min, na.rm = TRUE), 2), "\n")

cat("\nFinal CO2 Concentration Summary:\n")
cat("Mean (ppm):", round(mean(results$Final_CO2_ppm, na.rm = TRUE), 1), "\n")
cat("Range (ppm):", round(range(results$Final_CO2_ppm, na.rm = TRUE), 1), "\n")

# Show results
print(head(results))

# Save results
write.csv(results, "co2_production_rates.csv", row.names = FALSE)
cat("Results saved to: co2_production_rates.csv\n")