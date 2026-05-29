# ================================
#
# Calculate respiration for YEP experiment
#
# Author: Vanessa Garayburu-Caruso (vanessa.garayburu-caruso@pnnl.gov)
#
# ==== Loading libraries =========
# Clear workspace and load libraries
rm(list = ls())
library(dplyr)
library(segmented)
library(ggplot2)

# ==== set wd ==========

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("./..")

# ===== PARAMETERS =====
plot_toggle <- TRUE

# ===== READ DATA =====
data <- read.csv('./YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv', skip = 2) %>%
  filter(grepl("YEP", Sample_Name)) %>%
  filter(DO_mg_per_L != -9999) %>%
  dplyr::select(Sample_Name, Elapsed_Seconds, DO_mg_per_L, Temperature_degreesC) %>%
  mutate(across(c(Elapsed_Seconds, DO_mg_per_L), as.numeric))

unique_samples <- unique(data$Sample_Name)
print(paste("Found", length(unique_samples), "samples"))

# ===== RESULTS DATAFRAME =====
results <- data.frame(
  Sample_Name = character(),
  
  # Segment 1 (before breakpoint)
  Slope_break_1 = numeric(),
  Slope_break_1_SE = numeric(),
  Intercept_break_1 = numeric(),
  R_squared_break_1 = numeric(),
  Adjusted_R_squared_break_1 = numeric(),
  P_value_break_1 = numeric(),
  N_observations_break_1 = numeric(),
  Residual_SE_break_1 = numeric(),
  Total_Incubation_Time_Min_break_1 = numeric(),
  Number_Points_In_Respiration_Regression_break_1 = numeric(),
  Number_Points_Removed_Respiration_Regression_break_1 = numeric(),
  DO_Concentration_At_Incubation_Time_Zero_break_1 = numeric(),
  Rate_mg_L_h_break_1 = numeric(),  # NEW: Rate in mg/L/h with original sign
  
  # Segment 2 (after breakpoint)
  Slope_break_2 = numeric(),
  Slope_break_2_SE = numeric(),
  Intercept_break_2 = numeric(),
  R_squared_break_2 = numeric(),
  Adjusted_R_squared_break_2 = numeric(),
  P_value_break_2 = numeric(),
  N_observations_break_2 = numeric(),
  Residual_SE_break_2 = numeric(),
  Total_Incubation_Time_Min_break_2 = numeric(),
  Number_Points_In_Respiration_Regression_break_2 = numeric(),
  Number_Points_Removed_Respiration_Regression_break_2 = numeric(),
  DO_Concentration_At_Incubation_Time_Zero_break_2 = numeric(),
  Rate_mg_L_h_break_2 = numeric(),  # NEW: Rate in mg/L/h with original sign
  
  # Overall info
  Breakpoint_Minutes = numeric(),
  Segment1_Steeper_Than_Segment2 = logical(),
  Segmented_Successful = logical(),
  
  stringsAsFactors = FALSE
)

# ===== MAIN ANALYSIS LOOP =====
for(i in 1:length(unique_samples)) {
  cat("\n=== Processing sample", i, "of", length(unique_samples), ":", unique_samples[i], "===\n")
  
  # Get sample data
  sample_data <- data %>%
    filter(Sample_Name == unique_samples[i]) %>%
    arrange(Elapsed_Seconds)
  
  if(nrow(sample_data) < 6) {
    cat("  Skipping - insufficient data points\n")
    next
  }
  
  # Calculate relative time from start
  start_time <- min(sample_data$Elapsed_Seconds)
  sample_data$relative_time_min <- (sample_data$Elapsed_Seconds - start_time) / 60
  sample_data$relative_time_sec <- sample_data$Elapsed_Seconds - start_time
  
  original_n <- nrow(sample_data)
  max_time <- max(sample_data$relative_time_min)
  
  cat("  Original data:", original_n, "points, time range: 0 to", round(max_time, 2), "min\n")
  
  # STEP 1: Remove first 2 minutes
  sample_data <- sample_data %>%
    filter(relative_time_min > 2)
  cat("  After removing first 2 min:", nrow(sample_data), "points\n")
  
  # STEP 2: Remove last 2 minutes  
  if(nrow(sample_data) > 0) {
    current_max_time <- max(sample_data$relative_time_min)
    sample_data <- sample_data %>%
      filter(relative_time_min < (current_max_time - 2))
    cat("  After removing last 2 min:", nrow(sample_data), "points\n")
  }
  
  # STEP 3: Remove points below 2 mg/L
  sample_data <- sample_data %>%
    filter(DO_mg_per_L >= 2)
  cat("  After removing DO < 2 mg/L:", nrow(sample_data), "points\n")
  
  # Check if enough points remain
  if(nrow(sample_data) < 4) {
    cat("  Skipping - insufficient points after filtering\n")
    next
  }
  
  # STEP 4: RESET TIME TO ZERO after filtering
  min_time_after_filtering <- min(sample_data$relative_time_min)
  sample_data$relative_time_min_reset <- sample_data$relative_time_min - min_time_after_filtering
  sample_data$relative_time_sec_reset <- sample_data$relative_time_sec - (min_time_after_filtering * 60)
  
  cat("  Time reset: new range 0 to", round(max(sample_data$relative_time_min_reset), 2), "min\n")
  
  # STEP 5: Perform segmented regression
  tryCatch({
    # Initial linear fit using reset time
    initial_fit <- lm(DO_mg_per_L ~ relative_time_sec_reset, data = sample_data)
    
    # Segmented regression
    mid_point <- (max(sample_data$relative_time_sec_reset) + min(sample_data$relative_time_sec_reset)) / 2
    seg_fit <- segmented(initial_fit, seg.Z = ~relative_time_sec_reset, psi = mid_point)
    
    # Get breakpoint
    breakpoint_sec <- seg_fit$psi[2]
    breakpoint_min <- breakpoint_sec / 60
    
    cat("  Segmented regression successful, breakpoint at:", round(breakpoint_min, 2), "min\n")
    
    # Split data into two segments
    segment1_data <- sample_data %>%
      filter(relative_time_sec_reset <= breakpoint_sec)
    
    segment2_data <- sample_data %>%
      filter(relative_time_sec_reset > breakpoint_sec)
    
    cat("  Segment 1:", nrow(segment1_data), "points, Segment 2:", nrow(segment2_data), "points\n")
    
    # Fit linear models to each segment
    if(nrow(segment1_data) >= 3 && nrow(segment2_data) >= 3) {
      
      # Segment 1 analysis
      fit1 <- lm(DO_mg_per_L ~ relative_time_sec_reset, data = segment1_data)
      summary1 <- summary(fit1)
      slope1 <- coef(fit1)[2]
      slope1_se <- summary1$coefficients[2, 2]
      intercept1 <- coef(fit1)[1]
      rate1_mg_L_h <- slope1 * 3600  # Convert to mg/L/h, keep original sign
      
      # Segment 2 analysis  
      fit2 <- lm(DO_mg_per_L ~ relative_time_sec_reset, data = segment2_data)
      summary2 <- summary(fit2)
      slope2 <- coef(fit2)[2]
      slope2_se <- summary2$coefficients[2, 2]
      intercept2 <- coef(fit2)[1]
      rate2_mg_L_h <- slope2 * 3600  # Convert to mg/L/h, keep original sign
      
      # Check if segment 1 is steeper (more negative) than segment 2
      segment1_steeper <- slope1 <= slope2  # More negative or equal
      
      cat("  Segment 1 rate:", round(rate1_mg_L_h, 3), "mg/L/h\n")
      cat("  Segment 2 rate:", round(rate2_mg_L_h, 3), "mg/L/h\n")
      cat("  Segment 1 steeper or equal?", segment1_steeper, "\n")
      
      # Calculate additional metrics
      total_time_1 <- max(segment1_data$relative_time_min_reset)
      total_time_2 <- max(segment2_data$relative_time_min_reset) - min(segment2_data$relative_time_min_reset)
      
      points_removed_1 <- original_n - nrow(segment1_data)
      points_removed_2 <- original_n - nrow(segment2_data)
      
      do_time_zero_1 <- segment1_data$DO_mg_per_L[1]
      do_time_zero_2 <- segment2_data$DO_mg_per_L[1]
      
      # Store results
      result_row <- data.frame(
        Sample_Name = unique_samples[i],
        
        # Segment 1
        Slope_break_1 = slope1,
        Slope_break_1_SE = slope1_se,
        Intercept_break_1 = intercept1,
        R_squared_break_1 = summary1$r.squared,
        Adjusted_R_squared_break_1 = summary1$adj.r.squared,
        P_value_break_1 = summary1$coefficients[2, 4],
        N_observations_break_1 = nrow(segment1_data),
        Residual_SE_break_1 = summary1$sigma,
        Total_Incubation_Time_Min_break_1 = total_time_1,
        Number_Points_In_Respiration_Regression_break_1 = nrow(segment1_data),
        Number_Points_Removed_Respiration_Regression_break_1 = points_removed_1,
        DO_Concentration_At_Incubation_Time_Zero_break_1 = do_time_zero_1,
        Rate_mg_L_h_break_1 = rate1_mg_L_h,
        
        # Segment 2
        Slope_break_2 = slope2,
        Slope_break_2_SE = slope2_se,
        Intercept_break_2 = intercept2,
        R_squared_break_2 = summary2$r.squared,
        Adjusted_R_squared_break_2 = summary2$adj.r.squared,
        P_value_break_2 = summary2$coefficients[2, 4],
        N_observations_break_2 = nrow(segment2_data),
        Residual_SE_break_2 = summary2$sigma,
        Total_Incubation_Time_Min_break_2 = total_time_2,
        Number_Points_In_Respiration_Regression_break_2 = nrow(segment2_data),
        Number_Points_Removed_Respiration_Regression_break_2 = points_removed_2,
        DO_Concentration_At_Incubation_Time_Zero_break_2 = do_time_zero_2,
        Rate_mg_L_h_break_2 = rate2_mg_L_h,
        
        # Overall
        Breakpoint_Minutes = breakpoint_min,
        Segment1_Steeper_Than_Segment2 = segment1_steeper,
        Segmented_Successful = TRUE,
        
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, result_row)
      
      # Create plot
      if(plot_toggle) {
        p <- ggplot(sample_data, aes(x = relative_time_min_reset, y = DO_mg_per_L)) +
          geom_point(size = 2, alpha = 0.7) +
          geom_smooth(data = segment1_data, method = "lm", se = TRUE, 
                      color = "blue", aes(fill = "Segment 1")) +
          geom_smooth(data = segment2_data, method = "lm", se = TRUE, 
                      color = "red", aes(fill = "Segment 2")) +
          geom_vline(xintercept = breakpoint_min, linetype = "dashed", 
                     color = "black", alpha = 0.7) +
          theme_bw() +
          labs(title = paste(unique_samples[i]),
               subtitle = paste("Segment 1 rate:", round(rate1_mg_L_h, 3), "mg/L/h",
                                "| Segment 2 rate:", round(rate2_mg_L_h, 3), "mg/L/h"),
               x = "Time (minutes)",
               y = "DO (mg/L)",
               fill = "Regression Segments") +
          annotate("text", x = breakpoint_min, y = max(sample_data$DO_mg_per_L),
                   label = paste("Breakpoint:", round(breakpoint_min, 1), "min"),
                   vjust = -0.5, hjust = 0.5)
        
        print(p)
        
        # Save plot
        clean_name <- gsub("[^A-Za-z0-9_-]", "_", unique_samples[i])
        ggsave(paste0("04_Plots/", clean_name, "_segmented.pdf"),
               plot = p, width = 8, height = 6, device = "pdf")
      }
      
    } else {
      cat("  Insufficient points in one or both segments\n")
    }
    
  }, error = function(e) {
    cat("  Segmented regression failed:", e$message, "\n")
    # No special handling - just continue to next sample
  })
}

# ===== RESULTS SUMMARY =====

# Show summary
print(head(results[, c("Sample_Name", "Rate_mg_L_h_break_1", "Rate_mg_L_h_break_2", 
                       "Breakpoint_Minutes", "Segment1_Steeper_Than_Segment2")]))
results2 = results %>%
  dplyr::select(Sample_Name,
                Respiration_Rate_mg_DO_per_L_per_H = Rate_mg_L_h_break_1,
                Respiration_R_Squared = "R_squared_break_1",
                Respiration_R_Squared_Adj = "Adjusted_R_squared_break_1",
                Respiration_p_value = "P_value_break_1" ,
                Total_Incubation_Time_Min = Total_Incubation_Time_Min_break_1,
                Number_Points_In_Respiration_Regression = "Number_Points_In_Respiration_Regression_break_1" ,
                Number_Points_Removed_Respiration_Regression = "Number_Points_Removed_Respiration_Regression_break_1",
                DO_Concentration_At_Incubation_Time_Zero = "DO_Concentration_At_Incubation_Time_Zero_break_1") %>%
  mutate(Methods_Deviation = "N/A",
         dplyr::across(where(is.numeric), ~ round(.x, 2)))
# Save results
write.csv(results2, "./YEP_Respiration.csv", row.names = FALSE)
