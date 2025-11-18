# Note this code does not include Maggi's 2 minute point fast rule:
# If the 2 minute concentration is less than the Dissolved Oxygen Threshold set to remove low values from the back end of the curve, remove anything greater than 2 minutes
# And the deviation code RATE_006 rule:
## If samples were put on roller the same time as a picture was taken (indicated by RATE_006 Method Deviation) and the dissolved oxygen concentration is greater than the set threshold at the 2 minute measurement, remove the 2 minute point 
# Instead I added a rule that if DO goes down and up within the first 2 mins then remove the 1st 2 mins and start after with the regression line
# Clear workspace
rm(list=ls())

# Load libraries
library(dplyr)
library(ggplot2)
library(segmented)
library(broom)
library(lmtest)

# ===== User parameters =====
min_points <- 2
slope_thresh <- -0.006  
do_thresh <- 2
high_do <- 14
time_thresh_min <- 4
fast_rate <- 5.5
bp_pvalue_thresh <- 0.1
high_slope_thresh <- -0.04
time_same <- 7
break_toggle <- 1.4
conc_range <- 1.4
plot_toggle <- TRUE

# ===== Read data =====
data = read.csv('YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv', skip = 2) %>%
  filter(grepl("YEP", Sample_Name)) %>%
  dplyr::select(-c(Field_Name, IGSN, Material, Methods_Deviation))

unique_samples <- unique(data$Sample_Name)
print(paste("Found", length(unique_samples), "samples"))

# Create plots directory if it doesn't exist
if(plot_toggle) {
  if(!dir.exists("plots")) {
    dir.create("plots")
    cat("Created 'plots' directory\n")
  }
}

# ===== rates dataframe =====
results <- data.frame(
  Sample_Name = character(),
  slope_mg_L_per_sec = numeric(),
  rate_mg_L_per_min = numeric(),
  rate_mg_L_per_h = numeric(),
  R_squared = numeric(),
  R_squared_adj = numeric(),
  p_value = numeric(),
  total_time_min = numeric(),
  n_points_final = numeric(),
  n_points_original = numeric(),
  n_points_removed = numeric(),
  bp_pvalue = numeric(),
  start_DO = numeric(),
  end_DO = numeric(),
  used_segmented = logical(),
  early_fluctuation_removed = logical(),
  late_fluctuation_removed = logical(),
  stringsAsFactors = FALSE
)

# ===== Calculate rates =====
for(i in 1:length(unique_samples)) {
  cat("\n=== Processing sample", i, "of", length(unique_samples), ":", unique_samples[i], "===\n")
  
  # STEP 1: Get data for this sample
  sample_data <- data %>%
    filter(Sample_Name == unique_samples[i]) %>%
    arrange(Elapsed_Seconds) %>%
    mutate(across(c(Elapsed_Seconds, DO_mg_per_L), as.numeric))
  
  # Calculate relative time
  start_time <- min(sample_data$Elapsed_Seconds)
  sample_data$relative_time_min <- (sample_data$Elapsed_Seconds - start_time) / 60
  sample_data$relative_time_sec <- sample_data$Elapsed_Seconds - start_time
  
  original_n <- nrow(sample_data)
  cat("Step 1 - Original data:", original_n, "points\n")
  cat("Time range:", round(min(sample_data$relative_time_min), 2), "to",
      round(max(sample_data$relative_time_min), 2), "minutes\n")
  cat("DO range:", round(min(sample_data$DO_mg_per_L), 2), "to",
      round(max(sample_data$DO_mg_per_L), 2), "mg/L\n")
  
  # STEP 2: Remove high DO outliers first
  data_step2 <- sample_data %>%
    filter(DO_mg_per_L < high_do)
  cat("Step 2 - Remove DO >", high_do, "mg/L:", nrow(data_step2), "points remaining\n")
  
  # STEP 3: Remove low DO values after time threshold
  data_step3 <- data_step2 %>%
    filter(!(relative_time_min > time_thresh_min & DO_mg_per_L < do_thresh))
  cat("Step 3 - Remove DO <", do_thresh, "mg/L after", time_thresh_min, "min:", nrow(data_step3), "points remaining\n")
  
  # STEP 4: Check for EARLY fluctuation (first 2 minutes)
  early_fluctuation_removed <- FALSE
  first_2min_data <- data_step3 %>% filter(relative_time_min <= 2)
  
  if(nrow(first_2min_data) >= 3) {
    do_values <- first_2min_data$DO_mg_per_L
    has_fluctuation <- FALSE
    
    for(j in 2:(length(do_values)-1)) {
      if(do_values[j] < do_values[j-1] && do_values[j+1] > do_values[j]) {
        #if(do_values[j+1] > do_values[j]) {
        has_fluctuation <- TRUE
        print(j)
        cat("  Detected EARLY DO fluctuation: ", round(do_values[j-1], 3), "→",
            round(do_values[j], 3), "→", round(do_values[j+1], 3), "mg/L\n")
        break
      }
    }
    
    if(has_fluctuation) {
      data_step3 <- data_step3 %>% filter(relative_time_min > 2)
      early_fluctuation_removed <- TRUE
      cat("  Removed first 2 minutes due to fluctuation\n")
    }
  }
  
  cat("Step 4 - After early fluctuation check:", nrow(data_step3), "points remaining\n")
  
  # STEP 5: NEW - Check for LATE increasing DO trend (last 2 minutes)
  late_fluctuation_removed <- FALSE
  
  if(nrow(data_step3) >= 3) {
    max_time <- max(data_step3$relative_time_min)
    last_2min_data <- data_step3 %>% 
      filter(relative_time_min >= (max_time - 2)) %>%
      arrange(relative_time_min)
    
    if(nrow(last_2min_data) >= 3) {
      do_values_end <- last_2min_data$DO_mg_per_L
      has_increasing_trend <- FALSE
      
      # Check for consecutive increases in DO (bad for respiration measurements)
      consecutive_increases <- 0
      for(j in 2:length(do_values_end)) {
        if(do_values_end[j] > do_values_end[j-1]) {
          consecutive_increases <- consecutive_increases + 1
        } else {
          consecutive_increases <- 0  # Reset if decrease is found
        }
        
        # If we find 2 or more consecutive increases, flag as problematic
        if(consecutive_increases >= 2) {
          has_increasing_trend <- TRUE
          cat("  Detected increasing DO trend in last 2 minutes: ", 
              paste(round(tail(do_values_end, 3), 3), collapse = "→"), "mg/L\n")
          break
        }
      }
      
      if(has_increasing_trend) {
        data_step3 <- data_step3 %>% filter(relative_time_min < (max_time - 2))
        late_fluctuation_removed <- TRUE
        cat("  Removed last 2 minutes due to increasing DO trend\n")
      }
    }
  }
  
  cat("Step 5 - After late increasing trend check:", nrow(data_step3), "points remaining\n")
  
  # STEP 6: Ensure minimum points
  if(nrow(data_step3) < min_points) {
    cat("  WARNING: Less than", min_points, "points remaining. Skipping sample.\n")
    next
  }
  
  # STEP 7: Fit initial linear model
  fit_initial <- lm(DO_mg_per_L ~ relative_time_sec, data = data_step3)
  initial_slope <- coef(fit_initial)[2]
  initial_r2 <- summary(fit_initial)$r.squared
  cat("Step 7 - Initial fit: slope =", round(initial_slope * 3600, 3), "mg/L/hour, R² =", round(initial_r2, 3), "\n")
  
  # Initialize final results with initial fit
  final_data <- data_step3
  final_fit <- fit_initial
  used_segmented <- FALSE
  
  # STEP 8: Segmented Regression Analysis
  do_range <- first(data_step3$DO_mg_per_L) - last(data_step3$DO_mg_per_L)
  
  if(nrow(data_step3) > 3 && initial_slope < slope_thresh && do_range > conc_range) {
    cat("  Trying segmented regression (slope=", round(initial_slope*3600, 3), ", range=", round(do_range, 2), ")\n")
    
    tryCatch({
      mid_time <- (max(data_step3$relative_time_sec) + min(data_step3$relative_time_sec)) / 2
      seg_fit <- segmented(fit_initial, seg.Z = ~relative_time_sec, psi = mid_time)
      breakpoint <- seg_fit$psi[2]
      
      data_segmented <- data_step3 %>%
        filter(relative_time_sec <= breakpoint)
      
      if(nrow(data_segmented) >= min_points) {
        fit_segmented <- lm(DO_mg_per_L ~ relative_time_sec, data = data_segmented)
        slope_ratio <- coef(fit_segmented)[2] / coef(fit_initial)[2]
        
        cat("    Breakpoint:", round(breakpoint/60, 2), "min, Slope ratio:", round(slope_ratio, 2), "\n")
        
        if(slope_ratio > break_toggle) {
          final_data <- data_segmented
          final_fit <- fit_segmented
          used_segmented <- TRUE
          cat("    Using segmented regression\n")
        }
      }
    }, error = function(e) {
      cat("    Segmented regression failed:", e$message, "\n")
    })
  }
  
  # STEP 9: Breusch-Pagan Test and Point Removal
  if(!used_segmented && initial_slope < slope_thresh && do_range > conc_range) {
    cat("  Applying Breusch-Pagan test for heteroscedasticity\n")
    
    current_data <- final_data
    for(iter in 1:60) {
      if(nrow(current_data) <= min_points) break
      
      current_fit <- lm(DO_mg_per_L ~ relative_time_sec, data = current_data)
      bp_test <- bptest(current_fit)
      bp_pval <- bp_test$p.value
      
      if(bp_pval >= bp_pvalue_thresh) {
        break
      }
      current_data <- current_data[-nrow(current_data), ]
    }
    
    if(nrow(current_data) >= min_points) {
      final_data <- current_data
      final_fit <- lm(DO_mg_per_L ~ relative_time_sec, data = final_data)
      cat("    BP test iterations:", iter, ", final p-value:", round(bp_test$p.value, 3), "\n")
    }
  }
  
  # Final check for positive slope
  if(coef(final_fit)[2] > 0 && initial_slope < 0) {
    cat("  WARNING: Final slope is positive, reverting to initial data\n")
    final_data <- data_step3
    final_fit <- fit_initial
    used_segmented <- FALSE
  }
  
  cat("Final analysis: slope =", round(coef(final_fit)[2] * 3600, 3),
      "mg/L/hour, R² =", round(summary(final_fit)$r.squared, 3), "\n")
  
  # Get final BP test result
  final_bp <- tryCatch({
    bptest(final_fit)$p.value
  }, error = function(e) NA)
  
  # Store results
  results <- rbind(results, data.frame(
    Sample_Name = unique_samples[i],
    slope_mg_L_per_sec = coef(final_fit)[2],
    rate_mg_L_per_min = abs(coef(final_fit)[2]) * 60,
    rate_mg_L_per_h = abs(coef(final_fit)[2]) * 3600,
    R_squared = summary(final_fit)$r.squared,
    R_squared_adj = summary(final_fit)$adj.r.squared,
    p_value = summary(final_fit)$coefficients[2,4],
    total_time_min = max(final_data$relative_time_min) - min(final_data$relative_time_min),
    n_points_final = nrow(final_data),
    n_points_original = original_n,
    n_points_removed = original_n - nrow(final_data),
    bp_pvalue = final_bp,
    start_DO = final_data$DO_mg_per_L[1],
    end_DO = final_data$DO_mg_per_L[nrow(final_data)],
    used_segmented = used_segmented,
    early_fluctuation_removed = early_fluctuation_removed,
    late_fluctuation_removed = late_fluctuation_removed,
    stringsAsFactors = FALSE
  ))
  
  # Create and save plots as PDFs
  if(plot_toggle) {
    p <- ggplot(final_data, aes(x = relative_time_min, y = DO_mg_per_L)) +
      geom_point(size = 2, color = if(used_segmented) "red" else "blue") +
      geom_smooth(method = "lm", se = TRUE) +
      labs(title = paste(unique_samples[i]),
           subtitle = paste("Rate:", round(abs(coef(final_fit)[2]) * 3600, 3), "mg/L/hour, R² =",
                            round(summary(final_fit)$r.squared, 3),
                            if(used_segmented) " (Segmented)" else "",
                            if(early_fluctuation_removed) " (Early fluct. removed)" else "",
                            if(late_fluctuation_removed) " (Late fluct. removed)" else ""),
           x = "Time (minutes from start)",
           y = "DO (mg/L)") +
      theme_bw()
    
    # Display the plot
    print(p)
    
    # Clean the sample name to remove characters that aren't allowed in filenames
    clean_name <- gsub("[^A-Za-z0-9_-]", "_", unique_samples[i])
    
    # Save as PDF in plots folder
    ggsave(filename = paste0("plots/", clean_name, "_respiration_plot.pdf"),
           plot = p,
           width = 8, 
           height = 6, 
           units = "in",
           device = "pdf")
    
    cat("  Plot saved as:", paste0("plots/", clean_name, "_respiration_plot.pdf"), "\n")
  }
}

# ===== FINAL RESULTS =====
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Processed", nrow(results), "samples successfully\n")
cat("Samples with early fluctuation removed:", sum(results$early_fluctuation_removed), "\n")
cat("Samples with late fluctuation removed:", sum(results$late_fluctuation_removed), "\n")
cat("Samples using segmented regression:", sum(results$used_segmented), "\n")

print(head(results))
cat("\nSummary of rates (mg/L/h):\n")
print(summary(results$rate_mg_L_per_h))

# Save results to CSV
write.csv(results, "sediment_respiration_results.csv", row.names = FALSE)
cat("Results saved to: sediment_respiration_results.csv\n")