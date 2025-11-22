# Clear workspace
rm(list=ls())

# Load libraries
library(dplyr)
library(ggplot2)
library(segmented)
library(broom)
library(lmtest)
library(reshape2)

# ===== User parameters =====
min_points <- 2
slope_thresh <- -0.36 #in mg/L/h or -0.006 in mg/L/min  
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

# ===== EXPANDED results dataframe for comparison =====
results <- data.frame(
  Sample_Name = character(),
  
  # Original Linear Regression Results
  original_slope_mg_L_per_sec = numeric(),
  original_rate_mg_L_per_h = numeric(),
  original_R_squared = numeric(),
  original_R_squared_adj = numeric(),
  original_p_value = numeric(),
  original_AIC = numeric(),
  original_BIC = numeric(),
  original_bp_pvalue = numeric(),
  original_n_points = numeric(),
  
  # Segmented Regression Results
  segmented_slope_mg_L_per_sec = numeric(),
  segmented_rate_mg_L_per_h = numeric(),
  segmented_R_squared = numeric(),
  segmented_R_squared_adj = numeric(),
  segmented_p_value = numeric(),
  segmented_AIC = numeric(),
  segmented_BIC = numeric(),
  segmented_bp_pvalue = numeric(),
  segmented_n_points = numeric(),
  breakpoint_min = numeric(),
  segmented_successful = logical(),
  
  # Comparison Metrics
  better_method = character(),  # "original" or "segmented" 
  R2_improvement = numeric(),   # segmented R² - original R²
  AIC_improvement = numeric(),  # original AIC - segmented AIC (positive = segmented better)
  
  # Data Processing Info
  total_time_min = numeric(),
  n_points_original = numeric(),
  n_points_removed = numeric(),
  start_DO = numeric(),
  end_DO = numeric(),
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
        has_fluctuation <- TRUE
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
  
  # STEP 5: Check for LATE increasing DO trend (last 2 minutes)
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
  
  # STEP 7: ALWAYS fit original linear model
  fit_original <- lm(DO_mg_per_L ~ relative_time_sec, data = data_step3)
  
  # Calculate original model statistics
  original_slope <- coef(fit_original)[2]
  original_r2 <- summary(fit_original)$r.squared
  original_r2_adj <- summary(fit_original)$adj.r.squared
  original_p <- summary(fit_original)$coefficients[2,4]
  original_aic <- AIC(fit_original)
  original_bic <- BIC(fit_original)
  original_bp <- tryCatch({bptest(fit_original)$p.value}, error = function(e) NA)
  
  cat("Step 7 - Original fit: slope =", round(original_slope * 3600, 3), "mg/L/hour, R² =", round(original_r2, 3), "\n")
  
  # STEP 8: ALWAYS attempt segmented regression (when >3 points)
  segmented_successful <- FALSE
  segmented_slope <- NA
  segmented_r2 <- NA
  segmented_r2_adj <- NA
  segmented_p <- NA
  segmented_aic <- NA
  segmented_bic <- NA
  segmented_bp <- NA
  segmented_n_points <- NA
  breakpoint_min <- NA
  fit_segmented <- NULL
  data_segmented <- NULL
  
  if(nrow(data_step3) > 3) {
    cat("  Attempting segmented regression...\n")
    
    tryCatch({
      # Try segmented regression
      mid_time <- (max(data_step3$relative_time_sec) + min(data_step3$relative_time_sec)) / 2
      seg_fit <- segmented(fit_original, seg.Z = ~relative_time_sec, psi = mid_time)
      
      breakpoint <- seg_fit$psi[2]
      breakpoint_min <- breakpoint / 60
      
      # Get data up to breakpoint
      data_segmented <- data_step3 %>% filter(relative_time_sec <= breakpoint)
      
      if(nrow(data_segmented) >= min_points) {
        fit_segmented <- lm(DO_mg_per_L ~ relative_time_sec, data = data_segmented)
        
        # Calculate segmented model statistics
        segmented_slope <- coef(fit_segmented)[2]
        segmented_r2 <- summary(fit_segmented)$r.squared
        segmented_r2_adj <- summary(fit_segmented)$adj.r.squared
        segmented_p <- summary(fit_segmented)$coefficients[2,4]
        segmented_aic <- AIC(fit_segmented)
        segmented_bic <- BIC(fit_segmented)
        segmented_bp <- tryCatch({bptest(fit_segmented)$p.value}, error = function(e) NA)
        segmented_n_points <- nrow(data_segmented)
        segmented_successful <- TRUE
        
        cat("    Segmented fit: breakpoint =", round(breakpoint_min, 2), "min, slope =", 
            round(segmented_slope * 3600, 3), "mg/L/hour, R² =", round(segmented_r2, 3), "\n")
      }
    }, error = function(e) {
      cat("    Segmented regression failed:", e$message, "\n")
    })
  }
  
  # STEP 9: Improved comparison and decision criteria
  better_method <- "original"  # default
  r2_improvement <- 0
  aic_improvement <- 0
  
  if(segmented_successful) {
    r2_improvement <- segmented_r2 - original_r2
    aic_improvement <- original_aic - segmented_aic
    
    # IMPROVED DECISION CRITERIA
    segmented_is_better <- (
      segmented_r2 > original_r2 &&                    # Better R²
        aic_improvement > 2 &&                           # Meaningful AIC improvement
        r2_improvement > 0.02 &&                         # Substantial R² improvement
        segmented_p < 0.05 &&                            # Significant slope
        original_p < 0.05 &&                             # Original also significant
        segmented_n_points >= (0.5 * nrow(data_step3)) && # Uses >50% of data
        breakpoint_min > 2 && breakpoint_min < (max(data_step3$relative_time_min) - 2) # Reasonable breakpoint
    )
    
    if(segmented_is_better) {
      better_method <- "segmented"
    }
    
    cat("  Comparison: R² improvement =", round(r2_improvement, 4), 
        ", AIC improvement =", round(aic_improvement, 2), "\n")
    cat("  Breakpoint at:", round(breakpoint_min, 2), "min, using", 
        round(100 * segmented_n_points / nrow(data_step3), 1), "% of data\n")
    cat("  Better method:", better_method, "\n")
  }
  
  # Store comprehensive results
  results <- rbind(results, data.frame(
    Sample_Name = unique_samples[i],
    
    # Original results
    original_slope_mg_L_per_sec = original_slope,
    original_rate_mg_L_per_h = abs(original_slope) * 3600,
    original_R_squared = original_r2,
    original_R_squared_adj = original_r2_adj,
    original_p_value = original_p,
    original_AIC = original_aic,
    original_BIC = original_bic,
    original_bp_pvalue = original_bp,
    original_n_points = nrow(data_step3),
    
    # Segmented results
    segmented_slope_mg_L_per_sec = ifelse(segmented_successful, segmented_slope, NA),
    segmented_rate_mg_L_per_h = ifelse(segmented_successful, abs(segmented_slope) * 3600, NA),
    segmented_R_squared = ifelse(segmented_successful, segmented_r2, NA),
    segmented_R_squared_adj = ifelse(segmented_successful, segmented_r2_adj, NA),
    segmented_p_value = ifelse(segmented_successful, segmented_p, NA),
    segmented_AIC = ifelse(segmented_successful, segmented_aic, NA),
    segmented_BIC = ifelse(segmented_successful, segmented_bic, NA),
    segmented_bp_pvalue = ifelse(segmented_successful, segmented_bp, NA),
    segmented_n_points = ifelse(segmented_successful, segmented_n_points, NA),
    breakpoint_min = breakpoint_min,
    segmented_successful = segmented_successful,
    
    # Comparisons
    better_method = better_method,
    R2_improvement = r2_improvement,
    AIC_improvement = aic_improvement,
    
    # Other info
    total_time_min = max(data_step3$relative_time_min) - min(data_step3$relative_time_min),
    n_points_original = original_n,
    n_points_removed = original_n - nrow(data_step3),
    start_DO = data_step3$DO_mg_per_L[1],
    end_DO = data_step3$DO_mg_per_L[nrow(data_step3)],
    early_fluctuation_removed = early_fluctuation_removed,
    late_fluctuation_removed = late_fluctuation_removed,
    
    stringsAsFactors = FALSE
  ))
  
  # Enhanced plotting with both rates displayed
  if(plot_toggle) {
    # Calculate rates for display
    original_rate <- abs(original_slope) * 3600
    segmented_rate_display <- ifelse(segmented_successful, abs(segmented_slope) * 3600, NA)
    
    # Create comparison plot
    p <- ggplot(data_step3, aes(x = relative_time_min, y = DO_mg_per_L)) +
      geom_point(size = 2.5, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", 
                  alpha = 0.3, aes(fill = "Original Linear")) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(title = paste("Sample:", unique_samples[i]),
           subtitle = paste(
             "Original: Rate =", round(original_rate, 3), "mg/L/h, R² =", round(original_r2, 3),
             ifelse(segmented_successful, 
                    paste("\nSegmented: Rate =", round(segmented_rate_display, 3), 
                          "mg/L/h, R² =", round(segmented_r2, 3),
                          "\nSelected Method:", toupper(better_method)), 
                    "\nSegmented: FAILED")),
           x = "Time (minutes from start)",
           y = "DO (mg/L)",
           fill = "Regression Type")
    
    # Add segmented fit if successful
    if(segmented_successful && !is.null(data_segmented)) {
      p <- p + 
        geom_smooth(data = data_segmented, method = "lm", se = TRUE, 
                    color = "red", linetype = "dashed", alpha = 0.3,
                    aes(fill = "Segmented Linear")) +
        geom_vline(xintercept = breakpoint_min, color = "red", 
                   linetype = "dotted", alpha = 0.7, size = 1) +
        annotate("text", x = breakpoint_min + 0.5, 
                 y = max(data_step3$DO_mg_per_L) * 0.9, 
                 label = paste("Breakpoint:", round(breakpoint_min, 1), "min"), 
                 color = "red", size = 3, hjust = 0)
    }
    
    # Highlight the chosen method
    if(better_method == "segmented" && segmented_successful) {
      p <- p + 
        ggtitle(paste("Sample:", unique_samples[i], "- SEGMENTED SELECTED"),
                subtitle = paste(
                  "Original: Rate =", round(original_rate, 3), "mg/L/h, R² =", round(original_r2, 3),
                  "\nSegmented: Rate =", round(segmented_rate_display, 3), 
                  "mg/L/h, R² =", round(segmented_r2, 3),
                  "\nImprovement: ΔR² =", round(r2_improvement, 3), 
                  ", ΔAIC =", round(aic_improvement, 1))) +
        theme(plot.title = element_text(color = "red", face = "bold"))
    } else {
      p <- p + 
        theme(plot.title = element_text(color = "blue", face = "bold"))
    }
    
    print(p)
    
    # Save plot
    clean_name <- gsub("[^A-Za-z0-9_-]", "_", unique_samples[i])
    ggsave(filename = paste0("plots/", clean_name, "_comparison_plot.pdf"),
           plot = p, width = 10, height = 7, units = "in", device = "pdf")
    cat("  Plot saved as:", paste0("plots/", clean_name, "_comparison_plot.pdf"), "\n")
  }
}

# ===== FINAL RESULTS AND SUMMARY =====
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Processed", nrow(results), "samples successfully\n")
cat("Segmented regression successful:", sum(results$segmented_successful), "samples\n")
cat("Samples where segmented is better:", sum(results$better_method == "segmented", na.rm = TRUE), "samples\n")
cat("Samples where original is better:", sum(results$better_method == "original", na.rm = TRUE), "samples\n")
cat("Samples with early fluctuation removed:", sum(results$early_fluctuation_removed), "\n")
cat("Samples with late fluctuation removed:", sum(results$late_fluctuation_removed), "\n")

# Summary statistics
cat("\n=== COMPARISON SUMMARY ===\n")
cat("Mean R² improvement from segmented:", round(mean(results$R2_improvement, na.rm = TRUE), 4), "\n")
cat("Mean AIC improvement from segmented:", round(mean(results$AIC_improvement, na.rm = TRUE), 2), "\n")

# Show first few results
cat("\n=== SAMPLE RESULTS ===\n")
print(head(results[, c("Sample_Name", "original_R_squared", "segmented_R_squared", 
                       "R2_improvement", "AIC_improvement", "better_method")]))

# Save results to CSV
write.csv(results, "sediment_respiration_comparison_results.csv", row.names = FALSE)
cat("Results saved to: sediment_respiration_comparison_results.csv\n")

# Create summary comparison plot
if(plot_toggle && nrow(results) > 0) {
  comparison_data <- results %>%
    filter(segmented_successful) %>%
    dplyr::select(Sample_Name, original_R_squared, segmented_R_squared) %>%
    melt(id.vars = "Sample_Name", variable.name = "Method", value.name = "R_squared")
  
  if(nrow(comparison_data) > 0) {
    p_summary <- ggplot(comparison_data, aes(x = Sample_Name, y = R_squared, fill = Method)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "R² Comparison: Original vs Segmented Regression",
           x = "Sample", y = "R²") +
      scale_fill_manual(values = c("original_R_squared" = "blue", "segmented_R_squared" = "red"),
                        labels = c("Original Linear", "Segmented Linear"))
    
    ggsave("plots/R2_comparison_summary.pdf", p_summary, width = 12, height = 6, device = "pdf")
    cat("Summary comparison plot saved to: plots/R2_comparison_summary.pdf\n")
    
    # Also create a rate comparison plot
    rate_comparison_data <- results %>%
      filter(segmented_successful) %>%
      dplyr::select(Sample_Name, original_rate_mg_L_per_h, segmented_rate_mg_L_per_h) %>%
      melt(id.vars = "Sample_Name", variable.name = "Method", value.name = "Rate_mg_L_per_h")
    
    p_rates <- ggplot(rate_comparison_data, aes(x = Sample_Name, y = Rate_mg_L_per_h, fill = Method)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Rate Comparison: Original vs Segmented Regression",
           x = "Sample", y = "Rate (mg/L/h)") +
      scale_fill_manual(values = c("original_rate_mg_L_per_h" = "blue", "segmented_rate_mg_L_per_h" = "red"),
                        labels = c("Original Linear", "Segmented Linear"))
    
    ggsave("plots/Rate_comparison_summary.pdf", p_rates, width = 12, height = 6, device = "pdf")
    cat("Rate comparison plot saved to: plots/Rate_comparison_summary.pdf\n")
  }
}

