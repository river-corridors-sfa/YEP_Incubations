# ----------------------------------------
# Simplified DO time series fitting - Linear and Combined only
# ----------------------------------------
rm(list=ls(all=T))
library(deSolve)   
library(FME)       
library(ggplot2)   
library(dplyr)
library(tidyr)

# -------------------------------
# Load data
# -------------------------------
file_path <- "C:/Users/gara009/PNNL/RC-SFA - Documents/Study_YEP/PRELIMINARY_YEP_Data_Package/YEP_Sample_Data/YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv"

mapping_file_path <- "C:/Users/gara009/PNNL/RC-SFA - Documents/Study_YEP/DO/05_PublishReadyData/Merged_Firesting_Mapping_2025-10-24.csv"

mapping_df <- read.csv(mapping_file_path, stringsAsFactors = FALSE)


df <- read.csv(file_path, skip = 14, header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("Field_Name","Sample_Name","IGSN","Material","DateTime",
                  "Elapsed_Seconds","Temperature_degreesC","DO_mg_per_L",
                  "Firesting_Serial_Number","Methods_Deviation")


df <- df %>%
  inner_join(mapping_df %>% select(Sample_Name, Start_Time_Elapsed_Seconds = Elapsed_Seconds_Start, End_Time_Elapsed_Seconds = Elapsed_Seconds_End), 
             by = "Sample_Name") %>%
  filter(Elapsed_Seconds >= Start_Time_Elapsed_Seconds, 
         Elapsed_Seconds <= End_Time_Elapsed_Seconds) %>%
  select(-Start_Time_Elapsed_Seconds, -End_Time_Elapsed_Seconds)

df <- df %>%
  mutate(Elapsed_Seconds = as.numeric(Elapsed_Seconds),
         DO_mg_per_L     = as.numeric(DO_mg_per_L)) %>%
  filter(!is.na(Elapsed_Seconds), !is.na(DO_mg_per_L)) %>%
  filter(DO_mg_per_L > 0)

# -------------------------------
# Constants
# -------------------------------
K_O2 <- 10
BM10 <- 1.0

# -------------------------------
# Helper functions
# -------------------------------
compute_info <- function(SSR, k, n) {
  AIC <- n * log(SSR / n) + 2 * k  # Fixed multiplication operators
  BIC <- n * log(SSR / n) + k * log(n)  # Fixed multiplication operators
  list(AIC = AIC, BIC = BIC)
}

# -------------------------------
# Combined model solver (only one needed)
# -------------------------------
solveCOM <- function(pars, times, Oxy0, BM10, K_O2) {
  derivs <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      r1   <- vmax1 * Oxy/(Oxy + K_O2 + 1e-10)
      dOxy <- - r1 * BM1 - va * Oxy  # Fixed multiplication operators
      dBM1 <-   r1 * BM1
      list(c(dOxy, dBM1))
    })
  }
  state <- c(Oxy = Oxy0, BM1 = BM10)
  
  tryCatch({
    out <- ode(y = state, times = times, func = derivs, parms = pars,
               rtol = 1e-4, atol = 1e-6, method = "lsoda")
    result <- as.data.frame(out[, c("time","Oxy")])
    
    if(any(is.na(result$Oxy)) || any(result$Oxy < 0) || any(result$Oxy > 2*Oxy0)) {
      stop("Invalid solution")
    }
    return(result)
  }, error = function(e) {
    # Analytical approximation for small times
    total_rate <- pars$va + pars$vmax1 * (Oxy0/(Oxy0 + K_O2)) * BM10  # Fixed multiplication
    result <- data.frame(
      time = times,
      Oxy = pmax(0, Oxy0 * exp(-total_rate * times))  # Fixed multiplication
    )
    return(result)
  })
}

# -------------------------------
# Simplified summary function
# -------------------------------
create_simple_summary <- function(Data, fit_COM, model_COM, lm_DO, sample_name, n, Oxy0) {
  # Helper function to calculate R-squared and RMSE
  calc_fit_metrics <- function(observed, predicted) {
    residuals <- observed - predicted
    ss_res <- sum(residuals^2)
    ss_tot <- sum((observed - mean(observed))^2)
    r_squared <- 1 - (ss_res / ss_tot)
    rmse <- sqrt(mean(residuals^2))
    list(r_squared = r_squared, rmse = rmse)
  }
  
  # Initialize summary
  summary_row <- data.frame(
    Sample_Name = sample_name,
    n_points = n,
    Initial_DO_mg_per_L = Oxy0,
    
    # Linear model
    Linear_Rate_mg_DO_per_L_per_H = NA,
    Linear_AIC = NA,
    Linear_R_squared = NA,
    Linear_RMSE_mg_per_L = NA,
    
    # Combined model
    Combined_vmax1_mg_per_L_per_H = NA,
    Combined_va_per_H = NA,
    Combined_Initial_Rate_mg_DO_per_L_per_H = NA,
    Combined_AIC = NA,
    Combined_R_squared = NA,
    Combined_RMSE_mg_per_L = NA,
    
    # Best model
    Best_Model = NA,
    Best_Rate_mg_DO_per_L_per_H = NA
  )
  
  # Linear model metrics
  if(!is.null(lm_DO)) {
    b <- coef(lm_DO)[["time"]]
    summary_row$Linear_Rate_mg_DO_per_L_per_H <- -b
    
    SSR_lm <- sum(residuals(lm_DO)^2)
    info_lm <- compute_info(SSR_lm, k = 2, n = n)
    summary_row$Linear_AIC <- info_lm$AIC
    
    lm_metrics <- calc_fit_metrics(Data$Oxy, fitted(lm_DO))
    summary_row$Linear_R_squared <- lm_metrics$r_squared
    summary_row$Linear_RMSE_mg_per_L <- lm_metrics$rmse
  }
  
  # Combined model metrics
  if(!is.null(fit_COM) && !is.null(model_COM)) {
    vmax1 <- fit_COM$par[["vmax1"]]
    va <- fit_COM$par[["va"]]
    summary_row$Combined_vmax1_mg_per_L_per_H <- vmax1
    summary_row$Combined_va_per_H <- va
    
    bio_rate <- vmax1 * (Oxy0/(Oxy0 + K_O2)) * BM10  # Fixed multiplication
    abi_rate <- va * Oxy0
    summary_row$Combined_Initial_Rate_mg_DO_per_L_per_H <- bio_rate + abi_rate
    
    info_com <- compute_info(fit_COM$ssr, k = 2, n = n)
    summary_row$Combined_AIC <- info_com$AIC
    
    com_metrics <- calc_fit_metrics(Data$Oxy, model_COM$Oxy)
    summary_row$Combined_R_squared <- com_metrics$r_squared
    summary_row$Combined_RMSE_mg_per_L <- com_metrics$rmse
  }
  
  # Determine best model by AIC
  aic_values <- c()
  if(!is.na(summary_row$Linear_AIC)) aic_values <- c(aic_values, Linear = summary_row$Linear_AIC)
  if(!is.na(summary_row$Combined_AIC)) aic_values <- c(aic_values, Combined = summary_row$Combined_AIC)
  
  if(length(aic_values) > 0) {
    best_model <- names(aic_values)[which.min(aic_values)]
    summary_row$Best_Model <- best_model
    
    if(best_model == "Linear") {
      summary_row$Best_Rate_mg_DO_per_L_per_H <- summary_row$Linear_Rate_mg_DO_per_L_per_H
    } else if(best_model == "Combined") {
      summary_row$Best_Rate_mg_DO_per_L_per_H <- summary_row$Combined_Initial_Rate_mg_DO_per_L_per_H
    }
  }
  
  return(summary_row)
}

# -------------------------------
# Simplified fitting function
# -------------------------------
fit_two_models <- function(dat, sample_name, out_dir = ".") {
  # Clean data
  Data <- dat %>%
    arrange(time_hr) %>%
    select(time = time_hr, Oxy = DO_mg_per_L) %>%
    group_by(time) %>%
    summarise(Oxy = mean(Oxy, na.rm = TRUE), .groups = 'drop') %>%
    filter(Oxy > 0, !is.na(Oxy), !is.na(time)) %>%
    arrange(time) %>%
    as.data.frame()
  
  n <- nrow(Data)
  if (n < 3) {
    message("Skipping ", sample_name, ": not enough points.")
    return(NULL)
  }
  
  Oxy0 <- Data$Oxy[1]
  
  # Estimate reasonable parameter ranges from data
  linear_rate <- abs(coef(lm(Oxy ~ time, data = Data))[["time"]])
  
  message("Fitting models for ", sample_name, " (n=", n, ", linear rate~", round(linear_rate, 2), " mg/L/h)")
  
  # Combined model objective
  Objective_COM <- function(x) {
    if(x[["vmax1"]] <= 0 || x[["vmax1"]] > 100 || x[["va"]] <= 0 || x[["va"]] > 50) return(1e10)
    pars <- list(vmax1 = x[["vmax1"]], va = x[["va"]])
    out <- solveCOM(pars, times = Data$time, Oxy0 = Oxy0, BM10 = BM10, K_O2 = K_O2)
    if(any(is.na(out$Oxy))) return(1e10)
    modCost(model = out, obs = Data)
  }
  
  # Fit Combined model
  fit_COM <- NULL
  try({
    fit_COM <- modFit(p = c(vmax1 = linear_rate/2, va = linear_rate/Oxy0/2), f = Objective_COM,
                      lower = c(vmax1 = 0.1, va = 0.01), upper = c(vmax1 = 50, va = 10))
    message("  Combined fit successful: vmax1=", round(fit_COM$par[["vmax1"]], 3), 
            ", va=", round(fit_COM$par[["va"]], 3))
  }, silent = TRUE)
  
  # Generate predictions
  model_COM <- NULL
  if(!is.null(fit_COM)) {
    pars_COM <- list(vmax1 = fit_COM$par[["vmax1"]], va = fit_COM$par[["va"]])
    model_COM <- solveCOM(pars_COM, times = Data$time, Oxy0 = Oxy0, BM10 = BM10, K_O2 = K_O2)
  }
  
  # Linear model
  lm_DO <- lm(Oxy ~ time, data = Data)
  R_lin <- -coef(lm_DO)[["time"]]
  message("  Linear fit: rate=", round(R_lin, 3), " mg/L/h")
  
  # Create simplified plots with rate values
  tryCatch({
    pdf_file <- file.path(out_dir, paste0("DO_fits_simple_", gsub("[^A-Za-z0-9_-]", "_", sample_name), ".pdf"))
    pdf(pdf_file, width = 8, height = 6)
    
    # Single plot with both models
    plot(Data$time, Data$Oxy, pch = 16, col = "black", cex = 1.2,
         xlab = "Time (hours)", ylab = "DO (mg O2 / L)",
         main = paste("Linear vs Combined model:", sample_name))
    
    # Plot linear model
    lines(Data$time, fitted(lm_DO), col = "blue", lwd = 3, lty = 1)
    
    # Plot combined model if successful
    if(!is.null(model_COM)) {
      lines(model_COM$time, model_COM$Oxy, col = "red", lwd = 3, lty = 2)
    }
    
    # Create rate text instead of AIC text
    rate_text <- "Respiration Rates:\n"
    rate_text <- paste0(rate_text, sprintf("Linear: %.2f mg/L/h\n", R_lin))
    
    if(!is.null(fit_COM)) {
      # Calculate combined initial rate
      vmax1 <- fit_COM$par[["vmax1"]]
      va <- fit_COM$par[["va"]]
      bio_rate <- vmax1 * (Oxy0/(Oxy0 + K_O2)) * BM10
      abi_rate <- va * Oxy0
      combined_rate <- bio_rate + abi_rate
      rate_text <- paste0(rate_text, sprintf("Combined: %.2f mg/L/h\n", combined_rate))
      rate_text <- paste0(rate_text, sprintf("  (Bio: %.2f, Abio: %.2f)", bio_rate, abi_rate))
    }
    
    # Position rate text in middle-left of plot
    text(x = max(Data$time) * 0.05, y = mean(range(Data$Oxy)),
         labels = rate_text, adj = c(0, 0.5),
         bg = "white", cex = 0.9,
         col = "black", family = "mono")
    
    # Legend
    legend_labels <- c("Observed", "Linear")
    legend_colors <- c("black", "blue")
    legend_pch <- c(16, NA)
    legend_lty <- c(NA, 1)
    legend_lwd <- c(NA, 3)
    
    if(!is.null(model_COM)) {
      legend_labels <- c(legend_labels, "Combined")
      legend_colors <- c(legend_colors, "red")
      legend_pch <- c(legend_pch, NA)
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 3)
    }
    
    legend("topright", legend_labels, 
           pch = legend_pch, lty = legend_lty, 
           col = legend_colors, lwd = legend_lwd, 
           bty = "n", cex = 1.0)
    
    dev.off()
    message("Wrote: ", pdf_file)
  }, error = function(e) {
    message("Error creating plots for ", sample_name, ": ", e$message)
  })
  
  # Create summary
  summary_row <- create_simple_summary(Data, fit_COM, model_COM, lm_DO, sample_name, n, Oxy0)
  
  list(summary = summary_row)
}

# -------------------------------
# Run analysis for all reactors
# -------------------------------
# Rebase time per Sample_Name and convert to hours
df2 <- df %>%
  group_by(Sample_Name) %>%
  arrange(Elapsed_Seconds, .by_group = TRUE) %>%
  mutate(time_hr = (Elapsed_Seconds - min(Elapsed_Seconds)) / 3600) %>%
  ungroup()

samples <- unique(df2$Sample_Name)
out_dir <- "."

message("Starting simplified analysis for ", length(samples), " samples...")

all_summaries <- lapply(samples, function(sname) {
  dat_s <- df2 %>% filter(Sample_Name == sname)
  fit_two_models(dat_s, sname, out_dir = out_dir)
})

# Combine summaries
summary_df <- bind_rows(lapply(all_summaries, function(x) if(!is.null(x)) x$summary else NULL))
summary_file <- file.path(out_dir, "DO_simple_model_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
message("Wrote summary CSV: ", summary_file)

# Print summary
if(nrow(summary_df) > 0) {
  message("\n=== SIMPLIFIED MODEL PERFORMANCE SUMMARY ===")
  
  # Show model preference
  model_summary <- summary_df %>%
    count(Best_Model, sort = TRUE) %>%
    mutate(Percentage = round(n / sum(n) * 100, 1))
  
  print(model_summary)
  
  # Show rate statistics
  message("\n=== RESPIRATION RATE STATISTICS ===")
  rate_stats <- summary_df %>%
    summarise(
      Mean_Linear_Rate = round(mean(Linear_Rate_mg_DO_per_L_per_H, na.rm = TRUE), 2),
      Mean_Combined_Rate = round(mean(Combined_Initial_Rate_mg_DO_per_L_per_H, na.rm = TRUE), 2),
      Min_Rate = round(min(Best_Rate_mg_DO_per_L_per_H, na.rm = TRUE), 2),
      Max_Rate = round(max(Best_Rate_mg_DO_per_L_per_H, na.rm = TRUE), 2),
      Mean_Best_Rate = round(mean(Best_Rate_mg_DO_per_L_per_H, na.rm = TRUE), 2)
    )
  
  print(rate_stats)
}

message("Analysis complete!")