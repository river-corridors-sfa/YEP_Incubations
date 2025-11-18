# This script has been modified from Maggi Laan's scripts. It applies the following rules to determine DO consuption rates
# Rules
# * Minimum Number of Points (min.points) = 2
# + Keep at least 2 points to generate a slope
# 
# * Slope Threshold (slope.thresh) = -0.006
# + Don't remove points from slopes if the initial fit has a slope greater than this value. This value was chosen manually. 
#   
# * Low Dissolved Oxygen Threshold (do.thresh) = 2
#   + Used in initial removal of data points below this value in conjunction with time threshold
#   
# * Time Threshold (time.thresh) = 4
#   + Used in conjunction with Dissolved Oxygen Threshold to remove points below Dissolved Oxygen Threshold at a certain time
#   
# * High Dissolved Oxygen Threshold (high.do) = 14
#   + Initial removal of data points above this value
# 
# * Fast Rate Threshold (fast) = 5.5 
#   + Use theoretical saturated 0 minute value if 2 minute value is below this threshold
#    + 5.5 has least number of differences between theoretical and real slopes between replicates
#   + Tested range: 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, No removal
#  
# * Breusch-Pagan Heteroscedasticity Test (bp.fit) = 0.1
#   + p-value used to determine normal distribution of residuals
#   + 0.1 middle of selection, no large differences in histograms
#   + Tested range: 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, No removals
# 
# * High Slope Threshold (high.slope.thresh) = -0.04
#   + Remove saturation point if slope is pretty flat     but starts low (might be able to remove)
# 
# * Time Same (time.same) = 7
#   + Keep 2 minute DO value even if it was put on the rollers at the same time a photo was taken if it is lower than this value
#   
# * Break Point Toggle (break.toggle) = 1.4 
#   + Ratio between slope before/after breaking data at predicted segmented regression point
#   + 1.4 was best for breaking more non-linear slopes upon manual inspection, however changing the ratio doesn't change histograms
# + Tested range: 1.0, 1.2, 1.4, 1.6, 1.8, 2.0
# 
# * Concentration Toggle = 1.4
# + Range between first and last DO concentration in data frame, used to keep more linear slopes from breaking/having points removed in heteroscedasticity 
# + 1.4 was best for not removing data//using segmented regression for more linear samples
# + Tested range: 0, 0.5, 1.0, 1.4, 1.5, 2.0

rm(list=ls())
# ===== Load libraries =====
library(dplyr); library(ggplot2);library(ggsignif)
library(ggpubr);library(reshape2);library(ggpmisc)
library(segmented);library(broom);library(lmtest);library(car); library(corrplot)
library(ggpmisc);library(lubridate); library(readxl); library(glmnet)
library(tidyverse);library(patchwork)
library(readr); library(pdftools)

# ======== User set parameters based on the rules =====

min.points = 2
slope.thresh = -0.006
do.thresh = 2
high.do = 14
time.thresh.sec = 4 * 60 # in seconds
fast.time.sec = 2 * 60              # 2 minutes = 120 seconds
fast = 5.5 
bpfit = 0.1
high.slope.thresh = -0.04 
time.same = 7
break.toggle = 1.4
conc.range = 1.4
plot.toggle = TRUE

# ===== Read in data =====

data = read.csv('YEP_DP_downloaded_11-18-25/YEP_Sample_Data/YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temperature.csv', skip = 2) %>%
  filter(grepl("YEP", Sample_Name)) %>%
  dplyr::select(-c(Field_Name, IGSN, Material, Methods_Deviation))

samples = data %>% 
  filter(DO_mg_per_L != -9999)


# ====== Data pre-cleaning steps  =======
# Make Data Frames
respiration <- as.data.frame(matrix(NA, ncol = 17, nrow =1))

colnames(respiration) = c("Sample_Name","slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj",  "p_value", "total_incubation_time_min", "number_of_points", "no_points_rem",  "breusch_p_value","break_point", "slope_ratio",  "0_min_concentration", "2_min_concentration", "last_concentration", "theoretical")

location = c("-H1", "-H2", "-H3","-H4", "-H5","-U1", "-U2", "-U3", "-U4", "-U5","-S1", "-S2", "-S3", "-S4", "-S5")

for (i in 1:length(location)){
  
  ## Subset by replicate
  
  data_location_subset = samples[grep(location[i],samples$Sample_Name),]
  
  unique.incubations = unique(data_location_subset$Sample_Name)
  
  rate = as.data.frame(matrix(NA, ncol = 14, nrow = length(unique(samples$Sample_Name))))
  
  colnames(rate) = c("Sample_Name", "slope_of_the_regression", "rate_mg_per_L_per_min", "rate_mg_per_L_per_h", "R_squared", "R_squared_adj", "p_value", "total_incubation_time_min", "breusch_p_value", "slope_ratio", "0_min_concentration", "2_min_concentration", "last_concentration", "theoretical")
  
  for (j in 1:length(unique.incubations)){
    
    # Subset by individual sites
    
    data_site_subset = subset(data_location_subset, data_location_subset$Sample_Name == unique.incubations[j])
    
    # Put in ascending order for minutes
    data_site_subset = data_site_subset[order(data_site_subset$Elapsed_Seconds, decreasing = FALSE),]
    data_site_subset$Elapsed_Seconds = as.numeric(data_site_subset$Elapsed_Seconds)
    
    # Fit linear model to all data
    
    fit_all = lm(DO_mg_per_L~Elapsed_Seconds, data = data_site_subset)
    
    # Remove points greater than High Dissolved Oxygen Threshold
    
    data_site_subset_low = data_site_subset %>% 
      filter(DO_mg_per_L < high.do) 
    
    # Fit linear model to slopes with High Dissolved Oxygen Points Removed
    
    fit_low = lm(DO_mg_per_L~Elapsed_Seconds, data = data_site_subset_low)
    
    low.slope = fit_low$coefficients[[2]]
   
     # Calculate the start time for this specific sample
    start_time = min(data_site_subset_low$Elapsed_Seconds)
    
    # Remove samples if at > 4 minutes FROM START, they are below the DO threshold
    data_site_subset_thresh = data_site_subset_low %>%
      filter(!((Elapsed_Seconds - start_time) > time.thresh.sec & DO_mg_per_L < do.thresh))
    
    data_site_subset_rem = data_site_subset_thresh
    
    # If the 2 minute concentration (FROM START) is less than the DO threshold, 
    # remove anything greater than 2 minutes FROM START
    two_min_target_time = start_time + fast.time.sec  # This is the actual "2 minute" timepoint
    two_min_data = data_site_subset_rem %>% 
      filter(abs(Elapsed_Seconds - two_min_target_time) == min(abs(Elapsed_Seconds - two_min_target_time)))
    
    if(nrow(two_min_data) > 0 && two_min_data$DO_mg_per_L[1] <= do.thresh){
      data_site_subset_rem = data_site_subset_rem %>%
        filter((Elapsed_Seconds - start_time) <= fast.time.sec)
    }
    
    data_site_subset_beg =  data_site_subset_low # for now as I figure out some stuff

# ====== Fit linear models ========
    
    # Fit linear model to cleaned data (low points removed, high points removed, first point removed if put on rollers at same time as picture) 
    
    fitog = lm(DO_mg_per_L~Elapsed_Seconds, data = data_site_subset_beg)
    slope_original = fitog$coefficients[[2]] #slope of original cleaned data
    
    # Calculate segmented regression slopes:
    
    if (nrow(data_site_subset_beg) > 3) {
      
      ## Only works if more than 3 rows of data, fits segmented regression and gets best estimate of where to break the data
      
      # psi gives estimated start point for the analysis
      
      midpoint <- (max(data_site_subset_beg$Elapsed_Seconds) + min(data_site_subset_beg$Elapsed_Seconds)) / 2
      
      segmentog = segmented(fitog, 
                            seg.Z = ~Elapsed_Seconds, 
                            psi = list(Elapsed_Seconds = midpoint))      
      
      fit_seg = numeric(length(data_site_subset_beg$Elapsed_Seconds)) * NA
      
      # gives DO estimates for segmented regression lm
      fit_seg[complete.cases(rowSums(cbind(data_site_subset_beg$DO_mg_per_L, data_site_subset_beg$Elapsed_Seconds)))] <- broken.line(segmentog)$fit
      
      data_seg = data.frame(DO_mg_per_L = data_site_subset_beg$DO_mg_per_L, Elapsed_Seconds = data_site_subset_beg$Elapsed_Seconds, fit = fit_seg)
      
      # 1: Initial Break Point, 2: Chosen Break Point, 3: Approx. Standard Error
      seg = segmentog$psi
      #Chosen Break Point
      est = seg[[2]]
      #Standard Error
      bp_se = seg[[3]]
      #Rounded Break Point
      data_site_subset_beg$break_point = 2*round(seg[[2]]/2)
      
      #Break Data at rounded value
      data_site_subset_break = subset(data_site_subset_beg, Elapsed_Seconds <= data_site_subset_beg$break_point[1])
      
    } else {
      
      # If number of rows < 3, don't fit segmented regression
      
      data_site_subset_beg = data_site_subset_beg
      
      est = "NA"
      
      data_site_subset_break = data_site_subset_beg
      
    }
    
    # Fit linear model to segmented data
    fit_break = lm(DO_mg_per_L~Elapsed_Seconds, data = data_site_subset_break)
    u = fit_break$coefficients
    b = u[[1]] #intercept
    c = u[[2]] #slope
    slope_beginning = c
    r = summary(fit_break)$r.squared #r squared of segmented data
    residuals = deviance(fit_break)
    r.adj = summary(fit_break)$adj.r.squared
    p = summary(fit_break)$coefficients[4]
    pog = p
    rog = r
    resog = residuals
    bp = bptest(fit_break)[[4]]
    bpog = bptest(fit_break)[[4]]
    
    #Calculate ratio of slopes before and after segmented regression
    
    slope.ratio = fit_break$coefficients[[2]]/fitog$coefficients[[2]]
    
    #If ratio of slopes before/after segmented regression, the original slope is not flat, and the range of the concentrations from the unbroken data is > 1.4, then use the segmented regression
    
    if(slope.ratio > break.toggle & slope_original < slope.thresh & (first(data_site_subset_beg$DO_mg_per_L) - last(data_site_subset_beg$DO_mg_per_L)) > conc.range) {
      
      data_site_subset_fin = data_site_subset_break
      
    } else {
      
      data_site_subset_break = data_site_subset_beg
      
      data_site_subset_fin = data_site_subset_break
      
    }
    
    # Fit final lm to data after deciding to use segmented regression or not
    fit_break = lm(DO_mg_per_L~Elapsed_Seconds, data = data_site_subset_fin)
    u = fit_break$coefficients
    b = u[[1]] #intercept
    c = u[[2]] #slope
    slope_beginning = c
    r = summary(fit_break)$r.squared
    residuals = deviance(fit_break)
    r.adj = summary(fit_break)$adj.r.squared
    p = summary(fit_break)$coefficients[4]
    pog = p
    rog = r
    resog = residuals
    bp = bptest(fit_break)[[4]]
    bpog = bptest(fit_break)[[4]]
    
    # Fit lm to data, this will change with loop
    fit = lm(data_site_subset_fin$DO_mg_per_L~data_site_subset_fin$Elapsed_Seconds)
    u = fit$coefficients
    b = u[[1]] #Intercept
    c = u[[2]] #rate mg/L min
    r = summary(fit)$r.squared
    r.adj = summary(fit)$adj.r.squared
    residuals = deviance(fit)
    p = summary(fit)$coefficients[4]
    r2 = r
    res2 = residuals
    bp = bptest(fit)[[4]]
    
    
    
    if (slope_beginning >= slope.thresh | ((first(data_site_subset_beg$DO_mg_per_L) - last(data_site_subset_beg$DO_mg_per_L)) < conc.range)) 
      
    {
      
      #if it has a flat slope, or original data has low concentration range in cleaned data, don't do anything else
      
      data_site_subset_fin = data_site_subset_fin
      
      fit = lm(data_site_subset_fin$DO_mg_per_L~data_site_subset_fin$Elapsed_Seconds)
      u = fit$coefficients
      b = u[[1]] #Intercept
      c = u[[2]] #rate mg/L min
      r = summary(fit)$r.squared
      r.adj = summary(fit)$adj.r.squared
      residuals = deviance(fit)
      p = summary(fit)$coefficients[4]
      r2 = r
      res2 = residuals
      bp = bptest(fit)[[4]]
      
    }
    
    else {
      
      #else start looping to remove data
      
      for (l in 1:60) {
        
        if (nrow(data_site_subset_fin)<= min.points){
          
          #if there are more 2 or less points fit final lm to final data:
          
          data_site_subset_fin = data_site_subset_fin
          
          fit = lm(data_site_subset_fin$DO_mg_per_L~data_site_subset_fin$Elapsed_Seconds)
          u = fit$coefficients
          b = u[[1]] #Intercept
          c = u[[2]] #rate mg/L min
          r = summary(fit)$r.squared
          r.adj = summary(fit)$adj.r.squared
          residuals = deviance(fit)
          p = summary(fit)$coefficients[4]
          r2 = r
          res2 = residuals
          bp = bptest(fit)[[4]]
          
        }
        
        else if (bp < bpfit & nrow(data_site_subset_fin) >= min.points){
          
          #Else if the p-value of the Breusch-Pagan Heteroscedasticity test is less than the threshold, start removing data from the back end until it is greater than threshold or there are only 2 points
          
          data_site_subset_fin = data_site_subset_fin[-nrow(data_site_subset_fin),]
          
          fit = lm(data_site_subset_fin$DO_mg_per_L~data_site_subset_fin$Elapsed_Seconds)
          u = fit$coefficients
          b = u[[1]] #Intercept
          c = u[[2]] #rate mg/L min
          r = summary(fit)$r.squared
          r.adj = summary(fit)$adj.r.squared
          residuals = deviance(fit)
          p = summary(fit)$coefficients[4]
          r2 = r
          res2 = residuals
          bp = bptest(fit)[[4]]
          
        }
        
      }
      
      if (slope_beginning < 0 & c > 0){
        
        # If the slope of the final is positive, use the slope of the cleaned data
        
        data_site_subset_fin = data_site_subset_beg
        
        fit = lm(data_site_subset_fin$DO_mg_per_L~data_site_subset_fin$Elapsed_Seconds)
        u = fit$coefficients
        b = u[[1]] #Intercept
        c = u[[2]] #rate mg/L min
        r = summary(fit)$r.squared
        r.adj = summary(fit)$adj.r.squared
        residuals = deviance(fit)
        p = summary(fit)$coefficients[4]
        r2 = r
        res2 = residuals
        bp = bptest(fit)[[4]]
        
      }
      
    }
    
    my.format <- "Slope: %s\nR2: %s\np: %s"  
    my.formula <- y ~ x
