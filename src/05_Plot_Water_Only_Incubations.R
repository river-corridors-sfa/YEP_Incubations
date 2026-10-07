# === 05 Plot water-only (no sediment) incubations =====
# Created by: VGC
# Plot DO vs time for the water-only control incubations, fit a simple
#          linear regression to each one, and show the regression stats inside
#          each panel. All samples are shown in one figure as separate panels.
#
rm(list = ls(all = TRUE))
library(dplyr)
library(ggplot2)

# ==== User input =====
do_file        <- "v2_data/Merged_DO_Firesting_Data_2026-10-07.csv"
water_pattern  <- "_INC-W"          # identifies water-only samples
trim_start_min <- 2                 # minutes removed from the start (set to 0 to keep all)
trim_end_min   <- 2                 # minutes removed from the end (set to 0 to keep all)
y_limits       <- c(0, 10)          # fixed y-axis range (mg/L) for all panels
n_columns      <- 3                 # number of panel columns
output_prefix  <- paste0("YEP_Water_Only_Incubations_", Sys.Date())

# ==== Read data =====
# The published data-package file has 2 extra lines above the column names
first_line <- readLines(do_file, n = 1)
skip_rows  <- if (startsWith(first_line, "#Columns")) 2 else 0

water <- read.csv(do_file, skip = skip_rows, check.names = FALSE) %>%
  filter(grepl(water_pattern, Sample_Name)) %>%          # also drops header/unit rows
  mutate(Elapsed_Seconds = as.numeric(Elapsed_Seconds),
         DO_mg_per_L     = as.numeric(DO_mg_per_L)) %>%
  filter(!is.na(DO_mg_per_L), DO_mg_per_L != -9999) %>%
  group_by(Sample_Name) %>%
  arrange(Elapsed_Seconds, .by_group = TRUE) %>%
  mutate(Time_h  = (Elapsed_Seconds - min(Elapsed_Seconds)) / 3600,
         Used_In_Regression = Time_h >  trim_start_min / 60 &
                              Time_h < (max(Time_h) - trim_end_min / 60)) %>%
  ungroup()

# ==== Linear regression for each sample =====
regression_results <- water %>%
  filter(Used_In_Regression) %>%
  group_by(Sample_Name) %>%
  group_modify(~ {
    fit <- summary(lm(DO_mg_per_L ~ Time_h, data = .x))
    data.frame(
      Slope_mg_DO_per_L_per_H = fit$coefficients[2, 1],
      Slope_SE                = fit$coefficients[2, 2],
      Intercept               = fit$coefficients[1, 1],
      R_Squared               = fit$r.squared,
      R_Squared_Adj           = fit$adj.r.squared,
      p_value                 = fit$coefficients[2, 4],
      Number_Points           = nrow(.x),
      Regression_Time_h       = max(.x$Time_h) - min(.x$Time_h),
      Mean_DO_mg_per_L        = mean(.x$DO_mg_per_L)
    )
  }) %>%
  ungroup()

# Text shown inside each panel
labels <- regression_results %>%
  mutate(label = paste0(
    "Slope: ", formatC(Slope_mg_DO_per_L_per_H, format = "f", digits = 3), " mg/L/h\n",
    "R\u00b2 = ", formatC(R_Squared, format = "f", digits = 2), "\n",
    "p ", ifelse(p_value < 0.001, "< 0.001", paste("=", formatC(p_value, format = "f", digits = 3))), "\n",
    "n = ", Number_Points
  ))

# ==== Plot =====
p <- ggplot(water, aes(x = Time_h, y = DO_mg_per_L)) +
  geom_point(data = filter(water, !Used_In_Regression),
             colour = "grey70", size = 0.6) +
  geom_point(data = filter(water, Used_In_Regression),
             colour = "steelblue", size = 0.6, alpha = 0.6) +
  geom_smooth(data = filter(water, Used_In_Regression),
              method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.7) +
  geom_label(data = labels, aes(label = label),
             x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.1,
             size = 3, label.size = 0, fill = "white", alpha = 0.85,
             inherit.aes = FALSE) +
  facet_wrap(~ Sample_Name, ncol = n_columns) +
  scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = 2)) +
  labs(x = "Incubation time (h)", y = "DO (mg/L)",
       ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        plot.subtitle = element_text(size = 9))

# ==== Export =====
n_rows <- ceiling(length(unique(water$Sample_Name)) / n_columns)
ggsave(paste0(output_prefix, ".pdf"), p, width = 3.5 * n_columns, height = 3 * n_rows + 1)
ggsave(paste0(output_prefix, ".png"), p, width = 3.5 * n_columns, height = 3 * n_rows + 1, dpi = 300)
write.csv(regression_results, paste0(output_prefix, "_Regression_Results.csv"), row.names = FALSE)

