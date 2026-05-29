# ================================
# FTICR-MS presence/absence NMDS, PERMANOVA, and envfit
# YEP v2 data package
# ================================

library(tidyverse)
library(vegan)
library(grid)

# Run from repo root or from src/.
if (!dir.exists("v2_data") && dir.exists("../v2_data")) {
  setwd("..")
}

set.seed(123)

# ----------------
# paths and settings
# ----------------
out_dir <- file.path("Figures", "fticr_nmds_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treatment_colors <- c(
  "Control" = "#0072B2",
  "Unburned + DOC" = "darkgreen",
  "High Burn + DOC" = "#8B4513"
)

moisture_shapes <- c(
  "Dry" = 16,
  "Wet" = 17
)

permutations <- 9999
envfit_alpha <- 0.05
label_points <- FALSE


fticr_data_path <- file.path(
  "v2_data", "v2_YEP_Sample_Data", "FTICR",
  "YEP_Sediment_CoreMS_Processed_ICR_Data.csv"
)
sample_metadata_path <- file.path("v2_data", "v2_YEP_Sample_Name_Metadata.csv")

npoc_path <- file.path("v2_data", "v2_YEP_Sample_Data", "YEP_Sediment_NPOC_TN.csv")
gas_path <- file.path("v2_data", "v2_YEP_Sample_Data", "v2_YEP_Sediment_CO2_CH4_N2O.csv")
ions_path <- file.path("v2_data", "v2_YEP_Sample_Data", "YEP_Sediment_Ions.csv")
resp_path <- file.path("v2_data", "v2_YEP_Sample_Data", "YEP_Sediment_Incubations_Respiration_Rates.csv")

#  envfit vectors.
env_var_groups <- list(
  npoc_tn = c(
    "Extractable_NPOC_mg_per_kg",
    "Extractable_TN_mg_per_kg"
  ),
  gases = c(
    "Rate_CO2_moles_per_L_per_hr",
    "Rate_CH4_moles_per_L_per_hr",
    "Rate_N2O_moles_per_L_per_hr"
  ),
  ions = c(
    "00940_Cl_mg_per_L",
    "00945_SO4_mg_per_L_as_SO4",
    "00618_NO3_mg_per_L_as_N",
    "00915_Ca_mg_per_L",
    "00925_Mg_mg_per_L",
    "00935_K_mg_per_L",
    "00930_Na_mg_per_L"
  ),
  respiration = c(
    "Normalized_Respiration_Rate_mg_DO_per_H_per_kg_dry_sediment"
  )
)

# ----------------
# Helpers
# ----------------
map_treatment <- function(x) {
  case_when(
    x == "Columbia River synthetic river water" ~ "Control",
    x == "Unburned Douglas fir leachate" ~ "Unburned + DOC",
    x == "Burned (high severity) Douglas fir leachate" ~ "High Burn + DOC",
    TRUE ~ NA_character_
  )
}

map_moisture <- function(x) {
  case_when(
    str_to_lower(str_squish(x)) == "dry" ~ "Dry",
    str_to_lower(str_squish(x)) == "wet" ~ "Wet",
    TRUE ~ NA_character_
  )
}

fmt_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return("< 0.001")
  }
  paste0("= ", sprintf("%.3f", p))
}

to_number <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^<", "")
  x[x %in% c("", "NA", "N/A", "-9999")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

mean_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

read_yep_table <- function(path) {
  read_csv(
    path,
    skip = 2,
    na = c("", "NA", "N/A"),
    show_col_types = FALSE,
    name_repair = "minimal",
    trim_ws = TRUE
  ) %>%
    filter(!is.na(Sample_Name), str_detect(Sample_Name, "^YEP"))
}

make_join_key <- function(treatment, moisture, day, replicate) {
  good_key <- !is.na(treatment) & !is.na(moisture) &
    !is.na(day) & !is.na(replicate)

  if_else(
    good_key,
    paste(treatment, moisture, day, replicate, sep = "|"),
    NA_character_
  )
}

pretty_env_name <- function(x) {
  recode(
    x,
    "Extractable_NPOC_mg_per_kg" = "NPOC",
    "Extractable_TN_mg_per_kg" = "TN",
    "Rate_CO2_moles_per_L_per_hr" = "CO2 rate",
    "Rate_CH4_moles_per_L_per_hr" = "CH4 rate",
    "Rate_N2O_moles_per_L_per_hr" = "N2O rate",
    "00940_Cl_mg_per_L" = "Cl",
    "00945_SO4_mg_per_L_as_SO4" = "SO4",
    "00618_NO3_mg_per_L_as_N" = "NO3",
    "00915_Ca_mg_per_L" = "Ca",
    "00925_Mg_mg_per_L" = "Mg",
    "00935_K_mg_per_L" = "K",
    "00930_Na_mg_per_L" = "Na",
    "Respiration_Rate_mg_DO_per_L_per_H" = "DO respiration",
    "DO_Concentration_At_Incubation_Time_Zero" = "Initial DO",
    .default = x
  )
}

prepare_env_table <- function(path, vars, sample_design) {
  dat <- read_yep_table(path)
  vars_present <- intersect(vars, names(dat))

  if (length(vars_present) == 0) {
    warning("No requested envfit variables found in ", path)
    return(tibble(join_key = character()))
  }

  missing_vars <- setdiff(vars, names(dat))
  if (length(missing_vars) > 0) {
    warning(
      "Missing requested envfit variables in ", basename(path), ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  dat %>%
    select(Sample_Name, all_of(vars_present)) %>%
    mutate(across(all_of(vars_present), to_number)) %>%
    left_join(
      sample_design %>%
        select(Sample_Name, Treatment, Moisture, Incubation_Day, Replicate_Identifier, join_key),
      by = "Sample_Name"
    ) %>%
    filter(!is.na(join_key)) %>%
    group_by(join_key) %>%
    summarise(across(all_of(vars_present), mean_or_na), .groups = "drop")
}

save_adonis <- function(model, path) {
  write_lines(capture.output(print(model)), path)
}

target_sample_material <- "Aqueous sample post-incubation (treatment solution and sediment)"
# ----------------
# Sample metadata
# ----------------
sample_design <- read_csv(
  sample_metadata_path,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    Sample_Material = str_squish(Sample_Material),
    Incubation_Day = str_squish(Incubation_Day),
    Replicate_Identifier = str_squish(Replicate_Identifier),
    Treatment = map_treatment(str_squish(Treatment)),
    Moisture = map_moisture(Sediment_Field_Moisture_Conditions),
    Treatment = factor(Treatment, levels = names(treatment_colors)),
    Moisture = factor(Moisture, levels = names(moisture_shapes)),
    Incubation_Day = factor(Incubation_Day, levels = c("Day_1", "Day_2", "Day_3")),
    join_key = make_join_key(Treatment, Moisture, Incubation_Day, Replicate_Identifier)
  )

# ----------------
# FTICR presence/absence matrix
# ----------------
fticr_raw <- read_csv(
  fticr_data_path,
  show_col_types = FALSE,
  name_repair = "minimal",
  trim_ws = TRUE
)

# the next line on the script is where we decide which samples to include in the NMDS based on the sample metadata and the presence of non-NA values in the FTICR data. We also convert the FTICR data to numeric and create a presence/absence matrix for the NMDS analysis.
sample_cols <- setdiff(names(fticr_raw), "Calibrated_Mass")

fticr_numeric <- fticr_raw %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), to_number))

pa_all <- t(ifelse(as.matrix(fticr_numeric) > 0, 1, 0))
pa_all[is.na(pa_all)] <- 0
storage.mode(pa_all) <- "numeric"
rownames(pa_all) <- sample_cols
colnames(pa_all) <- paste0("mol_", seq_len(nrow(fticr_raw)))

fticr_sample_info <- tibble(Sample_Name = rownames(pa_all)) %>%
  left_join(sample_design, by = "Sample_Name") %>%
  filter(
    Sample_Material == target_sample_material,
    !is.na(Treatment),
    Moisture %in% names(moisture_shapes)
  ) %>%
  arrange(Sample_Name)

pa <- pa_all[fticr_sample_info$Sample_Name, , drop = FALSE]
pa <- pa[, colSums(pa) > 0, drop = FALSE]
pa <- pa[rowSums(pa) > 0, , drop = FALSE]

fticr_sample_info <- fticr_sample_info %>%
  filter(Sample_Name %in% rownames(pa)) %>%
  arrange(match(Sample_Name, rownames(pa)))

stopifnot(identical(fticr_sample_info$Sample_Name, rownames(pa)))

write_csv(
  fticr_sample_info,
  file.path(out_dir, "fticr_nmds_sample_metadata.csv")
)

cat("NMDS sample count:", nrow(pa), "\n")
cat("Presence/absence molecule count:", ncol(pa), "\n")

# ----------------
# NMDS and PERMANOVA
# ----------------
dist_jaccard <- vegdist(pa, method = "jaccard", binary = TRUE)

nmds <- metaMDS(
  pa,
  distance = "jaccard",
  binary = TRUE,
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

site_scores <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(fticr_sample_info, by = "Sample_Name")

set.seed(123)
permanova_main <- adonis2(
  dist_jaccard ~ Treatment + Moisture + Incubation_Day,
  data = fticr_sample_info,
  permutations = permutations,
  by = "margin"
)

set.seed(123)
permanova_interaction <- adonis2(
  dist_jaccard ~ Treatment * Moisture + Incubation_Day,
  data = fticr_sample_info,
  permutations = permutations,
  by = "terms"
)

save_adonis(
  permanova_main,
  file.path(out_dir, "fticr_pa_permanova_main_effects.txt")
)
save_adonis(
  permanova_interaction,
  file.path(out_dir, "fticr_pa_permanova_interaction_model.txt")
)

write_csv(
  as.data.frame(permanova_main) %>% rownames_to_column("Term"),
  file.path(out_dir, "fticr_pa_permanova_main_effects.csv")
)
write_csv(
  as.data.frame(permanova_interaction) %>% rownames_to_column("Term"),
  file.path(out_dir, "fticr_pa_permanova_interaction_model.csv")
)

permanova_tab <- as.data.frame(permanova_main)
label_main <- paste0(
  "Stress = ", round(nmds$stress, 3),
  "\nTreatment: R2 = ", round(permanova_tab["Treatment", "R2"], 3),
  ", p ", fmt_p(permanova_tab["Treatment", "Pr(>F)"]),
  "\nMoisture: R2 = ", round(permanova_tab["Moisture", "R2"], 3),
  ", p ", fmt_p(permanova_tab["Moisture", "Pr(>F)"])
)

base_nmds_plot <- ggplot(
  site_scores,
  aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Moisture)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  scale_shape_manual(values = moisture_shapes, drop = FALSE) +
  labs(
    title = "FTICR-MS presence/absence NMDS",
    subtitle = label_main,
    x = "NMDS1",
    y = "NMDS2",
    color = "Treatment",
    shape = "Sediment field moisture"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

if (label_points) {
  base_nmds_plot <- base_nmds_plot +
    geom_text(aes(label = Sample_Name), size = 2.5, vjust = -0.8, check_overlap = TRUE)
}

print(base_nmds_plot)

ggsave(
  file.path(out_dir, "fticr_pa_nmds_treatment_moisture.png"),
  base_nmds_plot,
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "fticr_pa_nmds_treatment_moisture.pdf"),
  base_nmds_plot,
  width = 8,
  height = 6
)

write_csv(site_scores, file.path(out_dir, "fticr_pa_nmds_site_scores.csv"))

# ----------------
# Envfit variables
# ----------------
env_tables <- list(
  prepare_env_table(npoc_path, env_var_groups$npoc_tn, sample_design),
  prepare_env_table(gas_path, env_var_groups$gases, sample_design),
  prepare_env_table(ions_path, env_var_groups$ions, sample_design),
  prepare_env_table(resp_path, env_var_groups$respiration, sample_design)
)

env_joined <- reduce(env_tables, full_join, by = "join_key")

env_input <- fticr_sample_info %>%
  select(Sample_Name, join_key) %>%
  left_join(env_joined, by = "join_key") %>%
  arrange(match(Sample_Name, rownames(pa)))

env_vars <- setdiff(names(env_input), c("Sample_Name", "join_key"))
env_mat <- env_input %>%
  select(all_of(env_vars)) %>%
  as.data.frame()
rownames(env_mat) <- env_input$Sample_Name

env_vars_keep <- names(env_mat)[
  map_lgl(env_mat, ~ sum(!is.na(.x)) >= 4 && sd(.x, na.rm = TRUE) > 0)
]

if (length(env_vars_keep) == 0) {
  stop("No envfit variables had at least 4 non-missing values and non-zero variance.")
}

env_mat <- env_mat[, env_vars_keep, drop = FALSE]

set.seed(123)
env_fit <- envfit(
  nmds,
  env_mat,
  permutations = permutations,
  na.rm = TRUE
)

write_csv(
  env_input,
  file.path(out_dir, "fticr_pa_envfit_input_variables.csv")
)
write_lines(
  capture.output(print(env_fit)),
  file.path(out_dir, "fticr_pa_envfit_results.txt")
)

env_vectors <- scores(env_fit, display = "vectors") %>%
  as.data.frame() %>%
  rownames_to_column("Variable")
names(env_vectors)[2:3] <- c("NMDS1", "NMDS2")

env_vectors <- env_vectors %>%
  mutate(
    r2 = env_fit$vectors$r,
    p_value = env_fit$vectors$pvals,
    Label = pretty_env_name(Variable)
  ) %>%
  arrange(p_value)

write_csv(env_vectors, file.path(out_dir, "fticr_pa_envfit_vectors.csv"))

env_vectors_plot <- env_vectors %>%
  filter(p_value <= envfit_alpha)

if (nrow(env_vectors_plot) == 0) {
  message(
    "No envfit vectors passed p <= ", envfit_alpha,
    ". Plotting the 8 lowest p-value vectors instead."
  )
  env_vectors_plot <- env_vectors %>% slice_head(n = min(8, n()))
}

arrow_multiplier <- ordiArrowMul(
  as.matrix(env_vectors_plot[, c("NMDS1", "NMDS2")])
)

env_vectors_plot <- env_vectors_plot %>%
  mutate(
    xend = NMDS1 * arrow_multiplier,
    yend = NMDS2 * arrow_multiplier,
    label_x = xend * 1.08,
    label_y = yend * 1.08,
    Label = paste0(Label, "\np=", sprintf("%.3f", p_value))
  )

envfit_plot <- base_nmds_plot +
  geom_segment(
    data = env_vectors_plot,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.5,
    color = "gray20"
  ) +
  geom_text(
    data = env_vectors_plot,
    aes(x = label_x, y = label_y, label = Label),
    inherit.aes = FALSE,
    size = 3,
    color = "gray20"
  ) +
  labs(
    title = "FTICR-MS presence/absence NMDS with envfit vectors",
    caption = paste0(
      "Envfit vectors shown at p <= ", envfit_alpha,
      "; if none pass, the 8 lowest p-value vectors are shown."
    )
  )

print(envfit_plot)

ggsave(
  file.path(out_dir, "fticr_pa_nmds_envfit.png"),
  envfit_plot,
  width = 8.5,
  height = 6.5,
  dpi = 300
)
ggsave(
  file.path(out_dir, "fticr_pa_nmds_envfit.pdf"),
  envfit_plot,
  width = 8.5,
  height = 6.5
)

cat("Done. Outputs written to:", normalizePath(out_dir, winslash = "/"), "\n")
