# YEP Incubations

This repo is organized around the current ST2A immediate impacts manuscript.
The active analysis inputs are in `v2_data/`; older downloaded data packages
and exploratory analyses are archived under `Archive_2/`.

## Active workflow

Run scripts from the repo root, in this order:

1. `src/Model_DO_biotic_abiotic_v2.R`
   - Generates model-fit outputs in `modeling_outputs/`.
   - Main downstream table: `modeling_outputs/YEP_Complete_Analysis_Following_Paper_Enhanced.csv`.
   - SI model comparison: `modeling_outputs/Model_Performance_Summary.pdf` and `modeling_outputs/YEP_Model_Performance_Table.csv`.

2. `src/Table_1.R`
   - Generates Table 1 summary CSVs in `Tables/`.
   - Uses v2 respiration, CO2, NPOC, and ion data.

3. `src/Plot_fitted_Vmax_kL_DryWet_Wilcox.R`
   - Generates Figure 1 outputs in `Figures/`.

4. `src/Plot_fitted_Vmax_kL_Treatment_KW_Dunn.R`
   - Generates Figure 2 outputs in `Figures/`.

5. `src/FTICR_NMDS_presence_absence_envfit.R`
   - Generates Figure 3 NMDS/PERMANOVA outputs in `Figures/fticr_nmds_v2/`.

6. `src/Plot_nitrate_v2.R`
   - Generates nitrate SI/Table 1 support outputs in `Figures/nitrate_v2/`.

## Archive

`Archive_2/manuscript_cleanup_2026-05-29/` contains scripts and outputs from
exploratory analyses that are not linked to the current manuscript figures,
Table 1, or SI placeholders.
