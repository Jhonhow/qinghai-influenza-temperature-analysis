# Qinghai influenza–temperature analysis

This repository contains R and supporting scripts for our ecological time-series study of influenza activity and ambient temperature across counties in Qinghai, China. The analysis covers data preparation, primary models, robustness analyses, and the generation of figures and tables for the manuscript.

The code is provided to document the analytical workflow. It does not include restricted surveillance records or the full set of meteorological input files. Please read the data and reproducibility notes below before attempting to run the scripts.

## Analysis workflow

The scripts are organized by stage. Run them from the repository root and retain the folder structure so that relative paths and sourced scripts resolve correctly.

1. **Configure and check inputs:** `R/00_config.R`, followed by the relevant scripts in `01_data_preparation/` (including the input audit and protocol records).
2. **Prepare county-level weather inputs:** use the weather preparation scripts in `01_data_preparation/` to create or combine county-day inputs, then run `04_build_county_week_panel_v4.R` to construct the county-week analysis panel.
3. **Fit the models:** run `02_primary_models/80_primary_models.R` for the primary analyses. Other scripts in this folder contain confirmatory, contextual, subtype, and lag analyses.
4. **Assess robustness:** scripts `92`–`96` in `03_robustness_validation/` cover the robustness and temporal-sensitivity analyses and related supplementary outputs.
5. **Generate outputs:** scripts in `04_figure_generation/` and `05_table_generation/` generate figures and tables.

The scripts may require local input and output paths to be configured before use. The workflow above describes the intended order based on the available scripts; it is not a claim that the repository has been tested as a complete, one-click pipeline in a clean environment.

## Data and access

- **Influenza surveillance data:** county-week aggregated surveillance records were used in the study. These records are not included because they are subject to data-access restrictions. Individual-level records are also not included. For information about authorized access, contact the corresponding author, Fenglu Zhao (`fengluzhao@cams.cn`) and Xuli Li (`qhcdcxull@126.com`). Access is subject to the relevant institutional approvals.
- **Meteorological data:** the study used daily 0.1° ChinaMet gridded meteorological data. The original gridded files are not included in this repository. The author-provided processing description is to average grid-cell values within administrative polygons, including cells whose centers fall inside each polygon, and then aggregate county-day values to epidemiological weeks.
- **Administrative boundaries:** the described processing used 2024 provincial, prefectural, and county boundaries from the National Platform for Common Geospatial Information Services (Tianditu), approval number GS(2024)0650. Boundary files are not bundled here; users should obtain them from the provider and check the applicable terms before use.

### Weather extraction script

`01_data_preparation/07_extract_daily_grid_to_admin_mean.R` implements the confirmed cell-center inclusion rule as a configurable grid-to-polygon daily-mean template. It was added to make the described processing approach reproducible; it is not the historical extraction script used to produce the study data. The original ChinaMet rasters and historical extraction log were not available for validation, so this template has not been tested against the study rasters. It passed a syntax check and a small synthetic test. Results from this template should therefore be checked against the study's authorized inputs before being used for research outputs.

## Software

The study analysis was conducted using R 4.5.1. Package versions from the original analysis environment were not available. When running the scripts, record the R and package versions in the active environment (for example, with `sessionInfo()`).

## Reuse

Please check the repository for a `LICENSE` file before reusing or redistributing code. In the absence of a license, this repository does not specify reuse permissions. Data access remains subject to the restrictions described above.
