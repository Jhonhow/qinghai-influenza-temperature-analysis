# Candidate analysis-code bundle (review before release)

This is a curated, **not-yet-publishable** copy of scripts selected as most relevant to the current GeoHealth manuscript. Source files in the project archive remain unchanged. Candidate copies have limited path/source-call portability edits; these edits were not tested by refitting models.

## Suggested workflow and dependencies

Run from this directory (the candidate root), after installing the packages used by the scripts and arranging authorized inputs:

1. `R/00_config.R` — shared configuration; then `01_data_preparation/01_audit_raw_data.R` through `03_freeze_analysis_protocol_inputs.R` for input checks/protocol records.
2. `01_data_preparation/07_extract_daily_grid_to_admin_mean.R` — reconstructed, configurable grid-to-polygon daily-mean template. It has not been run against the original ChinaMet files and is not the recovered historical script. Then `05_merge_daily_weather_inputs.py` and `06_build_daily_weather_workbook.py` merge/build county-day weather inputs.
3. `01_data_preparation/04_build_county_week_panel_v4.R` — build the weekly panel from county-day inputs and the prior panel. This is the downstream aggregation boundary represented in the available code.
4. `02_primary_models/80_primary_models.R` — primary models; `81`–`83` are confirmatory/context analyses. `61` and `75` contain subtype/lag analyses.
5. `03_robustness_validation/92`–`96` — robustness, temporal sensitivity, and supplementary result tables.
6. `04_figure_generation/` and `05_table_generation/` — figure and table outputs. Confirm each script against the final manuscript/figure files before release.

Some robustness scripts source other scripts, so run from the candidate root and retain the folder structure. Inputs, output folders, and permissions may still require author-specific configuration.

## Paths and inputs

- `QH_PROJECT_ROOT`: project/data root. Defaults to the current working directory where used.
- `QH_MET_WORKBOOK`: optional county-day meteorology workbook; default is `<QH_PROJECT_ROOT>/data/updated_met_t2m_pressure.xlsx`.
- `QH_COUNTY_SHAPEFILE`: optional county boundary shapefile; default is `data/qinghai_shp/china_county.shp` relative to the working directory.
- For script 07, set `QH_GRID_MANIFEST`, `QH_BOUNDARY_SHP`, and `QH_ADMIN_ID_FIELD`; optionally set `QH_ADMIN_NAME_FIELD`, `QH_ADMIN_LEVEL`, `QH_BOUNDARY_FILTER_FIELD`, and `QH_BOUNDARY_FILTER_VALUE`. The bundled county shapefile contains 2,891 national county features, so filter to Qinghai or use a Qinghai-only boundary file before extraction.
- Surveillance count data and the county-week analysis panel are restricted and are not included in this code bundle. Do not add them to a public repository.
- The analyst must supply the corresponding ChinaMet county-day derivatives and other covariates. Original gridded files may be large and redistribution rights are unconfirmed.

## Scope and known gaps

The author reports that daily 0.1° ChinaMet grid values were averaged within the 2024 Tianditu provincial, prefectural, and county administrative polygons (approval no. GS(2024)0650) to produce daily summaries, then county daily data were aggregated to epidemiological weeks. The author has confirmed that cells are included when their centers fall inside a county polygon. The historical extraction code/log and exact NoData handling remain unavailable. Script `07_extract_daily_grid_to_admin_mean.R` is a transparent reconstruction of the confirmed unweighted cell-center rule, not the recovered historical script; no ChinaMet raster was available for research-data validation. See `docs/DATA_PROVENANCE_NOTE.md`. No monthly meteorological exposure workflow was located.

The author confirms that the analysis used R 4.5.1. The script writes a sessionInfo() record for the environment in which the template is run; this is not evidence of the historical package versions.

The new script passed an R 4.5.1 syntax check and a small synthetic two-polygon smoke test. No ChinaMet raster was available in this package, so no research extraction or output was produced.

The local `.mjs` table-formatting helpers are deliberately not included because they depend on a machine-local Codex `@oai/artifact-tool` runtime; use a portable, documented workflow if table formatting is part of the release. This candidate has not been run end-to-end, dependencies/environment locks and a license are missing, and some script names/lineage (notably Result 3) need author confirmation. **Do not make this repository public until those items and the authors' sharing approval are addressed.**
