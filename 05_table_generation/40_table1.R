local_lib <- normalizePath("R_libs_q1", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

source(file.path("R", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(scales)
  library(grid)
  library(gridExtra)
  library(openxlsx)
  library(svglite)
  library(ragg)
})

version_id <- "v3_20260728"
figure_dir <- file.path(project_root, "result0715", "main_text_figures")
table_dir <- file.path(project_root, "result0715", "main_text_tables")
source_dir <- file.path(
  project_root,
  "result0715",
  "source_data",
  "figure1_table1_final_20260805"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

panel_path <- file.path(
  project_root,
  "model_outputs",
  paste0("reanalysis_", version_id),
  "county_week_panel_v3_2016_2025.csv"
)
bundle_path <- file.path(
  project_root,
  "model_outputs",
  paste0("reanalysis_", version_id),
  "regularized_altitude_interaction_20260730",
  "regularized_analysis_bundle.rds"
)

if (!file.exists(panel_path)) stop("Missing v3 county-week panel.")
if (!file.exists(bundle_path)) stop("Missing frozen primary-model bundle.")

panel <- read.csv(
  panel_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) |>
  mutate(week_start = as.Date(week_start))

bundle <- readRDS(bundle_path)
primary <- bundle$full_scaled$train |>
  mutate(week_start = as.Date(week_start))

primary_counties <- sort(unique(primary$county))
if (length(primary_counties) != 45L) {
  stop("Frozen primary model does not contain exactly 45 counties.")
}

study_panel <- panel |>
  filter(
    county %in% primary_counties,
    iso_year >= 2017,
    iso_year <= 2024
  )

met_panel <- study_panel |>
  filter(
    met_valid_5day_main,
    !is.na(temperature),
    !is.na(altitude)
  )

lab_observed <- met_panel |>
  filter(
    lab_available,
    !is.na(tests),
    tests >= 1,
    !is.na(positive),
    positive >= 0,
    positive <= tests
  )

ili_candidates <- study_panel |>
  filter(
    ili_available,
    !is.na(ILI),
    !is.na(visits),
    ILI >= 0,
    visits > 0
  )

invalid_ili <- ili_candidates |>
  filter(ILI > visits) |>
  select(county, city, iso_year, iso_week, week_start, ILI, visits, ili_days, ili_rows, ili_valid_rows)
write_csv_utf8(
  invalid_ili,
  file.path(source_dir, "excluded_logically_invalid_ili_county_weeks.csv")
)

ili_observed <- ili_candidates |>
  filter(ILI <= visits)

if (nrow(primary) != 9758L || sum(primary$tests) != 77930L ||
    sum(primary$positive) != 11072L) {
  stop("Frozen primary-model denominator no longer matches the results manifest.")
}

cold_threshold <- -15

altitude_summary <- primary |>
  distinct(county, city, altitude)

met_county <- met_panel |>
  group_by(county) |>
  summarise(
    met_valid_weeks = n(),
    cold_weeks = sum(temperature <= cold_threshold),
    cold_week_pct = 100 * cold_weeks / met_valid_weeks,
    mean_temperature = mean(temperature),
    .groups = "drop"
  )

lab_county <- lab_observed |>
  group_by(county) |>
  summarise(
    lab_observed_weeks = n(),
    total_tests = sum(tests),
    total_positive = sum(positive),
    positivity = total_positive / total_tests,
    tests_per_lab_week = total_tests / lab_observed_weeks,
    .groups = "drop"
  )

county_context <- altitude_summary |>
  left_join(met_county, by = "county") |>
  left_join(lab_county, by = "county") |>
  mutate(
    lab_observed_weeks = coalesce(lab_observed_weeks, 0L),
    total_tests = coalesce(total_tests, 0),
    total_positive = coalesce(total_positive, 0),
    lab_coverage_pct = 100 * lab_observed_weeks / met_valid_weeks,
    tests_per_lab_week = ifelse(lab_observed_weeks > 0, tests_per_lab_week, NA_real_)
  ) |>
  arrange(altitude)

if (any(is.na(county_context$met_valid_weeks)) ||
    any(is.na(county_context$cold_week_pct))) {
  stop("Meteorological coverage is missing for a primary-analysis county.")
}

write_csv_utf8(
  county_context,
  file.path(source_dir, "figure1_county_context.csv")
)

# Administrative geometry -------------------------------------------------

shp_files <- list.files(data_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
county_shp_path <- shp_files[basename(shp_files) == "中国_县.shp"][1]
if (is.na(county_shp_path)) {
  county_shp_path <- shp_files[!grepl("境界线", basename(shp_files))][1]
}
if (is.na(county_shp_path) || !file.exists(county_shp_path)) {
  stop("Cannot find the Qinghai county polygon shapefile.")
}

qinghai <- st_read(
  county_shp_path,
  quiet = TRUE,
  options = "ENCODING=UTF-8"
) |>
  filter(.data[["省名"]] == "青海省") |>
  mutate(map_county = normalize_county(.data[["县名"]])) |>
  st_make_valid() |>
  group_by(map_county) |>
  summarise(.groups = "drop")

unmatched_primary <- anti_join(
  county_context |> select(county),
  st_drop_geometry(qinghai) |> select(map_county),
  by = c("county" = "map_county")
)
if (nrow(unmatched_primary) > 0) {
  write_csv_utf8(unmatched_primary, file.path(source_dir, "unmatched_primary_counties.csv"))
  stop("At least one primary-analysis county could not be matched to the shapefile.")
}

map_data <- qinghai |>
  left_join(county_context, by = c("map_county" = "county")) |>
  mutate(analytic_county = !is.na(altitude)) |>
  st_transform(3857)

map_join_audit <- st_drop_geometry(map_data) |>
  transmute(
    shapefile_county = map_county,
    analytic_county,
    altitude,
    met_valid_weeks,
    lab_observed_weeks
  ) |>
  arrange(desc(analytic_county), shapefile_county)
write_csv_utf8(map_join_audit, file.path(source_dir, "figure1_map_join_audit.csv"))

# Weekly surveillance ------------------------------------------------------

study_start <- as.Date("2017-01-02")
study_end <- as.Date("2024-12-23")
week_grid <- tibble(week_start = seq(study_start, study_end, by = "7 days"))

lab_weekly <- lab_observed |>
  group_by(week_start) |>
  summarise(
    laboratory_tests = sum(tests),
    laboratory_positive = sum(positive),
    laboratory_positivity_pct = 100 * laboratory_positive / laboratory_tests,
    laboratory_counties = n_distinct(county),
    .groups = "drop"
  )

ili_weekly <- ili_observed |>
  group_by(week_start) |>
  summarise(
    ILI = sum(ILI),
    visits = sum(visits),
    ILI_rate_pct = 100 * ILI / visits,
    ILI_counties = n_distinct(county),
    .groups = "drop"
  )

weekly_surveillance <- week_grid |>
  left_join(lab_weekly, by = "week_start") |>
  left_join(ili_weekly, by = "week_start")
write_csv_utf8(
  weekly_surveillance,
  file.path(source_dir, "figure1_weekly_surveillance.csv")
)

# Table 1 -----------------------------------------------------------------

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_num <- function(x, digits = 1) format(round(x, digits), nsmall = digits, trim = TRUE)
fmt_pct <- function(x, digits = 1) paste0(fmt_num(x, digits), "%")
fmt_date_en <- function(x) {
  paste(as.integer(format(x, "%d")), month.abb[as.integer(format(x, "%m"))], format(x, "%Y"))
}
fmt_median_iqr <- function(x, digits = 1) {
  q <- quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE)
  paste0(fmt_num(q[[2]], digits), " (", fmt_num(q[[1]], digits), "-", fmt_num(q[[3]], digits), ")")
}

altitude_values <- altitude_summary$altitude
ili_counties <- sort(unique(ili_observed$county))
ili_met_panel <- met_panel |> filter(county %in% ili_counties)

table1 <- bind_rows(
  tibble(
    section = "Study setting and meteorological exposure",
    characteristic = c(
      "Study period",
      "Primary-analysis counties, n",
      "County-government-seat altitude, median (IQR); range, m",
      "Meteorologically valid county-weeks, n",
      "Cold weeks at or below -15 °C, county median (IQR), %"
    ),
    value = c(
      "2 Jan 2017 to 23 Dec 2024",
      fmt_int(length(primary_counties)),
      paste0(fmt_median_iqr(altitude_values, 0), "; ", fmt_int(min(altitude_values)), "-", fmt_int(max(altitude_values))),
      fmt_int(nrow(met_panel)),
      fmt_median_iqr(county_context$cold_week_pct, 1)
    ),
    definition = c(
      "ISO-week panel used for descriptive meteorological and surveillance context",
      "Counties represented in the frozen primary laboratory model",
      "Altitude of the county government seat; modelled continuously, with observed range shown",
      "At least five valid daily temperature observations in the week",
      "Calculated from all meteorologically valid weeks, independent of laboratory testing"
    )
  ),
  tibble(
    section = "Descriptive laboratory surveillance",
    characteristic = c(
      "Laboratory-observed county-weeks, n",
      "Laboratory coverage, % of met-valid county-weeks",
      "Laboratory tests, n",
      "Laboratory-positive specimens, n",
      "Pooled laboratory positivity, %",
      "Primary model, county-weeks / tests / positives, n"
    ),
    value = c(
      fmt_int(nrow(lab_observed)),
      fmt_pct(100 * nrow(lab_observed) / nrow(met_panel), 1),
      fmt_int(sum(lab_observed$tests)),
      fmt_int(sum(lab_observed$positive)),
      fmt_pct(100 * sum(lab_observed$positive) / sum(lab_observed$tests), 1),
      paste(fmt_int(nrow(primary)), fmt_int(sum(primary$tests)), fmt_int(sum(primary$positive)), sep = " / ")
    ),
    definition = c(
      "County-weeks with valid temperature and at least one laboratory test",
      "Laboratory-observed weeks divided by meteorologically valid weeks across 45 counties",
      "All valid tests in the descriptive laboratory panel",
      "All positives in the descriptive laboratory panel",
      "Laboratory-positive specimens divided by laboratory tests",
      "Complete-case primary-model sample; 99.5% of descriptive laboratory-observed weeks retained"
    )
  ),
  tibble(
    section = "ILI surveillance",
    characteristic = c(
      "ILI-covered counties, n",
      "ILI-observed county-weeks, n",
      "ILI coverage in 18-county footprint, %",
      "ILI consultations, n",
      "Outpatient visits, n",
      "Pooled ILI consultation rate, %"
    ),
    value = c(
      fmt_int(length(ili_counties)),
      fmt_int(nrow(ili_observed)),
      fmt_pct(100 * nrow(ili_observed) / nrow(ili_met_panel), 1),
      fmt_int(sum(ili_observed$ILI)),
      fmt_int(sum(ili_observed$visits)),
      fmt_pct(100 * sum(ili_observed$ILI) / sum(ili_observed$visits), 2)
    ),
    definition = c(
      "Counties contributing at least one valid ILI county-week",
      "County-weeks with valid ILI consultations and outpatient visits",
      "ILI-observed weeks divided by meteorologically valid weeks in the 18 ILI-covered counties",
      "ILI consultations reported by the ILI surveillance panel",
      "Outpatient visits forming the ILI denominator",
      "ILI consultations divided by outpatient visits"
    )
  )
)

write_csv_utf8(table1, file.path(source_dir, "table1_final_source.csv"))

denominator_audit <- tibble(
  analysis_set = c(
    "Meteorologically valid panel",
    "Descriptive laboratory panel",
    "Frozen primary laboratory model",
    "ILI panel",
    "ILI-footprint meteorological panel"
  ),
  counties = c(
    n_distinct(met_panel$county),
    n_distinct(lab_observed$county),
    n_distinct(primary$county),
    n_distinct(ili_observed$county),
    n_distinct(ili_met_panel$county)
  ),
  county_weeks = c(
    nrow(met_panel),
    nrow(lab_observed),
    nrow(primary),
    nrow(ili_observed),
    nrow(ili_met_panel)
  ),
  tests_or_visits = c(
    NA_real_,
    sum(lab_observed$tests),
    sum(primary$tests),
    sum(ili_observed$visits),
    NA_real_
  ),
  positive_or_ILI = c(
    NA_real_,
    sum(lab_observed$positive),
    sum(primary$positive),
    sum(ili_observed$ILI),
    NA_real_
  )
)
write_csv_utf8(denominator_audit, file.path(source_dir, "denominator_audit.csv"))

# Figure 1 ----------------------------------------------------------------

palette <- c(
  ink = "#25282B",
  muted = "#6D7680",
  line = "#3C4650",
  grid = "#E9EDF0",
  outside = "#ECEFF1",
  blue = "#2C6DA4",
  blue_fill = "#DCEAF3",
  red = "#C85B4B",
  red_fill = "#F1DCD7",
  teal = "#27847F",
  teal_fill = "#DCECE9"
)

theme_map_final <- function(base_size = 7.0) {
  theme_void(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(
        size = base_size + 0.7,
        face = "bold",
        colour = palette[["ink"]],
        margin = margin(b = 3)
      ),
      legend.title = element_text(size = base_size - 0.1, colour = palette[["ink"]]),
      legend.text = element_text(size = base_size - 0.5, colour = palette[["ink"]]),
      plot.margin = margin(3, 5, 2, 5)
    )
}

theme_ts_final <- function(base_size = 7.0, show_x = FALSE) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = palette[["line"]]),
      axis.ticks = element_line(linewidth = 0.25, colour = palette[["line"]]),
      axis.text = element_text(size = base_size - 0.6, colour = palette[["ink"]]),
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      axis.title = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = palette[["grid"]]),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(
        size = base_size + 0.8,
        face = "bold",
        colour = palette[["ink"]],
        margin = margin(b = 2)
      ),
      plot.margin = margin(1, 5, 1, 5)
    )
}

base_map <- ggplot(map_data) +
  geom_sf(fill = palette[["outside"]], colour = "white", linewidth = 0.16) +
  coord_sf(datum = NA)

p_a <- base_map +
  geom_sf(aes(fill = altitude), colour = "white", linewidth = 0.16) +
  scale_fill_gradientn(
    colours = c("#B9E3CF", "#5CB7AC", "#3979AE", "#2B2342"),
    name = "Altitude (m)",
    breaks = c(2000, 3000, 4000),
    na.value = palette[["outside"]],
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(26, "mm"),
      barheight = unit(2.2, "mm")
    )
  ) +
  labs(title = "a  Altitude gradient") +
  theme_map_final() +
  theme(legend.position = "bottom")

p_b <- base_map +
  geom_sf(aes(fill = cold_week_pct), colour = "white", linewidth = 0.16) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#BFD7EA", "#5B91BD", "#174A75"),
    name = "Weeks at or below -15 °C (%)",
    breaks = c(10, 30, 50),
    limits = c(0, 60),
    oob = squish,
    na.value = palette[["outside"]],
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(26, "mm"),
      barheight = unit(2.2, "mm")
    )
  ) +
  labs(title = "b  Absolute cold exposure") +
  theme_map_final() +
  theme(legend.position = "bottom")

p_c <- base_map +
  geom_sf(aes(fill = lab_coverage_pct), colour = "white", linewidth = 0.16) +
  scale_fill_gradientn(
    colours = c("#F7F2E8", "#DFC39F", "#B87848", "#754126"),
    name = "Met-valid weeks with lab data (%)",
    breaks = c(25, 50, 75, 100),
    limits = c(0, 100),
    oob = squish,
    na.value = palette[["outside"]],
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(26, "mm"),
      barheight = unit(2.2, "mm")
    )
  ) +
  labs(title = "c  Laboratory surveillance coverage") +
  theme_map_final() +
  theme(legend.position = "bottom")

year_breaks <- c(as.Date("2017-01-02"), as.Date(paste0(2018:2024, "-01-01")))
x_limits <- as.Date(c("2017-01-02", "2024-12-30"))

make_ts <- function(data, y, label, line_colour, fill_colour, title = NULL, show_x = FALSE) {
  ggplot(data, aes(x = week_start, y = .data[[y]])) +
    geom_area(fill = fill_colour, alpha = 0.55, linewidth = 0, na.rm = TRUE) +
    geom_line(colour = line_colour, linewidth = 0.32, na.rm = TRUE) +
    annotate(
      "text",
      x = as.Date("2017-01-09"),
      y = Inf,
      label = label,
      hjust = 0,
      vjust = 1.25,
      size = 2.25,
      family = "Arial",
      fontface = "bold",
      colour = palette[["ink"]]
    ) +
    scale_x_date(
      breaks = year_breaks,
      date_labels = "%Y",
      limits = x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.035, 0.08))) +
    labs(title = title) +
    theme_ts_final(show_x = show_x)
}

p_d1 <- make_ts(
  weekly_surveillance,
  "laboratory_tests",
  "Laboratory tests",
  palette[["blue"]],
  palette[["blue_fill"]],
  title = "d  Weekly surveillance structure",
  show_x = FALSE
)
p_d2 <- make_ts(
  weekly_surveillance,
  "laboratory_positivity_pct",
  "Laboratory positivity (%)",
  palette[["red"]],
  palette[["red_fill"]],
  show_x = FALSE
)
p_d3 <- make_ts(
  weekly_surveillance,
  "ILI_rate_pct",
  "ILI consultation rate (18 counties, %)",
  palette[["teal"]],
  palette[["teal_fill"]],
  show_x = TRUE
)

top_maps <- p_a | p_b | p_c
time_series <- p_d1 / p_d2 / p_d3 + plot_layout(heights = c(1, 1, 1))

figure1 <- top_maps / time_series +
  plot_layout(heights = c(0.82, 1.05)) +
  plot_annotation(
    title = "Figure 1. Altitude, cold exposure and influenza surveillance structure across Qinghai",
    subtitle = "County-level surveillance context, 2017-2024",
    theme = theme(
      plot.title = element_text(
        family = "Arial",
        face = "bold",
        size = 10.0,
        colour = palette[["ink"]],
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        family = "Arial",
        size = 7.4,
        colour = palette[["muted"]],
        margin = margin(b = 4)
      )
    )
  )

figure_stem <- file.path(
  figure_dir,
  "Figure1_Altitude_cold_exposure_surveillance_FINAL_20260805"
)
width_mm <- 183
height_mm <- 150
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(figure_stem, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure1)
dev.off()

grDevices::cairo_pdf(
  paste0(figure_stem, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Arial",
  bg = "white"
)
print(figure1)
dev.off()

ragg::agg_tiff(
  paste0(figure_stem, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(figure1)
dev.off()

ragg::agg_png(
  paste0(figure_stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(figure1)
dev.off()

legend_lines <- c(
  "**Figure 1. Altitude, cold exposure and influenza surveillance structure across Qinghai.**",
  "**a,** County-government-seat altitude for the 45 counties represented in the frozen primary laboratory model.",
  "**b,** Percentage of all meteorologically valid county-weeks with weekly mean temperature at or below -15 C. The denominator includes all weeks with at least five valid daily temperature observations and is independent of laboratory testing.",
  "**c,** Percentage of meteorologically valid county-weeks with a valid laboratory observation and at least one test. Grey polygons in panels a-c denote county-level areas outside the 45-county primary analytic geography.",
  "**d,** Weekly laboratory tests and pooled laboratory positivity among descriptive laboratory-observed county-weeks, and pooled ILI consultation rate among the 18 ILI-covered counties. Laboratory positivity is the weekly sum of positive specimens divided by the weekly sum of tests; the ILI rate is the weekly sum of ILI consultations divided by the weekly sum of outpatient visits.",
  "The panels describe exposure and surveillance structure and do not estimate a causal temperature effect."
)
write_lines_utf8(
  legend_lines,
  file.path(figure_dir, "Figure1_FINAL_20260805_legend.md")
)

# Formatted Table 1 workbook ----------------------------------------------

workbook <- createWorkbook(creator = "Qinghai influenza study")
addWorksheet(workbook, "Table 1", gridLines = FALSE)
addWorksheet(workbook, "Source and QA", gridLines = FALSE)

title_style <- createStyle(
  fontName = "Arial",
  fontSize = 15,
  fontColour = "#111111",
  textDecoration = "bold",
  fgFill = "#FFFFFF",
  halign = "left",
  valign = "center"
)
subtitle_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#444444",
  fgFill = "#FFFFFF",
  halign = "left",
  valign = "center",
  wrapText = TRUE
)
header_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#111111",
  textDecoration = "bold",
  fgFill = "#FFFFFF",
  halign = "left",
  valign = "center",
  wrapText = TRUE,
  border = c("top", "bottom"),
  borderColour = "#111111",
  borderStyle = "medium"
)
section_style <- createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#111111",
  textDecoration = "bold",
  fgFill = "#F0F0F0",
  halign = "left",
  valign = "center"
)
body_style <- createStyle(
  fontName = "Arial",
  fontSize = 9.5,
  fontColour = "#111111",
  valign = "center",
  wrapText = TRUE,
  border = "bottom",
  borderColour = "#D6DEE3",
  borderStyle = "hair"
)
note_style <- createStyle(
  fontName = "Arial",
  fontSize = 9,
  fontColour = "#333333",
  fgFill = "#FFFFFF",
  valign = "top",
  wrapText = TRUE
)

writeData(
  workbook,
  "Table 1",
  "Table 1. Study setting, meteorological exposure and influenza surveillance characteristics across Qinghai, 2017-2024",
  startRow = 1,
  startCol = 1
)
mergeCells(workbook, "Table 1", cols = 1:2, rows = 1)
addStyle(workbook, "Table 1", title_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setRowHeights(workbook, "Table 1", rows = 1, heights = 30)

writeData(
  workbook,
  "Table 1",
  "Descriptive meteorological, laboratory and ILI panels are reported separately from the frozen complete-case primary model sample.",
  startRow = 2,
  startCol = 1
)
mergeCells(workbook, "Table 1", cols = 1:2, rows = 2)
addStyle(workbook, "Table 1", subtitle_style, rows = 2, cols = 1:2, gridExpand = TRUE)
setRowHeights(workbook, "Table 1", rows = 2, heights = 28)

writeData(
  workbook,
  "Table 1",
  x = data.frame(
    Characteristic = character(),
    Value = character(),
    check.names = FALSE
  ),
  startRow = 4,
  startCol = 1,
  colNames = TRUE
)
addStyle(workbook, "Table 1", header_style, rows = 4, cols = 1:2, gridExpand = TRUE)
setRowHeights(workbook, "Table 1", rows = 4, heights = 28)

current_row <- 5
for (section_name in unique(table1$section)) {
  writeData(workbook, "Table 1", section_name, startRow = current_row, startCol = 1)
  mergeCells(workbook, "Table 1", cols = 1:2, rows = current_row)
  addStyle(workbook, "Table 1", section_style, rows = current_row, cols = 1:2, gridExpand = TRUE)
  setRowHeights(workbook, "Table 1", rows = current_row, heights = 24)
  current_row <- current_row + 1

  section_data <- table1 |>
    filter(section == section_name) |>
    transmute(
      Characteristic = characteristic,
      Value = value
    )
  start_body <- current_row
  writeData(
    workbook,
    "Table 1",
    section_data,
    startRow = start_body,
    startCol = 1,
    colNames = FALSE
  )
  end_body <- start_body + nrow(section_data) - 1
  addStyle(
    workbook,
    "Table 1",
    body_style,
    rows = start_body:end_body,
    cols = 1:2,
    gridExpand = TRUE,
    stack = TRUE
  )
  setRowHeights(workbook, "Table 1", rows = start_body:end_body, heights = 34)
  current_row <- end_body + 1
}

note_row <- current_row + 1
table_note <- paste0(
  "Note: Values are n, %, or median (IQR) as indicated. ILI, influenza-like illness; IQR, interquartile range; met-valid, meteorologically valid. ",
  "Cold-week exposure was calculated from all met-valid county-weeks and was not conditioned on laboratory testing. ",
  "The descriptive laboratory panel required valid temperature and at least one test. The primary model sample additionally applied complete-case and lag-basis requirements. ",
  "For the combined primary-model row, values are county-weeks, tests and positive specimens, respectively; 99.5% of descriptive laboratory-observed weeks were retained. ",
  "ILI coverage was evaluated only within the 18 counties contributing ILI observations. ",
  "One county-week with ILI consultations exceeding outpatient visits was excluded as logically invalid."
)
writeData(workbook, "Table 1", table_note, startRow = note_row, startCol = 1)
mergeCells(workbook, "Table 1", cols = 1:2, rows = note_row)
addStyle(workbook, "Table 1", note_style, rows = note_row, cols = 1:2, gridExpand = TRUE)
setRowHeights(workbook, "Table 1", rows = note_row, heights = 55)

setColWidths(workbook, "Table 1", cols = 1, widths = 72)
setColWidths(workbook, "Table 1", cols = 2, widths = 30)
freezePane(workbook, "Table 1", firstActiveRow = 5)
pageSetup(workbook, "Table 1", orientation = "portrait", fitToWidth = 1, fitToHeight = 0)

writeData(workbook, "Source and QA", denominator_audit, startRow = 1, startCol = 1)
writeData(
  workbook,
  "Source and QA",
  table1,
  startRow = nrow(denominator_audit) + 4,
  startCol = 1
)
addStyle(workbook, "Source and QA", header_style, rows = c(1, nrow(denominator_audit) + 4), cols = 1:5, gridExpand = TRUE)
setColWidths(workbook, "Source and QA", cols = 1:5, widths = c(39, 18, 18, 18, 18))
freezePane(workbook, "Source and QA", firstActiveRow = 2)

table_xlsx <- file.path(
  table_dir,
  "Table1_Study_setting_and_surveillance_FINAL_20260805.xlsx"
)
saveWorkbook(workbook, table_xlsx, overwrite = TRUE)

# Table 1 preview rendered in R -------------------------------------------

table_display <- table1 |>
  select(section, characteristic, value)

display_rows <- list()
for (section_name in unique(table_display$section)) {
  display_rows[[length(display_rows) + 1]] <- tibble(
    row_type = "section",
    characteristic = section_name,
    value = ""
  )
  body <- table_display |>
    filter(section == section_name) |>
    transmute(
      row_type = "body",
      characteristic,
      value
    )
  display_rows[[length(display_rows) + 1]] <- body
}
preview_data <- bind_rows(display_rows)
preview_data$row <- rev(seq_len(nrow(preview_data)))

preview_plot <- ggplot(preview_data) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = row - 0.48, ymax = row + 0.48, fill = row_type),
    colour = "#DEDEDE",
    linewidth = 0.15
  ) +
  geom_text(
    aes(x = 0.015, y = row, label = characteristic, fontface = ifelse(row_type == "section", "bold", "plain")),
    hjust = 0,
    size = 2.35,
    family = "Arial",
    colour = "#111111"
  ) +
  geom_text(
    data = preview_data |> filter(row_type == "body"),
    aes(x = 0.76, y = row, label = value),
    hjust = 0,
    size = 2.35,
    family = "Arial",
    colour = "#111827"
  ) +
  annotate("rect", xmin = 0, xmax = 1, ymin = max(preview_data$row) + 0.92, ymax = max(preview_data$row) + 1.9, fill = "white", colour = "#111111", linewidth = 0.35) +
  annotate("text", x = 0.015, y = max(preview_data$row) + 1.4, label = "Characteristic", hjust = 0, fontface = "bold", family = "Arial", size = 2.7, colour = "#111111") +
  annotate("text", x = 0.76, y = max(preview_data$row) + 1.4, label = "Value", hjust = 0, fontface = "bold", family = "Arial", size = 2.7, colour = "#111111") +
  scale_fill_manual(values = c(section = "#F0F0F0", body = "#FFFFFF"), guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.4, max(preview_data$row) + 2.0), clip = "off") +
  labs(
    title = "Table 1. Study setting, meteorological exposure and influenza surveillance characteristics\nacross Qinghai, 2017-2024",
    subtitle = "Descriptive panels and the frozen primary-model sample are reported with explicit denominators",
    caption = stringr::str_wrap(table_note, width = 175)
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.4, face = "bold", colour = "#111111", lineheight = 1.05, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 7.2, colour = "#444444", margin = margin(b = 7)),
    plot.caption = element_text(size = 5.6, colour = "#333333", hjust = 0, lineheight = 1.05, margin = margin(t = 6)),
    plot.margin = margin(8, 8, 8, 8)
  )

preview_stem <- file.path(
  table_dir,
  "Table1_Study_setting_and_surveillance_FINAL_20260805"
)
ragg::agg_png(
  paste0(preview_stem, ".png"),
  width = 183 / 25.4,
  height = 165 / 25.4,
  units = "in",
  res = 300,
  background = "white"
)
print(preview_plot)
dev.off()

grDevices::cairo_pdf(
  paste0(preview_stem, ".pdf"),
  width = 183 / 25.4,
  height = 165 / 25.4,
  family = "Arial",
  bg = "white"
)
print(preview_plot)
dev.off()

# QA ----------------------------------------------------------------------

qa_checks <- tibble(
  check = c(
    "Primary counties",
    "Primary county-weeks",
    "Primary tests",
    "Primary positives",
    "Descriptive lab positives do not exceed tests",
    "Cold exposure independent of lab filter",
    "ILI counties",
    "Logically invalid ILI county-weeks excluded",
    "Map join complete",
    "Figure year labels end at 2024"
  ),
  expected = c("45", "9758", "77930", "11072", "TRUE", "TRUE", "18", "1", "TRUE", "TRUE"),
  observed = c(
    as.character(n_distinct(primary$county)),
    as.character(nrow(primary)),
    as.character(sum(primary$tests)),
    as.character(sum(primary$positive)),
    as.character(all(lab_observed$positive <= lab_observed$tests)),
    as.character(nrow(met_panel) > nrow(lab_observed)),
    as.character(n_distinct(ili_observed$county)),
    as.character(nrow(invalid_ili)),
    as.character(nrow(unmatched_primary) == 0),
    as.character(max(year_breaks) == as.Date("2024-01-01"))
  )
) |>
  mutate(pass = expected == observed)
write_csv_utf8(qa_checks, file.path(source_dir, "figure1_table1_QA_checks.csv"))

if (!all(qa_checks$pass)) {
  stop("At least one Figure 1 / Table 1 QA check failed.")
}

qa_summary <- c(
  "# Figure 1 and Table 1 final QA",
  "",
  paste0("- Primary model: ", nrow(primary), " county-weeks; ", sum(primary$tests), " tests; ", sum(primary$positive), " positives."),
  paste0("- Meteorological panel: ", nrow(met_panel), " county-weeks across ", n_distinct(met_panel$county), " counties."),
  paste0("- Descriptive laboratory panel: ", nrow(lab_observed), " county-weeks; coverage ", fmt_pct(100 * nrow(lab_observed) / nrow(met_panel), 1), "."),
  paste0("- ILI panel: ", nrow(ili_observed), " county-weeks across ", n_distinct(ili_observed$county), " counties; within-footprint coverage ", fmt_pct(100 * nrow(ili_observed) / nrow(ili_met_panel), 1), "."),
  paste0("- Logically invalid ILI county-weeks excluded: ", nrow(invalid_ili), "."),
  paste0("- Cold exposure threshold: ", cold_threshold, " °C, calculated from all meteorologically valid weeks."),
  "- No primary-analysis counties were unmatched to the shapefile.",
  "- The weekly x-axis uses explicit 2017-2024 breaks and does not display 2025.",
  "- Figure exports: editable SVG/PDF, 600-dpi TIFF, and 300-dpi PNG preview.",
  "- Table exports: formatted XLSX, CSV source, PNG/PDF preview, and denominator audit."
)
write_lines_utf8(qa_summary, file.path(source_dir, "Figure1_Table1_FINAL_QA.md"))

message("Final Figure 1 and Table 1 completed.")
