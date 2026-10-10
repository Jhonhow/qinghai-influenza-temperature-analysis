suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
})

project_root <- normalizePath(Sys.getenv("QH_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir_model <- file.path(project_root, "model_outputs")
dir_panel <- file.path(project_root, "03_panel_construction")

as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  xx <- str_trim(as.character(x))
  xx[xx %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  suppressWarnings(as.Date(parse_date_time(
    xx, orders = c("ymd HMS", "ymd HM", "ymd", "Ymd HMS", "Ymd", "mdy", "dmy"),
    tz = "Asia/Shanghai"
  )))
}
county_alias <- c(
  "大通县" = "大通回族土族自治县", "乐都县" = "乐都区",
  "门源县" = "门源回族自治县", "同仁县" = "同仁市",
  "河南县" = "河南蒙古族自治县", "大柴旦行委" = "大柴旦行政委员会"
)
normalize_county <- function(x) {
  z <- str_squish(as.character(x))
  z <- str_replace(z, "^[]】）)]+", "")
  z <- dplyr::recode(z, !!!county_alias, .default = z)
  z[z %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  z
}
add_iso_week <- function(df, date_col = "date") {
  dt <- df[[date_col]]
  df %>% mutate(
    iso_year = isoyear(dt), iso_week = isoweek(dt),
    week_start = floor_date(dt, unit = "week", week_start = 1),
    week_id = sprintf("%04d-W%02d", iso_year, iso_week)
  )
}
write_csv_utf8 <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

version_id <- "v4_t2m_20260825"
out_dir <- file.path(dir_model, paste0("reanalysis_", version_id))
audit_dir <- file.path(dir_panel, paste0("reanalysis_", version_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

old_panel_path <- file.path(
  dir_model, "reanalysis_v3_20260728", "county_week_panel_v3_2016_2025.rds"
)
new_met_path <- Sys.getenv(
  "QH_MET_WORKBOOK",
  unset = file.path(project_root, "data", "updated_met_t2m_pressure.xlsx")
)
if (!file.exists(old_panel_path)) stop("Missing frozen v3 panel: ", old_panel_path)
if (!file.exists(new_met_path)) stop("Missing updated meteorological workbook: ", new_met_path)

first_nonmissing <- function(x) {
  z <- x[!is.na(x) & x != ""]
  if (length(z)) z[[1]] else NA_character_
}

old_panel <- readRDS(old_panel_path)
met_raw <- readxl::read_excel(new_met_path, sheet = "Qinghai_PREC_Long")

required_met <- c("date", "City", "County")
if (!all(required_met %in% names(met_raw))) {
  stop("Updated meteorological workbook lacks: ",
       paste(setdiff(required_met, names(met_raw)), collapse = ", "))
}

met_daily <- met_raw %>%
  transmute(
    date = parse_date_safe(.data[["date"]]),
    city = normalize_county(.data[["City"]]),
    county = normalize_county(.data[["County"]]),
    # Position-based access prevents Windows locale conversion from changing
    # the five Chinese header labels in the merged workbook.
    temperature = as_num(met_raw[[14]]),
    relative_humidity = as_num(met_raw[[10]]),
    pressure = as_num(met_raw[[15]]),
    precipitation = as_num(met_raw[[9]]),
    wind_speed = as_num(met_raw[[12]])
  ) %>%
  filter(!is.na(date), !is.na(county)) %>%
  mutate(
    # The updated pressure is already expressed in hPa.
    pressure_hpa = pressure,
    saturation_vapor_pressure_hpa = 6.112 * exp(
      (17.67 * temperature) / (temperature + 243.5)
    ),
    vapor_pressure_hpa = relative_humidity / 100 * saturation_vapor_pressure_hpa,
    absolute_humidity_g_m3 = 216.7 * vapor_pressure_hpa / (temperature + 273.15),
    specific_humidity_g_kg = 1000 * (
      0.622 * vapor_pressure_hpa /
        (pressure_hpa - 0.378 * vapor_pressure_hpa)
    )
  ) %>%
  add_iso_week("date")

if (anyDuplicated(met_daily[c("county", "date")]) > 0) {
  stop("Duplicate county-date keys in updated meteorological data.")
}

met_week <- met_daily %>%
  group_by(county, iso_year, iso_week, week_start, week_id) %>%
  summarise(
    city_met = first_nonmissing(city),
    met_days = n_distinct(date),
    temp_valid_days = sum(!is.na(temperature)),
    rh_valid_days = sum(!is.na(relative_humidity)),
    pressure_valid_days = sum(!is.na(pressure)),
    precip_valid_days = sum(!is.na(precipitation)),
    wind_valid_days = sum(!is.na(wind_speed)),
    ah_valid_days = sum(!is.na(absolute_humidity_g_m3)),
    q_valid_days = sum(!is.na(specific_humidity_g_kg)),
    temperature = mean(temperature, na.rm = TRUE),
    relative_humidity = mean(relative_humidity, na.rm = TRUE),
    pressure = mean(pressure, na.rm = TRUE),
    pressure_hpa = mean(pressure_hpa, na.rm = TRUE),
    precipitation = sum(precipitation, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    vapor_pressure_hpa = mean(vapor_pressure_hpa, na.rm = TRUE),
    absolute_humidity_g_m3 = mean(absolute_humidity_g_m3, na.rm = TRUE),
    specific_humidity_g_kg = mean(specific_humidity_g_kg, na.rm = TRUE),
    .groups = "drop"
  )

weather_cols <- c(
  "met_days", "temp_valid_days", "rh_valid_days", "pressure_valid_days",
  "precip_valid_days", "wind_valid_days", "ah_valid_days", "q_valid_days",
  "temperature", "relative_humidity", "pressure", "pressure_hpa",
  "precipitation", "wind_speed", "vapor_pressure_hpa",
  "absolute_humidity_g_m3", "specific_humidity_g_kg", "met_available",
  "met_valid_5day_main", "met_valid_4day_sensitivity",
  "county_week_temperature_climatology", "temperature_percentile_county",
  "temperature_anomaly", "temperature_anomaly_z"
)

panel <- old_panel %>%
  select(-any_of(weather_cols)) %>%
  left_join(
    met_week,
    by = c("county", "iso_year", "iso_week", "week_start", "week_id")
  ) %>%
  mutate(
    city = coalesce(city, city_met),
    met_available = !is.na(temperature),
    met_valid_5day_main = temp_valid_days >= 5 &
      rh_valid_days >= 5 & pressure_valid_days >= 5 &
      precip_valid_days >= 5 & wind_valid_days >= 5,
    met_valid_4day_sensitivity = temp_valid_days >= 4 &
      rh_valid_days >= 4 & pressure_valid_days >= 4 &
      precip_valid_days >= 4 & wind_valid_days >= 4,
    temperature_definition = "weekly mean of daily mean 2-m air temperature",
    pressure_definition = "weekly mean of daily mean surface pressure (hPa)"
  ) %>%
  select(-city_met) %>%
  arrange(county, week_start)

temperature_climatology <- panel %>%
  filter(iso_year >= 2016, iso_year <= 2024, met_available) %>%
  group_by(county, iso_week) %>%
  summarise(
    county_week_temperature_climatology = mean(temperature),
    .groups = "drop"
  )

panel <- panel %>%
  left_join(temperature_climatology, by = c("county", "iso_week")) %>%
  group_by(county) %>%
  mutate(
    temperature_percentile_county =
      (rank(temperature, ties.method = "average", na.last = "keep") - 0.5) /
      sum(!is.na(temperature)),
    temperature_anomaly = temperature - county_week_temperature_climatology,
    temperature_anomaly_z = temperature_anomaly / sd(temperature_anomaly, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(county, week_start)

if (nrow(panel) != nrow(old_panel)) stop("Panel row count changed.")
if (anyDuplicated(panel[c("county", "iso_year", "iso_week")]) > 0) {
  stop("Duplicate county-week keys detected in v4 panel.")
}

key <- c("county", "iso_year", "iso_week")
paired <- old_panel %>%
  select(all_of(key), old_temperature = temperature, old_pressure_hpa = pressure_hpa) %>%
  inner_join(
    panel %>% select(all_of(key), new_temperature = temperature,
                     new_pressure_hpa = pressure_hpa, altitude, week_start),
    by = key
  ) %>%
  filter(iso_year >= 2017, iso_year <= 2024)

main_lab <- panel %>%
  filter(
    iso_year >= 2017, iso_year <= 2024, lab_available, tests > 0,
    positive <= tests, met_valid_5day_main, !is.na(temperature), !is.na(altitude)
  )

paired_main <- paired %>%
  semi_join(main_lab, by = key) %>%
  mutate(
    difference_new_minus_old = new_temperature - old_temperature,
    altitude_band = cut(
      altitude,
      breaks = c(-Inf, quantile(main_lab$altitude, c(.1, .5, .9)), Inf),
      include.lowest = TRUE, dig.lab = 5
    ),
    month = month(week_start)
  )

comparison_overall <- paired_main %>%
  summarise(
    county_weeks = n(),
    old_mean = mean(old_temperature),
    new_mean = mean(new_temperature),
    mean_difference = mean(difference_new_minus_old),
    median_difference = median(difference_new_minus_old),
    rmse = sqrt(mean(difference_new_minus_old^2)),
    correlation = cor(old_temperature, new_temperature),
    old_min = min(old_temperature), old_max = max(old_temperature),
    new_min = min(new_temperature), new_max = max(new_temperature)
  )

comparison_by_month <- paired_main %>%
  group_by(month) %>%
  summarise(
    n = n(), old_mean = mean(old_temperature), new_mean = mean(new_temperature),
    mean_difference = mean(difference_new_minus_old),
    correlation = cor(old_temperature, new_temperature), .groups = "drop"
  )

comparison_by_altitude <- paired_main %>%
  group_by(altitude_band) %>%
  summarise(
    n = n(), altitude_median = median(altitude),
    old_mean = mean(old_temperature), new_mean = mean(new_temperature),
    mean_difference = mean(difference_new_minus_old),
    correlation = cor(old_temperature, new_temperature), .groups = "drop"
  )

temperature_quantiles <- bind_rows(
  tibble(source = "old_land_surface_temperature", value = paired_main$old_temperature),
  tibble(source = "new_2m_air_temperature", value = paired_main$new_temperature)
) %>%
  group_by(source) %>%
  summarise(
    n = sum(!is.na(value)), min = min(value, na.rm = TRUE),
    p01 = quantile(value, .01, na.rm = TRUE),
    p05 = quantile(value, .05, na.rm = TRUE),
    p10 = quantile(value, .10, na.rm = TRUE),
    p25 = quantile(value, .25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    p75 = quantile(value, .75, na.rm = TRUE),
    p90 = quantile(value, .90, na.rm = TRUE),
    p95 = quantile(value, .95, na.rm = TRUE),
    p99 = quantile(value, .99, na.rm = TRUE), max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

thresholds <- c(-25, -20, -15, -10, -5, 0, 5)
threshold_coverage <- bind_rows(lapply(thresholds, function(z) {
  bind_rows(
    tibble(source = "old_land_surface_temperature", threshold = z,
           proportion = mean(paired_main$old_temperature <= z)),
    tibble(source = "new_2m_air_temperature", threshold = z,
           proportion = mean(paired_main$new_temperature <= z))
  )
}))

altitude_anchors <- as.numeric(quantile(main_lab$altitude, c(.1, .5, .9)))
local_support <- bind_rows(lapply(altitude_anchors, function(a) {
  z <- main_lab %>% filter(abs(altitude - a) <= 250)
  tibble(
    altitude_anchor = a, altitude_window_m = 250, n = nrow(z),
    min = min(z$temperature), p01 = quantile(z$temperature, .01),
    p05 = quantile(z$temperature, .05), p10 = quantile(z$temperature, .10),
    median = median(z$temperature), p90 = quantile(z$temperature, .90),
    p95 = quantile(z$temperature, .95), p99 = quantile(z$temperature, .99),
    max = max(z$temperature)
  )
}))

county_threshold_coverage <- bind_rows(lapply(thresholds, function(z) {
  main_lab %>%
    group_by(county, city, altitude) %>%
    summarise(
      threshold = z, valid_weeks = n(), weeks_at_or_below = sum(temperature <= z),
      percentage = 100 * weeks_at_or_below / valid_weeks, .groups = "drop"
    )
}))

coverage_summary <- tibble(
  analytic_panel = "v4_primary_lab_2017_2024_met5_routine_outpatient",
  county_weeks = nrow(main_lab), counties = n_distinct(main_lab$county),
  start = min(main_lab$week_start), end = max(main_lab$week_start),
  total_tests = sum(main_lab$tests), total_positive = sum(main_lab$positive),
  total_fluA_only = sum(main_lab$fluA_only),
  total_fluB_only = sum(main_lab$fluB_only)
)

saveRDS(panel, file.path(out_dir, "county_week_panel_v4_t2m_2016_2025.rds"))
write_csv_utf8(panel, file.path(out_dir, "county_week_panel_v4_t2m_2016_2025.csv"))
write_csv_utf8(comparison_overall, file.path(audit_dir, "v4_old_new_temperature_comparison_overall.csv"))
write_csv_utf8(comparison_by_month, file.path(audit_dir, "v4_old_new_temperature_comparison_by_month.csv"))
write_csv_utf8(comparison_by_altitude, file.path(audit_dir, "v4_old_new_temperature_comparison_by_altitude.csv"))
write_csv_utf8(temperature_quantiles, file.path(audit_dir, "v4_temperature_quantiles.csv"))
write_csv_utf8(threshold_coverage, file.path(audit_dir, "v4_threshold_coverage_overall.csv"))
write_csv_utf8(local_support, file.path(audit_dir, "v4_local_temperature_support_by_altitude.csv"))
write_csv_utf8(county_threshold_coverage, file.path(audit_dir, "v4_county_threshold_coverage.csv"))
write_csv_utf8(coverage_summary, file.path(audit_dir, "v4_panel_coverage_summary.csv"))

message("Completed v4 T2m panel and exposure audit: ", out_dir)
