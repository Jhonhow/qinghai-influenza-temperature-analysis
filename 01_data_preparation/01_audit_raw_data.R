source(file.path("R", "00_config.R"))

sheet_inventory <- function(path) {
  sheets <- excel_sheets(path)
  bind_rows(lapply(sheets, function(sh) {
    dat <- read_excel(path, sheet = sh, n_max = 0)
    tibble(
      file = basename(path),
      sheet = sh,
      n_cols = ncol(dat),
      columns = paste(names(dat), collapse = " | ")
    )
  }))
}

file_inventory <- bind_rows(lapply(
  c(file_lab, file_ili, file_met, file_alt, file_gdp, file_pop_density, file_age),
  sheet_inventory
))
write_csv_utf8(file_inventory, file.path(dir_audit, "raw_file_sheet_inventory.csv"))

write_csv_utf8(
  tibble(raw_name = names(county_alias), standardized_name = unname(county_alias)),
  file.path(dir_audit, "county_name_aliases_applied.csv")
)

lab_raw <- read_excel(file_lab, sheet = "总表")
ili_raw <- read_excel(file_ili, sheet = "哨点医院")
met_raw <- read_excel(file_met, sheet = 1)
alt_raw <- read_excel(file_alt, sheet = "Sheet1")
gdp_raw <- read_excel(file_gdp, sheet = 1)
pop_raw <- read_excel(file_pop_density, sheet = 1)
age_raw <- read_excel(file_age, sheet = 1, col_names = FALSE)

lab <- lab_raw %>%
  mutate(
    date = parse_date_safe(.data[["发病日期"]]),
    county = normalize_county(.data[["区县匹配"]]),
    city = normalize_county(.data[["地市匹配"]]),
    positive_raw = as_num(.data[["阴阳性"]])
  )

ili <- ili_raw %>%
  mutate(
    date = parse_date_safe(.data[["时间"]]),
    county = normalize_county(.data[["区县"]]),
    city = normalize_county(.data[["地市州"]]),
    ili_cases = as_num(.data[["ILI合计"]]),
    visits = as_num(.data[["门诊就诊总数"]])
  )

met <- met_raw %>%
  mutate(
    date = parse_date_safe(.data[["date"]]),
    county = normalize_county(.data[["County"]]),
    city = normalize_county(.data[["City"]]),
    temperature = as_num(.data[["温度"]]),
    relative_humidity = as_num(.data[["相对湿度"]]),
    pressure = as_num(.data[["表面气压"]]),
    precipitation = as_num(.data[["降水量"]]),
    wind_speed = as_num(.data[["风速"]])
  )

alt <- alt_raw %>%
  transmute(
    county = normalize_county(.data[["区县"]]),
    altitude = as_num(.data[["海拔"]])
  ) %>%
  filter(!is.na(county))

gdp <- gdp_raw %>%
  transmute(county = normalize_county(.data[["区县"]])) %>%
  filter(!is.na(county))

pop_density <- pop_raw %>%
  transmute(county = normalize_county(.data[["区县"]])) %>%
  filter(!is.na(county))

age_counties <- age_raw[-c(1, 2), ][[1]]
age <- tibble(county = normalize_county(age_counties)) %>%
  filter(!is.na(county), county != "青海省")

source_overview <- bind_rows(
  lab %>% summarise(
    source = "laboratory",
    records = n(),
    counties = n_distinct(county, na.rm = TRUE),
    date_min = min(date, na.rm = TRUE),
    date_max = max(date, na.rm = TRUE),
    missing_county = sum(is.na(county)),
    missing_date = sum(is.na(date)),
    total_tests_or_visits = n(),
    total_positive_or_ili = sum(positive_raw == 1, na.rm = TRUE)
  ),
  ili %>% summarise(
    source = "ILI_sentinel",
    records = n(),
    counties = n_distinct(county, na.rm = TRUE),
    date_min = min(date, na.rm = TRUE),
    date_max = max(date, na.rm = TRUE),
    missing_county = sum(is.na(county)),
    missing_date = sum(is.na(date)),
    total_tests_or_visits = sum(visits, na.rm = TRUE),
    total_positive_or_ili = sum(ili_cases, na.rm = TRUE)
  ),
  met %>% summarise(
    source = "meteorology_daily",
    records = n(),
    counties = n_distinct(county, na.rm = TRUE),
    date_min = min(date, na.rm = TRUE),
    date_max = max(date, na.rm = TRUE),
    missing_county = sum(is.na(county)),
    missing_date = sum(is.na(date)),
    total_tests_or_visits = NA_real_,
    total_positive_or_ili = NA_real_
  )
)
write_csv_utf8(source_overview, file.path(dir_audit, "source_overview.csv"))

missingness <- bind_rows(
  tibble(
    source = "laboratory",
    variable = c("发病日期", "区县匹配", "阴阳性", "分型结果"),
    missing_n = c(sum(is.na(lab$date)), sum(is.na(lab$county)),
                  sum(is.na(lab[["阴阳性"]])), sum(is.na(lab[["分型结果"]]))),
    total_n = nrow(lab)
  ),
  tibble(
    source = "ILI_sentinel",
    variable = c("时间", "区县", "ILI合计", "门诊就诊总数"),
    missing_n = c(sum(is.na(ili$date)), sum(is.na(ili$county)),
                  sum(is.na(ili$ili_cases)), sum(is.na(ili$visits))),
    total_n = nrow(ili)
  ),
  tibble(
    source = "meteorology_daily",
    variable = c("date", "County", "温度", "相对湿度", "表面气压", "降水量", "风速"),
    missing_n = c(sum(is.na(met$date)), sum(is.na(met$county)),
                  sum(is.na(met$temperature)), sum(is.na(met$relative_humidity)),
                  sum(is.na(met$pressure)), sum(is.na(met$precipitation)),
                  sum(is.na(met$wind_speed))),
    total_n = nrow(met)
  )
) %>%
  mutate(missing_pct = missing_n / total_n)
write_csv_utf8(missingness, file.path(dir_audit, "key_variable_missingness.csv"))

lab_result_distribution <- lab %>%
  count(分型结果, 阴阳性, name = "n", sort = TRUE)
write_csv_utf8(lab_result_distribution, file.path(dir_audit, "laboratory_result_distribution.csv"))

county_sets <- tibble(
  county = sort(unique(na.omit(c(
    lab$county, ili$county, met$county, alt$county,
    gdp$county, pop_density$county, age$county
  ))))
) %>%
  mutate(
    in_lab = county %in% lab$county,
    in_ili = county %in% ili$county,
    in_meteorology = county %in% met$county,
    in_altitude = county %in% alt$county,
    in_gdp = county %in% gdp$county,
    in_population_density = county %in% pop_density$county,
    in_age_structure = county %in% age$county
  )
write_csv_utf8(county_sets, file.path(dir_audit, "county_source_presence_matrix.csv"))

date_failures <- bind_rows(
  lab %>%
    filter(is.na(date) & !is.na(.data[["发病日期"]])) %>%
    transmute(source = "laboratory", raw_date = as.character(.data[["发病日期"]]), county),
  ili %>%
    filter(is.na(date) & !is.na(.data[["时间"]])) %>%
    transmute(source = "ILI_sentinel", raw_date = as.character(.data[["时间"]]), county),
  met %>%
    filter(is.na(date) & !is.na(.data[["date"]])) %>%
    transmute(source = "meteorology_daily", raw_date = as.character(.data[["date"]]), county)
)
write_csv_utf8(date_failures, file.path(dir_audit, "date_parse_failures.csv"))

county_coverage <- bind_rows(
  lab %>%
    filter(!is.na(county), !is.na(date)) %>%
    group_by(source = "laboratory", county) %>%
    summarise(records = n(), date_min = min(date), date_max = max(date), .groups = "drop"),
  ili %>%
    filter(!is.na(county), !is.na(date)) %>%
    group_by(source = "ILI_sentinel", county) %>%
    summarise(records = n(), date_min = min(date), date_max = max(date), .groups = "drop"),
  met %>%
    filter(!is.na(county), !is.na(date)) %>%
    group_by(source = "meteorology_daily", county) %>%
    summarise(records = n(), date_min = min(date), date_max = max(date), .groups = "drop")
)
write_csv_utf8(county_coverage, file.path(dir_audit, "county_date_coverage_by_source.csv"))

report <- c(
  "# Preliminary Data Audit",
  "",
  paste0("Data directory: `", basename(data_dir), "`"),
  "",
  "## Core interpretation",
  "",
  "- Laboratory data are individual/record-level and can be aggregated to county-ISO-week using onset date.",
  "- ILI sentinel data are daily and can be aggregated to county-ISO-week using the date column.",
  "- Meteorological data are daily county-level observations, so weekly aggregation is methodologically acceptable.",
  "- Altitude, population density, GDP and age structure are county-level or county-year attributes and should not be treated as time-varying daily exposures.",
  "",
  "## Generated audit files",
  "",
  "- `raw_file_sheet_inventory.csv`",
  "- `source_overview.csv`",
  "- `key_variable_missingness.csv`",
  "- `laboratory_result_distribution.csv`",
  "- `county_source_presence_matrix.csv`",
  "- `date_parse_failures.csv`",
  "- `county_date_coverage_by_source.csv`",
  "",
  "## Boundary for next steps",
  "",
  "The next script should construct outcome-specific county-week panels. It should not force laboratory and ILI counties to the same intersection, and it should keep weekly meteorological valid-day counts for the main 5-day threshold and 4-day sensitivity analysis."
)
write_lines_utf8(report, file.path(dir_audit, "preliminary_data_audit_report.md"))

message("Audit complete. Outputs written to: ", dir_audit)
