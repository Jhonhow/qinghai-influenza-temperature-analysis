source(file.path("R", "00_config.R"))

updated_lab_path <- file.path(
  project_root,
  "2016-2025流感核酸监测信息20260703_已匹配区县.xlsx"
)
previous_lab_path <- file.path(
  data_dir,
  "2016-2025流感核酸监测信息20260703_已匹配区县.xlsx"
)
out_dir <- file.path(dir_audit, "lab_subtype_update_20260727")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

updated_raw <- read_excel(updated_lab_path, sheet = "总表")
previous_raw <- read_excel(previous_lab_path, sheet = "总表")

clean_text <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

classify_previous_subtype <- function(x) {
  x <- clean_text(x)
  case_when(
    str_detect(x, regex("甲型|甲流|H1|H3|A型|A\\(", ignore_case = TRUE)) ~ 1L,
    str_detect(x, regex("乙型|乙流|Victoria|Yamagata|B型|B\\(", ignore_case = TRUE)) ~ 2L,
    TRUE ~ NA_integer_
  )
}

updated <- updated_raw %>%
  transmute(
    in_updated = TRUE,
    record_id = as.character(.data[["序号"]]),
    date = parse_date_safe(.data[["发病日期"]]),
    reported_year = as_num(.data[["年"]]),
    reported_week = as_num(.data[["周"]]),
    reported_year_week = clean_text(.data[["年周"]]),
    city = normalize_county(.data[["地市匹配"]]),
    county = normalize_county(.data[["区县匹配"]]),
    address_code = clean_text(.data[["现住址国标"]]),
    specimen_type = clean_text(.data[["标本种类"]]),
    specimen_source = clean_text(.data[["标本来源"]]),
    receiving_unit = clean_text(.data[["收样单位"]]),
    raw_pcr_result = clean_text(.data[["原始标本PCR鉴定结果"]]),
    subtype_raw = as_num(.data[["甲乙流分型（1=甲流，2=乙流）"]]),
    positive_raw = as_num(.data[["阴阳性"]])
  ) %>%
  mutate(
    positive = case_when(
      positive_raw == 1 ~ 1L,
      positive_raw == 0 ~ 0L,
      str_detect(raw_pcr_result, "阴") ~ 0L,
      str_detect(
        raw_pcr_result,
        regex("阳|甲|乙|H1|H3|Victoria|Yamagata|A型|B型", ignore_case = TRUE)
      ) ~ 1L,
      TRUE ~ NA_integer_
    ),
    fluA = if_else(positive == 1L & subtype_raw == 1, 1L, 0L),
    fluB = if_else(positive == 1L & subtype_raw == 2, 1L, 0L),
    typed_positive = positive == 1L & subtype_raw %in% c(1, 2),
    iso_year_derived = isoyear(date),
    iso_week_derived = isoweek(date),
    epi_year_derived = epiyear(date),
    epi_week_derived = epiweek(date),
    week_start = floor_date(date, unit = "week", week_start = 1)
  )

previous <- previous_raw %>%
  transmute(
    in_previous = TRUE,
    record_id = as.character(.data[["序号"]]),
    date = parse_date_safe(.data[["发病日期"]]),
    reported_year = as_num(.data[["年"]]),
    city = normalize_county(.data[["地市匹配"]]),
    county = normalize_county(.data[["区县匹配"]]),
    address_code = clean_text(.data[["现住址国标"]]),
    previous_result = clean_text(.data[["分型结果"]]),
    positive_raw = as_num(.data[["阴阳性"]]),
    previous_subtype = classify_previous_subtype(.data[["分型结果"]])
  )

updated_id <- updated %>% distinct(record_id, .keep_all = TRUE)
previous_id <- previous %>% distinct(record_id, .keep_all = TRUE)

version_comparison <- full_join(
  previous_id %>%
    rename(
      previous_date = date,
      previous_city = city,
      previous_county = county,
      previous_address_code = address_code,
      previous_positive = positive_raw
    ),
  updated_id %>%
    rename(
      updated_date = date,
      updated_city = city,
      updated_county = county,
      updated_address_code = address_code,
      updated_positive = positive_raw,
      updated_subtype = subtype_raw
    ),
  by = "record_id"
) %>%
  mutate(
    in_previous = coalesce(in_previous, FALSE),
    in_updated = coalesce(in_updated, FALSE),
    record_presence = case_when(
      in_previous & in_updated ~ "both",
      in_previous ~ "previous_only",
      TRUE ~ "updated_only"
    ),
    date_changed = record_presence == "both" & (
      is.na(previous_date) != is.na(updated_date) |
        (!is.na(previous_date) & !is.na(updated_date) & previous_date != updated_date)
    ),
    city_changed = record_presence == "both" & (
      is.na(previous_city) != is.na(updated_city) |
        (!is.na(previous_city) & !is.na(updated_city) & previous_city != updated_city)
    ),
    county_changed = record_presence == "both" & (
      is.na(previous_county) != is.na(updated_county) |
        (!is.na(previous_county) & !is.na(updated_county) & previous_county != updated_county)
    ),
    positive_changed = record_presence == "both" & (
      is.na(previous_positive) != is.na(updated_positive) |
        (!is.na(previous_positive) & !is.na(updated_positive) &
          previous_positive != updated_positive)
    ),
    subtype_changed = record_presence == "both" &
      !is.na(previous_subtype) &
      !is.na(updated_subtype) &
      previous_subtype != updated_subtype
  )

audit_summary <- tibble::tribble(
  ~metric, ~value, ~interpretation,
  "updated_rows", nrow(updated), "Rows in the updated laboratory workbook",
  "previous_rows", nrow(previous), "Rows in the prior laboratory workbook",
  "updated_columns", ncol(updated_raw), "Columns in the updated workbook",
  "previous_columns", ncol(previous_raw), "Columns in the prior workbook",
  "updated_duplicate_record_ids", sum(duplicated(updated$record_id)), "Duplicate 序号 values",
  "previous_duplicate_record_ids", sum(duplicated(previous$record_id)), "Duplicate 序号 values",
  "updated_missing_dates", sum(is.na(updated$date)), "Records without parsed onset date",
  "updated_missing_counties", sum(is.na(updated$county)), "Records without matched county",
  "updated_missing_positive", sum(is.na(updated$positive)), "Records without classifiable positivity",
  "updated_positive", sum(updated$positive == 1, na.rm = TRUE), "Influenza-positive records",
  "updated_negative", sum(updated$positive == 0, na.rm = TRUE), "Influenza-negative records",
  "updated_fluA", sum(updated$fluA, na.rm = TRUE), "Influenza A-positive records",
  "updated_fluB", sum(updated$fluB, na.rm = TRUE), "Influenza B-positive records",
  "positive_untyped", sum(updated$positive == 1 & is.na(updated$subtype_raw), na.rm = TRUE),
    "Positive records without type 1/2",
  "negative_with_subtype", sum(updated$positive == 0 & updated$subtype_raw %in% c(1, 2), na.rm = TRUE),
    "Negative records carrying a type code",
  "invalid_subtype_code", sum(!is.na(updated$subtype_raw) & !updated$subtype_raw %in% c(0, 1, 2)),
    "Subtype values outside expected 0/1/2 coding",
  "subtype_without_positive", sum(updated$subtype_raw %in% c(1, 2) & updated$positive != 1, na.rm = TRUE),
    "Typed records not coded positive",
  "reported_iso_year_mismatch", sum(
    !is.na(updated$reported_year) &
      !is.na(updated$iso_year_derived) &
      updated$reported_year != updated$iso_year_derived
  ), "Reported year differs from ISO year derived from onset date",
  "reported_iso_week_mismatch", sum(
    !is.na(updated$reported_week) &
      !is.na(updated$iso_week_derived) &
      updated$reported_week != updated$iso_week_derived
  ), "Reported week differs from ISO week derived from onset date",
  "reported_epi_year_mismatch", sum(
    !is.na(updated$reported_year) &
      !is.na(updated$epi_year_derived) &
      updated$reported_year != updated$epi_year_derived
  ), "Reported year differs from epidemiological year derived from onset date",
  "reported_epi_week_mismatch", sum(
    !is.na(updated$reported_week) &
      !is.na(updated$epi_week_derived) &
      updated$reported_week != updated$epi_week_derived
  ), "Reported week differs from epidemiological week derived from onset date",
  "records_in_both_versions", sum(version_comparison$record_presence == "both"),
    "Unique 序号 values shared by both versions",
  "records_previous_only", sum(version_comparison$record_presence == "previous_only"),
    "Unique 序号 values only in the prior version",
  "records_updated_only", sum(version_comparison$record_presence == "updated_only"),
    "Unique 序号 values only in the updated version",
  "shared_date_changes", sum(version_comparison$date_changed, na.rm = TRUE),
    "Shared records with changed onset date",
  "shared_county_changes", sum(version_comparison$county_changed, na.rm = TRUE),
    "Shared records with changed matched county",
  "shared_positive_changes", sum(version_comparison$positive_changed, na.rm = TRUE),
    "Shared records with changed positivity",
  "shared_subtype_disagreements", sum(version_comparison$subtype_changed, na.rm = TRUE),
    "Shared typed records whose old text and new code disagree"
) %>%
  mutate(value = as.character(value))

subtype_distribution <- updated %>%
  mutate(
    subtype_label = case_when(
      subtype_raw == 1 ~ "Influenza A",
      subtype_raw == 2 ~ "Influenza B",
      is.na(subtype_raw) ~ "Missing",
      TRUE ~ paste0("Other: ", subtype_raw)
    )
  ) %>%
  count(positive_raw, positive, subtype_label, name = "n", sort = TRUE)

pcr_distribution <- updated %>%
  count(raw_pcr_result, positive_raw, subtype_raw, name = "n", sort = TRUE)

old_new_subtype_crosswalk <- version_comparison %>%
  filter(record_presence == "both") %>%
  count(
    previous_result,
    previous_positive,
    previous_subtype,
    raw_pcr_result,
    updated_positive,
    updated_subtype,
    name = "n",
    sort = TRUE
  )

version_differences <- version_comparison %>%
  filter(
    record_presence != "both" |
      date_changed |
      city_changed |
      county_changed |
      positive_changed |
      subtype_changed
  ) %>%
  select(
    record_id,
    record_presence,
    previous_date,
    updated_date,
    previous_city,
    updated_city,
    previous_county,
    updated_county,
    previous_positive,
    updated_positive,
    previous_result,
    previous_subtype,
    raw_pcr_result,
    updated_subtype,
    date_changed,
    city_changed,
    county_changed,
    positive_changed,
    subtype_changed
  )

reported_week_mismatch_examples <- updated %>%
  filter(
    (!is.na(reported_year) & reported_year != iso_year_derived) |
      (!is.na(reported_week) & reported_week != iso_week_derived)
  ) %>%
  select(
    record_id,
    date,
    reported_year,
    reported_week,
    reported_year_week,
    iso_year_derived,
    iso_week_derived,
    epi_year_derived,
    epi_week_derived,
    county
  ) %>%
  slice_head(n = 200)

yearly_coverage <- updated %>%
  filter(!is.na(date)) %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    tests = n(),
    positive = sum(positive == 1, na.rm = TRUE),
    influenza_A = sum(fluA, na.rm = TRUE),
    influenza_B = sum(fluB, na.rm = TRUE),
    positive_untyped = sum(positive == 1 & is.na(subtype_raw), na.rm = TRUE),
    positivity = positive / tests,
    typed_fraction_among_positive = if_else(
      positive > 0,
      (influenza_A + influenza_B) / positive,
      NA_real_
    ),
    .groups = "drop"
  )

county_coverage <- updated %>%
  filter(!is.na(county), !is.na(date), year(date) %in% 2017:2024) %>%
  group_by(city, county) %>%
  summarise(
    tests = n(),
    positive = sum(positive == 1, na.rm = TRUE),
    influenza_A = sum(fluA, na.rm = TRUE),
    influenza_B = sum(fluB, na.rm = TRUE),
    positive_untyped = sum(positive == 1 & is.na(subtype_raw), na.rm = TRUE),
    weeks_tested = n_distinct(week_start),
    weeks_with_A = n_distinct(week_start[fluA == 1]),
    weeks_with_B = n_distinct(week_start[fluB == 1]),
    typed_fraction_among_positive = if_else(
      positive > 0,
      (influenza_A + influenza_B) / positive,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(tests))

county_week_subtypes <- updated %>%
  filter(
    !is.na(county),
    !is.na(date),
    iso_year_derived %in% 2017:2024
  ) %>%
  group_by(
    county,
    city,
    iso_year = iso_year_derived,
    iso_week = iso_week_derived,
    week_start
  ) %>%
  summarise(
    tests = n(),
    positive = sum(positive == 1, na.rm = TRUE),
    influenza_A = sum(fluA, na.rm = TRUE),
    influenza_B = sum(fluB, na.rm = TRUE),
    positive_untyped = sum(positive == 1 & is.na(subtype_raw), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    positivity = positive / tests,
    influenza_A_positivity = influenza_A / tests,
    influenza_B_positivity = influenza_B / tests
  )

weekly_feasibility <- county_week_subtypes %>%
  summarise(
    county_weeks = n(),
    counties = n_distinct(county),
    total_tests = sum(tests),
    total_positive = sum(positive),
    total_influenza_A = sum(influenza_A),
    total_influenza_B = sum(influenza_B),
    county_weeks_with_positive = sum(positive > 0),
    county_weeks_with_A = sum(influenza_A > 0),
    county_weeks_with_B = sum(influenza_B > 0),
    zero_fraction_positive = mean(positive == 0),
    zero_fraction_A = mean(influenza_A == 0),
    zero_fraction_B = mean(influenza_B == 0)
  )

existing_panel <- readRDS(
  file.path(dir_model, "county_week_panel_2016_2025.rds")
)

panel_week_comparison <- full_join(
  existing_panel %>%
    filter(iso_year %in% 2017:2024, lab_available) %>%
    transmute(
      county,
      iso_year,
      iso_week,
      in_existing = TRUE,
      existing_tests = tests,
      existing_positive = positive,
      existing_A = fluA,
      existing_B = fluB
    ),
  county_week_subtypes %>%
    transmute(
      county,
      iso_year,
      iso_week,
      in_updated = TRUE,
      updated_tests = tests,
      updated_positive = positive,
      updated_A = influenza_A,
      updated_B = influenza_B
    ),
  by = c("county", "iso_year", "iso_week")
) %>%
  mutate(
    in_existing = coalesce(in_existing, FALSE),
    in_updated = coalesce(in_updated, FALSE),
    week_presence = case_when(
      in_existing & in_updated ~ "both",
      in_existing ~ "existing_only",
      TRUE ~ "updated_only"
    )
  )

panel_impact_summary <- panel_week_comparison %>%
  summarise(
    existing_county_weeks = sum(in_existing),
    updated_county_weeks = sum(in_updated),
    shared_county_weeks = sum(week_presence == "both"),
    existing_only_county_weeks = sum(week_presence == "existing_only"),
    updated_only_county_weeks = sum(week_presence == "updated_only"),
    weeks_with_test_count_change = sum(
      week_presence == "both" & existing_tests != updated_tests,
      na.rm = TRUE
    ),
    weeks_with_positive_count_change = sum(
      week_presence == "both" & existing_positive != updated_positive,
      na.rm = TRUE
    ),
    weeks_with_A_count_change = sum(
      week_presence == "both" & existing_A != updated_A,
      na.rm = TRUE
    ),
    weeks_with_B_count_change = sum(
      week_presence == "both" & existing_B != updated_B,
      na.rm = TRUE
    ),
    existing_tests = sum(existing_tests, na.rm = TRUE),
    updated_tests = sum(updated_tests, na.rm = TRUE),
    existing_positive = sum(existing_positive, na.rm = TRUE),
    updated_positive = sum(updated_positive, na.rm = TRUE),
    existing_A = sum(existing_A, na.rm = TRUE),
    updated_A = sum(updated_A, na.rm = TRUE),
    existing_B = sum(existing_B, na.rm = TRUE),
    updated_B = sum(updated_B, na.rm = TRUE)
  )

analytic_impact_summary <- existing_panel %>%
  filter(
    iso_year %in% 2017:2024,
    lab_available,
    met_valid_5day_main,
    !is.na(temperature),
    !is.na(altitude),
    !is.na(altitude_1000_c)
  ) %>%
  select(
    county,
    iso_year,
    iso_week,
    existing_tests = tests,
    existing_positive = positive,
    existing_A = fluA,
    existing_B = fluB
  ) %>%
  left_join(
    county_week_subtypes %>%
      select(
        county,
        iso_year,
        iso_week,
        updated_tests = tests,
        updated_positive = positive,
        updated_A = influenza_A,
        updated_B = influenza_B
      ),
    by = c("county", "iso_year", "iso_week")
  ) %>%
  summarise(
    analytic_county_weeks = n(),
    counties = n_distinct(county),
    weeks_with_test_count_change = sum(existing_tests != updated_tests, na.rm = TRUE),
    weeks_with_positive_count_change = sum(
      existing_positive != updated_positive,
      na.rm = TRUE
    ),
    weeks_with_A_count_change = sum(existing_A != updated_A, na.rm = TRUE),
    weeks_with_B_count_change = sum(existing_B != updated_B, na.rm = TRUE),
    existing_tests = sum(existing_tests),
    updated_tests = sum(updated_tests, na.rm = TRUE),
    existing_positive = sum(existing_positive),
    updated_positive = sum(updated_positive, na.rm = TRUE),
    existing_A = sum(existing_A),
    updated_A = sum(updated_A, na.rm = TRUE),
    existing_B = sum(existing_B),
    updated_B = sum(updated_B, na.rm = TRUE)
  )

write_csv_utf8(audit_summary, file.path(out_dir, "updated_lab_audit_summary.csv"))
write_csv_utf8(subtype_distribution, file.path(out_dir, "subtype_by_positivity_distribution.csv"))
write_csv_utf8(pcr_distribution, file.path(out_dir, "raw_pcr_subtype_distribution.csv"))
write_csv_utf8(old_new_subtype_crosswalk, file.path(out_dir, "old_new_subtype_crosswalk.csv"))
write_csv_utf8(version_differences, file.path(out_dir, "version_differences.csv"))
write_csv_utf8(
  reported_week_mismatch_examples,
  file.path(out_dir, "reported_week_mismatch_examples.csv")
)
write_csv_utf8(yearly_coverage, file.path(out_dir, "yearly_subtype_coverage.csv"))
write_csv_utf8(county_coverage, file.path(out_dir, "county_subtype_coverage_2017_2024.csv"))
write_csv_utf8(county_week_subtypes, file.path(out_dir, "county_week_subtypes_2017_2024.csv"))
write_csv_utf8(weekly_feasibility, file.path(out_dir, "subtype_model_feasibility.csv"))
write_csv_utf8(
  panel_week_comparison,
  file.path(out_dir, "existing_vs_updated_county_week_comparison.csv")
)
write_csv_utf8(
  panel_impact_summary,
  file.path(out_dir, "existing_panel_impact_summary.csv")
)
write_csv_utf8(
  analytic_impact_summary,
  file.path(out_dir, "existing_analytic_sample_impact_summary.csv")
)

report <- c(
  "# Updated Laboratory Subtype Audit",
  "",
  paste0("Updated workbook: `", basename(updated_lab_path), "`"),
  paste0("Previous workbook: `0718开始分析/", basename(previous_lab_path), "`"),
  "",
  "## Scope",
  "",
  "- Compare the revised laboratory workbook with the version used by the existing pipeline.",
  "- Validate the new influenza A/B code against positivity and the previous text result.",
  "- Quantify subtype coverage by year, county and county-week for 2017-2024.",
  "- Do not overwrite the previous county-week panel or fitted models during this audit.",
  "",
  "## Required interpretation boundary",
  "",
  "- Influenza A and B outcomes should use all laboratory tests as denominators.",
  "- A-versus-B composition among positives is a separate secondary analysis.",
  "- Missing subtype among positive records must not be recoded as influenza A or B.",
  "- The updated source should not replace the existing analytic input until version consistency checks pass.",
  "",
  "See the CSV files in this directory for numerical results."
)
write_lines_utf8(report, file.path(out_dir, "updated_lab_subtype_audit_report.md"))

message("Updated laboratory subtype audit complete: ", out_dir)
