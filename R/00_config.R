suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(splines)
  library(dlnm)
  library(broom)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- list.dirs(project_root, recursive = FALSE, full.names = TRUE)
data_dir <- data_dir[grepl("0718", basename(data_dir))][1]
if (is.na(data_dir) || !dir.exists(data_dir)) {
  stop("Cannot find the data directory whose name starts with 0718.")
}

dir_audit <- file.path(project_root, "01_data_audit")
dir_clean <- file.path(project_root, "02_data_cleaning")
dir_panel <- file.path(project_root, "03_panel_construction")
dir_desc <- file.path(project_root, "04_descriptive")
dir_dlnm <- file.path(project_root, "05_frequentist_dlnm")
dir_fig <- file.path(project_root, "08_figures")
dir_tab <- file.path(project_root, "09_tables")
dir_model <- file.path(project_root, "model_outputs")

for (d in c(dir_audit, dir_clean, dir_panel, dir_desc, dir_dlnm,
            dir_fig, dir_tab, dir_model)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

find_file <- function(pattern) {
  files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
  hit <- files[grepl(pattern, basename(files))]
  if (length(hit) == 0) stop("Cannot find file matching: ", pattern)
  hit[1]
}

file_lab <- find_file("核酸")
file_ili <- find_file("ILI")
file_met <- find_file("气象")
file_alt <- find_file("海拔")
file_gdp <- find_file("GDP")
file_pop_density <- find_file("人口密度")
file_age <- find_file("年龄结构")

as_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

county_alias <- c(
  "大通县" = "大通回族土族自治县",
  "乐都县" = "乐都区",
  "门源县" = "门源回族自治县",
  "同仁县" = "同仁市",
  "河南县" = "河南蒙古族自治县",
  "大柴旦行委" = "大柴旦行政委员会"
)

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    return(as.Date(x, origin = "1899-12-30"))
  }
  xx <- str_trim(as.character(x))
  xx[xx %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  numeric_value <- suppressWarnings(as.numeric(xx))
  is_excel_serial <- !is.na(numeric_value) &
    str_detect(xx, "^[0-9]+(?:\\.[0-9]+)?$") &
    numeric_value >= 1 &
    numeric_value <= 100000
  result <- as.Date(rep(NA_character_, length(xx)))
  result[is_excel_serial] <- as.Date(
    numeric_value[is_excel_serial],
    origin = "1899-12-30"
  )
  needs_parsing <- !is_excel_serial & !is.na(xx)
  if (any(needs_parsing)) {
    parsed <- suppressWarnings(parse_date_time(
      xx[needs_parsing],
      orders = c(
        "ymd HMS",
        "ymd HM",
        "ymd",
        "Ymd HMS",
        "Ymd",
        "mdy",
        "dmy"
      ),
      tz = "Asia/Shanghai"
    ))
    result[needs_parsing] <- as.Date(parsed)
  }
  result
}

normalize_county <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace(x, "^[]】）)]+", "")
  x <- dplyr::recode(x, !!!county_alias, .default = x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

add_iso_week <- function(df, date_col = "date") {
  dt <- df[[date_col]]
  df %>%
    mutate(
      iso_year = isoyear(dt),
      iso_week = isoweek(dt),
      week_start = floor_date(dt, unit = "week", week_start = 1),
      week_id = sprintf("%04d-W%02d", iso_year, iso_week)
    )
}

write_csv_utf8 <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

write_lines_utf8 <- function(x, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(x, con = con, useBytes = TRUE)
}

source_flag <- function(x) {
  ifelse(is.na(x), "missing", "observed")
}
