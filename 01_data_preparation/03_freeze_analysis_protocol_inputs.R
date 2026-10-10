source(file.path("R", "00_config.R"))

freeze_dir <- file.path(dir_audit, "analysis_protocol_freeze_20260728")
if (!dir.exists(freeze_dir)) dir.create(freeze_dir, recursive = TRUE)

root_lab <- file.path(
  project_root,
  "2016-2025流感核酸监测信息20260703_已匹配区县.xlsx"
)
legacy_lab <- file.path(
  data_dir,
  "2016-2025流感核酸监测信息20260703_已匹配区县.xlsx"
)
root_ili <- file.path(
  project_root,
  "2016-2025青海省ILI20260703_已匹配区县.xlsx"
)
legacy_ili <- file.path(
  data_dir,
  "2016-2025青海省ILI20260703_已匹配区县.xlsx"
)

source_files <- tibble::tribble(
  ~role, ~path, ~planned_use,
  "laboratory surveillance", root_lab, "canonical primary outcome input",
  "laboratory surveillance (previous)", legacy_lab, "version comparison only",
  "ILI surveillance", root_ili, "candidate canonical secondary outcome input",
  "ILI surveillance (previous)", legacy_ili, "version comparison only",
  "meteorology", file.path(data_dir, "青海省16-24年气象数据.xlsx"), "canonical exposure and covariate input",
  "altitude", file.path(data_dir, "各区县对应的海拔.xlsx"), "canonical continuous effect-modifier input",
  "GDP", file.path(data_dir, "青海各区县16-24GDP数据.xlsx"), "descriptive and sensitivity input",
  "population density", file.path(data_dir, "青海各区县人口密度数据.xlsx"), "descriptive and sensitivity input",
  "age structure", file.path(data_dir, "青海省各区县年龄结构数据.xlsx"), "descriptive and sensitivity input"
)

md5_file <- function(path) {
  unname(tools::md5sum(path))
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

file_manifest <- source_files %>%
  mutate(
    exists = file.exists(path),
    file_name = basename(path),
    directory = dirname(path),
    size_bytes = if_else(exists, file.size(path), NA_real_),
    modified_time = if_else(
      exists,
      format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S"),
      NA_character_
    ),
    md5 = if_else(exists, unname(tools::md5sum(path)), NA_character_),
    sha256 = if_else(
      exists,
      vapply(path, sha256_file, character(1)),
      NA_character_
    )
  ) %>%
  select(
    role,
    planned_use,
    file_name,
    directory,
    size_bytes,
    modified_time,
    md5,
    sha256,
    exists
  )

write_csv_utf8(
  file_manifest,
  file.path(freeze_dir, "input_file_manifest.csv")
)

inventory_one_workbook <- function(path, role) {
  if (!file.exists(path)) {
    return(tibble(
      role = role,
      file_name = basename(path),
      sheet = NA_character_,
      n_rows = NA_integer_,
      n_cols = NA_integer_,
      columns = NA_character_,
      read_error = "file not found"
    ))
  }

  purrr::map_dfr(excel_sheets(path), function(sheet_name) {
    result <- tryCatch(
      {
        dat <- read_excel(path, sheet = sheet_name)
        tibble(
          role = role,
          file_name = basename(path),
          sheet = sheet_name,
          n_rows = nrow(dat),
          n_cols = ncol(dat),
          columns = paste(names(dat), collapse = " | "),
          read_error = NA_character_
        )
      },
      error = function(e) {
        tibble(
          role = role,
          file_name = basename(path),
          sheet = sheet_name,
          n_rows = NA_integer_,
          n_cols = NA_integer_,
          columns = NA_character_,
          read_error = conditionMessage(e)
        )
      }
    )
    result
  })
}

workbook_inventory <- purrr::map2_dfr(
  source_files$path,
  source_files$role,
  inventory_one_workbook
)

write_csv_utf8(
  workbook_inventory,
  file.path(freeze_dir, "workbook_sheet_inventory.csv")
)

compare_workbook_sheets <- function(current_path, previous_path, source_label) {
  current_sheets <- excel_sheets(current_path)
  previous_sheets <- excel_sheets(previous_path)
  common_sheets <- intersect(current_sheets, previous_sheets)

  purrr::map_dfr(common_sheets, function(sheet_name) {
    current <- read_excel(current_path, sheet = sheet_name)
    previous <- read_excel(previous_path, sheet = sheet_name)
    common_columns <- intersect(names(current), names(previous))
    same_dimensions <- identical(dim(current), dim(previous))
    same_column_names <- identical(names(current), names(previous))
    same_common_values <- identical(
      current[, common_columns, drop = FALSE],
      previous[, common_columns, drop = FALSE]
    )

    tibble(
      source = source_label,
      sheet = sheet_name,
      current_rows = nrow(current),
      previous_rows = nrow(previous),
      current_columns = ncol(current),
      previous_columns = ncol(previous),
      common_columns = length(common_columns),
      same_dimensions = same_dimensions,
      same_column_names = same_column_names,
      same_common_values = same_common_values
    )
  })
}

version_comparison <- bind_rows(
  compare_workbook_sheets(root_lab, legacy_lab, "laboratory"),
  compare_workbook_sheets(root_ili, legacy_ili, "ILI")
)

write_csv_utf8(
  version_comparison,
  file.path(freeze_dir, "duplicate_workbook_comparison.csv")
)

cat("Protocol input audit completed:", freeze_dir, "\n")
