# Reconstructed template for daily ChinaMet grid -> administrative-area means.
# This is NOT the recovered historical extraction script. The manuscript
# authors report using daily polygon means and Tianditu 2024 boundaries
# (GS(2024)0650); the author confirms the cell-center inclusion rule. The exact
# original code and NoData handling were not found. This implementation uses
# an unweighted mean of raster cells whose centers fall inside each polygon.
#
# Required environment variables:
#   QH_GRID_MANIFEST       CSV: variable,date,raster_path[,layer]
#   QH_BOUNDARY_SHP        path to a polygon shapefile
#   QH_ADMIN_ID_FIELD      field identifying the chosen admin unit
# Optional:
#   QH_ADMIN_NAME_FIELD    field with human-readable unit name
#   QH_ADMIN_LEVEL         label written to output (default: county)
#   QH_BOUNDARY_FILTER_FIELD + QH_BOUNDARY_FILTER_VALUE
#                          exact-match filter, e.g. to retain study-area rows
#   QH_GRID_MEAN_OUTPUT    output CSV (default: admin_daily_means.csv)
#
# Run from the candidate-code root. If the boundary file covers all China,
# provide only the intended study-area polygons or configure the exact filter.
# For prefecture/province means, rerun with the corresponding ID/name field;
# dissolved geometry is created from repeated IDs if needed.

required <- c("terra", "sf", "dplyr", "readr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required R packages: ", paste(missing, collapse = ", "))
}

get_required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

manifest_path <- normalizePath(get_required_env("QH_GRID_MANIFEST"),
                               winslash = "/", mustWork = TRUE)
boundary_path <- normalizePath(get_required_env("QH_BOUNDARY_SHP"),
                               winslash = "/", mustWork = TRUE)
id_field <- get_required_env("QH_ADMIN_ID_FIELD")
name_field <- Sys.getenv("QH_ADMIN_NAME_FIELD", unset = "")
admin_level <- Sys.getenv("QH_ADMIN_LEVEL", unset = "county")
output_path <- Sys.getenv("QH_GRID_MEAN_OUTPUT", unset = "admin_daily_means.csv")
filter_field <- Sys.getenv("QH_BOUNDARY_FILTER_FIELD", unset = "")
filter_value <- Sys.getenv("QH_BOUNDARY_FILTER_VALUE", unset = "")

manifest <- readr::read_csv(manifest_path, show_col_types = FALSE,
                            col_types = readr::cols(.default = readr::col_character()))
needed <- c("variable", "date", "raster_path")
if (!all(needed %in% names(manifest))) {
  stop("Manifest must contain columns: ", paste(needed, collapse = ", "))
}
manifest$date <- as.Date(manifest$date)
if (anyNA(manifest$date)) stop("Manifest date values must parse as YYYY-MM-DD.")
if (!"layer" %in% names(manifest)) manifest$layer <- ""
manifest$layer[is.na(manifest$layer)] <- ""
if (anyNA(manifest$variable) || any(!nzchar(manifest$variable))) {
  stop("Manifest variable values must be non-empty.")
}
if (anyNA(manifest$raster_path) || any(!nzchar(manifest$raster_path))) {
  stop("Manifest raster_path values must be non-empty.")
}
if (anyDuplicated(manifest[c("variable", "date")])) {
  stop("Manifest must have at most one raster per variable and date; use layer to select a band.")
}

manifest_dir <- dirname(manifest_path)
is_absolute_path <- function(x) grepl("^([A-Za-z]:[/\\\\]|[/\\\\]{2})", x)
manifest$raster_path <- vapply(manifest$raster_path, function(p) {
  if (!is_absolute_path(p)) p <- file.path(manifest_dir, p)
  normalizePath(p, winslash = "/", mustWork = TRUE)
}, character(1))

admin <- sf::st_read(boundary_path, quiet = TRUE, stringsAsFactors = FALSE)
if (is.na(sf::st_crs(admin))) stop("Boundary file has no CRS; do not assign one by guesswork.")
if (!id_field %in% names(admin)) {
  stop("ID field not found. Available fields: ", paste(names(admin), collapse = ", "))
}
group_fields <- id_field
if (nzchar(name_field)) {
  if (!name_field %in% names(admin)) stop("Name field not found: ", name_field)
  group_fields <- c(group_fields, name_field)
}
if (nzchar(filter_field) || nzchar(filter_value)) {
  if (!nzchar(filter_field) || !nzchar(filter_value)) {
    stop("Set both QH_BOUNDARY_FILTER_FIELD and QH_BOUNDARY_FILTER_VALUE, or neither.")
  }
  if (!filter_field %in% names(admin)) stop("Filter field not found: ", filter_field)
  admin <- admin[as.character(admin[[filter_field]]) == filter_value, , drop = FALSE]
  if (!nrow(admin)) stop("Boundary filter retained zero features.")
}

admin <- sf::st_make_valid(admin)
admin <- dplyr::group_by(admin, dplyr::across(dplyr::all_of(group_fields)))
admin <- dplyr::summarise(admin, .groups = "drop")
if (!nrow(admin)) stop("No administrative polygons remain after filtering/dissolving.")
admin_key <- data.frame(
  polygon_index = seq_len(nrow(admin)),
  admin_id = as.character(admin[[id_field]]),
  stringsAsFactors = FALSE
)
if (nzchar(name_field)) admin_key$admin_name <- as.character(admin[[name_field]])

results <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  raster <- terra::rast(row$raster_path)
  if (nzchar(row$layer)) {
    layer_index <- suppressWarnings(as.integer(row$layer))
    layer <- if (!is.na(layer_index)) layer_index else row$layer
    raster <- raster[[layer]]
  } else if (terra::nlyr(raster) != 1L) {
    stop("Raster has multiple layers; set the manifest 'layer' column: ", row$raster_path)
  }
  if (terra::nlyr(raster) != 1L) stop("Exactly one layer must be selected per manifest row.")
  raster_crs <- terra::crs(raster)
  if (!nzchar(raster_crs)) stop("Raster CRS is missing: ", row$raster_path)

  admin_projected <- sf::st_transform(admin, raster_crs)
  # touches=FALSE + small=FALSE means only raster cells whose centers lie
  # inside each polygon are included; no cell-area weighting is applied.
  extracted <- terra::extract(
    raster, terra::vect(admin_projected), ID = TRUE,
    touches = FALSE, small = FALSE
  )
  if (ncol(extracted) < 2L) stop("No raster value column extracted: ", row$raster_path)
  value_name <- names(extracted)[2L]
  values <- data.frame(
    polygon_index = extracted$ID,
    grid_value = as.numeric(extracted[[value_name]])
  )
  joined <- merge(admin_key, values, by = "polygon_index", all.x = TRUE, sort = FALSE)
  group_cols <- c("admin_id", if ("admin_name" %in% names(joined)) "admin_name")
  summary <- dplyr::summarise(
    dplyr::group_by(joined, dplyr::across(dplyr::all_of(group_cols))),
    mean_value = if (sum(!is.na(grid_value)) > 0L) mean(grid_value, na.rm = TRUE) else NA_real_,
    n_valid_cells = sum(!is.na(grid_value)),
    .groups = "drop"
  )
  summary$admin_level <- admin_level
  summary$date <- row$date
  summary$variable <- row$variable
  summary$mean_method <- "unweighted mean; cell centers within polygon"
  summary$raster_path <- row$raster_path
  summary$boundary_path <- boundary_path
  results[[i]] <- summary
  message(sprintf("Processed %s %s (%d/%d)", row$variable, row$date, i, nrow(manifest)))
}

output <- dplyr::bind_rows(results)
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(output, output_path, na = "")
writeLines(capture.output(utils::sessionInfo()), paste0(output_path, ".sessionInfo.txt"))
message("Wrote ", normalizePath(output_path, winslash = "/", mustWork = TRUE))
