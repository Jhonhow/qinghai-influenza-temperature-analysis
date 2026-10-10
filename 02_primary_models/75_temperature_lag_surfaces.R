local_lib <- normalizePath("R_libs_q1", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

source(file.path("R", "00_config.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(dlnm)
  library(glmmTMB)
})

set.seed(20260814)

version_id <- "v3_20260728"
analysis_id <- "result3_harmonized_portable_temperature_lag_surfaces_20260814"
model_root <- file.path(dir_model, paste0("reanalysis_", version_id))
bundle_dir <- file.path(
  model_root,
  "regularized_altitude_interaction_20260730"
)
primary_fit_dir <- file.path(
  model_root,
  "portable_calendar_full_fit_20260730"
)
subtype_fit_dir <- file.path(
  model_root,
  "figure4_harmonized_portable_subtypes_20260809"
)
source_dir <- file.path(
  project_root,
  "result0715",
  "source_data",
  analysis_id
)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
write_lines_utf8 <- function(x, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(x, con = con, useBytes = TRUE)
}

bundle <- readRDS(file.path(bundle_dir, "regularized_analysis_bundle.rds"))
model_dat <- bundle$model_data
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
raw_int_names <- bundle$raw_int_names
z_int_names <- bundle$z_int_names

# These are the frozen portable-calendar fits used by the primary overall
# model and the harmonized subtype analysis. No model is selected or refitted
# according to the displayed surfaces.
fits <- list(
  `Overall positivity` = readRDS(
    file.path(primary_fit_dir, "full_data_map.rds")
  ),
  `Influenza A-only positivity` = readRDS(
    file.path(subtype_fit_dir, "fluA_only__portable_map.rds")
  ),
  `Influenza B-only positivity` = readRDS(
    file.path(subtype_fit_dir, "fluB_only__portable_map.rds")
  )
)

interaction_center <- vapply(
  model_dat[raw_int_names], mean, numeric(1), na.rm = TRUE
)
interaction_scale <- vapply(
  model_dat[raw_int_names], stats::sd, numeric(1), na.rm = TRUE
)
if (any(!is.finite(interaction_scale) | interaction_scale <= 0)) {
  stop("Invalid frozen altitude-interaction scaling.")
}

figure2_curve_path <- file.path(
  project_root,
  "result0715",
  "source_data",
  "figure2_exposure_distribution_supported_contrasts_20260805",
  "panel_b_mixed_reference_curves.csv"
)
figure4_curve_path <- file.path(
  project_root,
  "result0715",
  "source_data",
  "figure4_harmonized_portable_subtypes_20260809",
  "figure4_harmonized_curves.csv"
)

figure2_curves <- read.csv(
  figure2_curve_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
figure4_curves <- read.csv(
  figure4_curve_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

p50_spec <- figure2_curves |>
  filter(
    altitude_anchor == "P50",
    reference_temperature == -5,
    within_p05_p95 %in% TRUE
  ) |>
  arrange(temperature)

if (nrow(p50_spec) == 0L) stop("No frozen P50 Figure 2 curve was found.")
if (n_distinct(p50_spec$altitude) != 1L ||
    n_distinct(p50_spec$altitude_1000_c) != 1L) {
  stop("The P50 Figure 2 model slice is not unique.")
}

reference_temperature <- -5
slice_altitude <- unique(p50_spec$altitude)
slice_altitude_1000_c <- unique(p50_spec$altitude_1000_c)
temperature_grid <- sort(unique(p50_spec$temperature))
lag_grid <- 0:2

extract_fixed <- function(fit) {
  list(
    beta = glmmTMB::fixef(fit)$cond,
    vcov = as.matrix(stats::vcov(fit)$cond)
  )
}

basis_array_lag <- function(temperatures, reference_temperature) {
  original_names <- colnames(cb_all)
  output <- array(
    0,
    dim = c(length(temperatures), length(lag_grid), length(original_names)),
    dimnames = list(
      as.character(temperatures),
      paste0("lag", lag_grid),
      cb_names
    )
  )
  zero_vcov <- matrix(
    0,
    nrow = length(original_names),
    ncol = length(original_names),
    dimnames = list(original_names, original_names)
  )
  for (j in seq_along(original_names)) {
    coefficient <- rep(0, length(original_names))
    names(coefficient) <- original_names
    coefficient[[j]] <- 1
    prediction <- crosspred(
      cb_all,
      coef = coefficient,
      vcov = zero_vcov,
      model.link = "logit",
      at = temperatures,
      cen = reference_temperature,
      cumul = FALSE
    )
    output[, , j] <- prediction$matfit[
      , paste0("lag", lag_grid), drop = FALSE
    ]
  }
  output
}

estimate_lag_or <- function(fit, basis_values) {
  fixed <- extract_fixed(fit)
  design <- rep(0, length(fixed$beta))
  names(design) <- names(fixed$beta)
  design[intersect(cb_names, names(design))] <-
    basis_values[intersect(cb_names, names(design))]
  interaction_values <-
    basis_values * slice_altitude_1000_c / interaction_scale
  names(interaction_values) <- z_int_names
  design[intersect(z_int_names, names(design))] <-
    interaction_values[intersect(z_int_names, names(design))]
  log_estimate <- sum(design * fixed$beta)
  standard_error <- sqrt(as.numeric(
    t(design) %*% fixed$vcov %*% design
  ))
  c(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

lag_basis <- basis_array_lag(temperature_grid, reference_temperature)
surface_rows <- list()
for (outcome_name in names(fits)) {
  fit <- fits[[outcome_name]]
  for (temperature_index in seq_along(temperature_grid)) {
    for (lag_index in seq_along(lag_grid)) {
      basis_values <- lag_basis[
        temperature_index, lag_index, , drop = TRUE
      ]
      names(basis_values) <- cb_names
      values <- estimate_lag_or(fit, basis_values)
      surface_rows[[length(surface_rows) + 1L]] <- tibble(
        outcome = outcome_name,
        temperature = temperature_grid[[temperature_index]],
        lag_week = lag_grid[[lag_index]],
        altitude_anchor = "P50",
        altitude = slice_altitude,
        altitude_1000_c = slice_altitude_1000_c,
        reference_temperature = reference_temperature,
        estimate = values[["estimate"]],
        low = values[["low"]],
        high = values[["high"]],
        log_estimate = values[["log_estimate"]],
        standard_error = values[["standard_error"]]
      )
    }
  }
}
surfaces <- bind_rows(surface_rows)

# A cumulative DLNM contrast is the sum of lag-specific log contrasts. This
# audit requires the three plotted lag surfaces to reproduce the already
# frozen cumulative curves used in Figures 2 and 4.
cumulative_from_surfaces <- surfaces |>
  group_by(outcome, temperature, altitude_anchor, altitude) |>
  summarise(
    reconstructed_cumulative_OR = exp(sum(log_estimate)),
    .groups = "drop"
  )

expected_curves <- figure4_curves |>
  filter(
    altitude_anchor == "P50",
    reference_temperature == reference_temperature,
    temperature %in% temperature_grid
  ) |>
  transmute(
    outcome = case_when(
      outcome == "Overall positivity" ~ "Overall positivity",
      outcome == "Influenza A only" ~ "Influenza A-only positivity",
      outcome == "Influenza B only" ~ "Influenza B-only positivity",
      TRUE ~ outcome
    ),
    temperature,
    frozen_cumulative_OR = estimate
  )

consistency <- cumulative_from_surfaces |>
  left_join(expected_curves, by = c("outcome", "temperature")) |>
  mutate(
    absolute_log_difference = abs(
      log(reconstructed_cumulative_OR) - log(frozen_cumulative_OR)
    )
  )

if (any(is.na(consistency$frozen_cumulative_OR))) {
  stop("At least one lag-surface point lacks a frozen cumulative comparator.")
}
max_consistency_error <- max(consistency$absolute_log_difference)
if (!is.finite(max_consistency_error) || max_consistency_error > 1e-8) {
  stop("Lag-specific surfaces do not reproduce the frozen cumulative curves.")
}

rhs_labels <- lapply(
  fits,
  function(fit) attr(terms(formula(fit)), "term.labels")
)
identical_rhs <- all(vapply(
  rhs_labels[-1], identical, logical(1), rhs_labels[[1]]
))
diagnostics <- lapply(names(fits), function(outcome_name) {
  fit <- fits[[outcome_name]]
  tibble(
    outcome = outcome_name,
    observations = nobs(fit),
    counties = n_distinct(model_dat$county),
    convergence = fit$fit$convergence,
    pdHess = fit$sdr$pdHess,
    identical_right_hand_side = identical_rhs,
    reference_temperature = reference_temperature,
    displayed_altitude = slice_altitude
  )
}) |>
  bind_rows()

qa <- tibble(
  check = c(
    "All three outcomes use 9,758 county-weeks",
    "All three outcomes use all 45 counties",
    "All three models have identical right-hand-side terms",
    "All models converged",
    "All fitted Hessians are positive definite",
    "Reference temperature is -5 degrees Celsius",
    "All surfaces use the frozen P50 altitude slice",
    "Lag-specific surfaces reproduce frozen cumulative curves"
  ),
  value = c(
    paste(diagnostics$observations, collapse = ", "),
    paste(diagnostics$counties, collapse = ", "),
    as.character(identical_rhs),
    paste(diagnostics$convergence, collapse = ", "),
    paste(diagnostics$pdHess, collapse = ", "),
    as.character(reference_temperature),
    as.character(slice_altitude),
    format(max_consistency_error, scientific = TRUE)
  ),
  passed = c(
    all(diagnostics$observations == 9758L),
    all(diagnostics$counties == 45L),
    identical_rhs,
    all(diagnostics$convergence == 0),
    all(diagnostics$pdHess),
    reference_temperature == -5,
    slice_altitude == 2980,
    max_consistency_error <= 1e-8
  )
)
if (!all(qa$passed)) stop("Harmonized Result 3 surface QA failed.")

write_csv_utf8(
  surfaces,
  file.path(source_dir, "panel_a_harmonized_temperature_lag_surfaces.csv")
)
write_csv_utf8(
  consistency,
  file.path(source_dir, "panel_a_cumulative_curve_consistency.csv")
)
write_csv_utf8(
  diagnostics,
  file.path(source_dir, "panel_a_model_diagnostics.csv")
)
write_csv_utf8(
  qa,
  file.path(source_dir, "panel_a_QA_checks.csv")
)

write_lines_utf8(
  c(
    "# Harmonized Result 3 temperature-lag surfaces",
    "",
    paste0(
      "Overall, influenza A-only and influenza B-only positivity use the same ",
      "9,758 county-weeks from 45 counties, the same hierarchical beta-binomial ",
      "DLNM, identical portable-calendar controls, continuous altitude ",
      "interaction and frozen MAP shrinkage prior."
    ),
    paste0(
      "All surfaces are evaluated at the prespecified P50 altitude slice ",
      "(", format(slice_altitude, big.mark = ","), " m) relative to -5 degrees ",
      "Celsius and display observed lags 0, 1 and 2 weeks."
    ),
    paste0(
      "The product of lag-specific odds ratios reproduces the frozen cumulative ",
      "lag 0-2 curves with maximum absolute log-scale difference ",
      format(max_consistency_error, scientific = TRUE), "."
    )
  ),
  file.path(source_dir, "panel_a_model_note.md")
)

message("Harmonized Result 3 lag surfaces completed.")
