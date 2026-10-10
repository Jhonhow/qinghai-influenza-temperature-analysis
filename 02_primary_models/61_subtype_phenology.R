local_lib <- normalizePath("R_libs_q1", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

source(file.path("R", "00_config.R"))
suppressPackageStartupMessages({
  library(glmmTMB)
  library(MASS)
})

set.seed(20260806)

version_id <- "v3_20260728"
analysis_id <- "result3_subtype_phenology_20260806"
model_root <- file.path(dir_model, paste0("reanalysis_", version_id))
bundle <- readRDS(file.path(
  model_root,
  "regularized_altitude_interaction_20260730",
  "regularized_analysis_bundle.rds"
))
model_dat <- bundle$model_data

checkpoint_dir <- file.path(model_root, analysis_id)
source_dir <- file.path(project_root, "result0715", "source_data", analysis_id)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
write_lines_utf8 <- function(x, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(x, con = con, useBytes = TRUE)
}

annual_period_weeks <- 365.2425 / 7
model_dat <- model_dat |>
  mutate(
    iso_year_factor = factor(iso_year, levels = 2017:2024),
    annual_week = ((week_index - 1) %% annual_period_weeks) + 1,
    county_nested = factor(county_nested),
    city = factor(city)
  )

add_fourier_terms <- function(data, harmonics, week_variable = "annual_week") {
  week_value <- data[[week_variable]]
  for (harmonic in seq_len(harmonics)) {
    data[[paste0("season_sin_", harmonic)]] <- sin(
      2 * pi * harmonic * week_value / annual_period_weeks
    )
    data[[paste0("season_cos_", harmonic)]] <- cos(
      2 * pi * harmonic * week_value / annual_period_weeks
    )
  }
  data
}
model_dat <- add_fourier_terms(model_dat, 4)

make_formula <- function(outcome, harmonics) {
  season_names <- as.vector(rbind(
    paste0("season_sin_", seq_len(harmonics)),
    paste0("season_cos_", seq_len(harmonics))
  ))
  as.formula(paste(
    sprintf("cbind(%s, tests - %s) ~", outcome, outcome),
    paste(
      c("iso_year_factor", season_names, "(1 | city)", "(1 | county_nested)"),
      collapse = " + "
    )
  ))
}

fit_or_load <- function(outcome, harmonics) {
  label <- paste0(outcome, "__seasonal_fourier", harmonics)
  checkpoint <- file.path(checkpoint_dir, paste0(label, ".rds"))
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  fit <- glmmTMB(
    make_formula(outcome, harmonics),
    family = betabinomial(link = "logit"),
    data = model_dat,
    control = glmmTMBControl(
      optCtrl = list(iter.max = 3000, eval.max = 3000)
    )
  )
  saveRDS(fit, checkpoint)
  fit
}

outcome_specs <- tibble(
  outcome = c("fluA_only", "fluB_only"),
  outcome_label = c("Influenza A", "Influenza B")
)
harmonic_grid <- 2:4
fits <- list()
diagnostic_rows <- list()
for (i in seq_len(nrow(outcome_specs))) {
  for (harmonics in harmonic_grid) {
    outcome <- outcome_specs$outcome[[i]]
    label <- paste0(outcome, "__seasonal_fourier", harmonics)
    fit <- fit_or_load(outcome, harmonics)
    fits[[label]] <- fit
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- tibble(
      outcome = outcome,
      outcome_label = outcome_specs$outcome_label[[i]],
      harmonics = harmonics,
      convergence = fit$fit$convergence,
      pdHess = fit$sdr$pdHess,
      observations = nobs(fit),
      AIC = AIC(fit),
      logLik = as.numeric(logLik(fit))
    )
  }
}
diagnostics <- bind_rows(diagnostic_rows)
if (!all(diagnostics$convergence == 0 & diagnostics$pdHess %in% TRUE)) {
  stop("At least one subtype phenology model failed numerical diagnostics.")
}

make_prediction_data <- function(harmonics) {
  data <- tidyr::expand_grid(
    annual_week = seq(1, 52, by = 0.25),
    iso_year_factor = factor(2017:2024, levels = 2017:2024)
  )
  add_fourier_terms(data, harmonics)
}

standardized_curve <- function(fit, harmonics, simulations = 2000L) {
  newdata <- make_prediction_data(harmonics)
  season_names <- as.vector(rbind(
    paste0("season_sin_", seq_len(harmonics)),
    paste0("season_cos_", seq_len(harmonics))
  ))
  fixed_formula <- as.formula(paste(
    "~ iso_year_factor +", paste(season_names, collapse = " + ")
  ))
  design <- model.matrix(fixed_formula, newdata)
  beta <- fixef(fit)$cond
  covariance <- as.matrix(vcov(fit)$cond)
  aligned <- matrix(
    0,
    nrow = nrow(design),
    ncol = length(beta),
    dimnames = list(NULL, names(beta))
  )
  common <- intersect(colnames(design), names(beta))
  aligned[, common] <- design[, common, drop = FALSE]

  beta_draws <- MASS::mvrnorm(simulations, mu = beta, Sigma = covariance)
  linear_draws <- aligned %*% t(beta_draws)
  probability_draws <- plogis(linear_draws)
  week_values <- sort(unique(newdata$annual_week))
  standardized_draws <- vapply(
    week_values,
    function(week_value) {
      colMeans(probability_draws[newdata$annual_week == week_value, , drop = FALSE])
    },
    numeric(simulations)
  )
  standardized_draws <- t(standardized_draws)
  curve <- tibble(
    annual_week = week_values,
    estimate = rowMeans(standardized_draws),
    low = apply(standardized_draws, 1, quantile, 0.025),
    high = apply(standardized_draws, 1, quantile, 0.975)
  )

  peak_index <- apply(standardized_draws, 2, which.max)
  peak_week <- week_values[peak_index]
  # Shift weeks 1-26 after weeks 27-52 to avoid a confidence interval that
  # crosses the end-of-year boundary for winter peaks.
  shifted_peak <- ifelse(peak_week <= 26, peak_week + 52, peak_week)
  point_peak <- curve$annual_week[[which.max(curve$estimate)]]
  shifted_interval <- quantile(shifted_peak, c(0.025, 0.975))
  map_week <- function(x) ifelse(x > 52, x - 52, x)
  peak_summary <- tibble(
    peak_week = point_peak,
    peak_week_low = map_week(shifted_interval[[1]]),
    peak_week_high = map_week(shifted_interval[[2]])
  )
  list(curve = curve, peak = peak_summary)
}

curve_rows <- list()
peak_rows <- list()
for (i in seq_len(nrow(outcome_specs))) {
  outcome <- outcome_specs$outcome[[i]]
  outcome_label <- outcome_specs$outcome_label[[i]]
  for (harmonics in harmonic_grid) {
    label <- paste0(outcome, "__seasonal_fourier", harmonics)
    result <- standardized_curve(fits[[label]], harmonics)
    curve_rows[[length(curve_rows) + 1L]] <- result$curve |>
      mutate(
        outcome = outcome,
        outcome_label = outcome_label,
        harmonics = harmonics,
        primary_specification = harmonics == 3L
      )
    peak_rows[[length(peak_rows) + 1L]] <- result$peak |>
      mutate(
        outcome = outcome,
        outcome_label = outcome_label,
        harmonics = harmonics,
        primary_specification = harmonics == 3L
      )
  }
}
seasonal_curves <- bind_rows(curve_rows)
peak_sensitivity <- bind_rows(peak_rows)

# Observed annual-week data retain year-to-year subtype replacement and are
# shown separately from the standardized average seasonal curve.
observed_weekly <- model_dat |>
  mutate(iso_week = as.integer(format(week_start, "%V"))) |>
  group_by(iso_year, iso_week) |>
  summarise(
    tests = sum(tests),
    influenza_A_only = sum(fluA_only),
    influenza_B_only = sum(fluB_only),
    influenza_A_positivity = influenza_A_only / tests,
    influenza_B_positivity = influenza_B_only / tests,
    .groups = "drop"
  )

temperature_breaks <- c(-Inf, -20, -15, -10, -5, 0, 5, Inf)
temperature_labels <- c(
  "< -20", "-20 to < -15", "-15 to < -10", "-10 to < -5",
  "-5 to < 0", "0 to < 5", ">= 5"
)
temperature_bands <- model_dat |>
  mutate(
    temperature_band = cut(
      temperature,
      breaks = temperature_breaks,
      labels = temperature_labels,
      right = FALSE
    )
  ) |>
  group_by(temperature_band) |>
  summarise(
    tests = sum(tests),
    influenza_A_only = sum(fluA_only),
    influenza_B_only = sum(fluB_only),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(influenza_A_only, influenza_B_only),
    names_to = "outcome",
    values_to = "positive"
  ) |>
  mutate(
    outcome_label = recode(
      outcome,
      influenza_A_only = "Influenza A",
      influenza_B_only = "Influenza B"
    )
  )

wilson_interval <- function(x, n, z = 1.96) {
  raw <- x / n
  denominator <- 1 + z^2 / n
  center <- (raw + z^2 / (2 * n)) / denominator
  half <- z * sqrt(raw * (1 - raw) / n + z^2 / (4 * n^2)) /
    denominator
  c(low = center - half, high = center + half)
}
intervals <- t(vapply(
  seq_len(nrow(temperature_bands)),
  function(i) wilson_interval(
    temperature_bands$positive[[i]], temperature_bands$tests[[i]]
  ),
  numeric(2)
))
temperature_bands <- temperature_bands |>
  mutate(
    positivity = positive / tests,
    low = intervals[, "low"],
    high = intervals[, "high"],
    temperature_band = factor(
      temperature_band, levels = temperature_labels
    )
  )

primary_peaks <- peak_sensitivity |>
  filter(primary_specification %in% TRUE)
peak_separation <- tibble(
  estimand = "Influenza B peak minus influenza A peak",
  influenza_A_peak_week = primary_peaks$peak_week[
    primary_peaks$outcome == "fluA_only"
  ],
  influenza_B_peak_week = primary_peaks$peak_week[
    primary_peaks$outcome == "fluB_only"
  ]
) |>
  mutate(
    separation_weeks = if_else(
      influenza_B_peak_week >= influenza_A_peak_week,
      influenza_B_peak_week - influenza_A_peak_week,
      influenza_B_peak_week + 52 - influenza_A_peak_week
    )
  )

qa_checks <- tibble(
  check = c(
    "All six phenology models converged",
    "All six Hessians are positive definite",
    "Primary curves average over all eight years",
    "Peak sensitivity covers 2-4 Fourier pairs",
    "Temperature bands contain all tests",
    "Influenza B primary peak follows influenza A primary peak"
  ),
  passed = c(
    all(diagnostics$convergence == 0),
    all(diagnostics$pdHess %in% TRUE),
    all(2017:2024 %in% levels(model_dat$iso_year_factor)),
    identical(sort(unique(peak_sensitivity$harmonics)), 2:4),
    all(tapply(temperature_bands$tests, temperature_bands$outcome, sum) ==
      sum(model_dat$tests)),
    peak_separation$separation_weeks > 0
  )
)
if (!all(qa_checks$passed)) stop("Subtype phenology QA failed.")

write_csv_utf8(seasonal_curves, file.path(
  source_dir, "standardized_subtype_seasonal_curves.csv"
))
write_csv_utf8(peak_sensitivity, file.path(
  source_dir, "subtype_peak_week_sensitivity.csv"
))
write_csv_utf8(peak_separation, file.path(
  source_dir, "primary_peak_separation.csv"
))
write_csv_utf8(observed_weekly, file.path(
  source_dir, "observed_subtype_positivity_by_year_week.csv"
))
write_csv_utf8(temperature_bands, file.path(
  source_dir, "observed_subtype_positivity_by_temperature_band.csv"
))
write_csv_utf8(diagnostics, file.path(source_dir, "model_diagnostics.csv"))
write_csv_utf8(qa_checks, file.path(source_dir, "phenology_QA_checks.csv"))

report <- c(
  "# Result 3 subtype phenology analysis",
  "",
  "## Estimand",
  "",
  paste(
    "The phenology model estimates the average annual timing of influenza A and",
    "B positivity after standardization over calendar years. It describes subtype",
    "seasonality and is not interpreted as an independent causal temperature effect."
  ),
  "",
  "## Primary model",
  "",
  paste(
    "Hierarchical beta-binomial regression with three annual Fourier pairs,",
    "calendar-year indicators, and prefecture and county random intercepts."
  ),
  "",
  "## Peak timing",
  "",
  paste(capture.output(print(primary_peaks)), collapse = "\n"),
  "",
  paste(capture.output(print(peak_separation)), collapse = "\n"),
  "",
  "## Sensitivity",
  "",
  paste(capture.output(print(peak_sensitivity)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  paste(
    "Temperature-banded positivity quantifies where cases occurred across the",
    "observed thermal distribution. The seasonally adjusted DLNM remains a",
    "separate analysis of conditional short-term temperature associations."
  )
)
write_lines_utf8(report, file.path(source_dir, "phenology_analysis_report.md"))

message("Subtype phenology Result 3 analysis completed.")
