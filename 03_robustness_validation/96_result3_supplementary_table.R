suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(dlnm)
})

set.seed(20260902)

bundle_path <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825",
  "v4_t2m_model_bundle.rds"
)
surface_path <- file.path(
  "result0715", "source_data", "v4_result3_t2m_harmonized_20260826",
  "panel_a_temperature_lag_surfaces.csv"
)
phenology_dir <- file.path(
  "result0715", "source_data", "result3_subtype_phenology_20260806"
)
output_dir <- file.path(
  "outputs", "result3_supplementary_20260902", "source"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
z_int_names <- bundle$z_int_names
interaction_scale <- bundle$interaction_scale
reference_temperature <- bundle$reference_temperature
p50 <- bundle$altitude_anchors |>
  filter(altitude_anchor == "P50")
if (nrow(p50) != 1L) stop("P50 display anchor is not unique.")
altitude <- p50$altitude[[1]]
altitude_1000_c <- p50$altitude_1000_c[[1]]

fits <- list(
  "Influenza A" = bundle$fits$fluA_only_full,
  "Influenza B" = bundle$fits$fluB_only_full
)

extract_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, variance = as.matrix(vcov(fit)$cond))
}
estimate_basis <- function(fit, basis, sign = 1) {
  fixed <- extract_fixed(fit)
  design <- matrix(
    0, nrow = nrow(basis), ncol = length(fixed$beta),
    dimnames = list(NULL, names(fixed$beta))
  )
  design[, cb_names] <- sign * basis
  design[, z_int_names] <-
    sweep(sign * basis, 2, interaction_scale, "/") * altitude_1000_c
  log_estimate <- as.numeric(design %*% fixed$beta)
  standard_error <- sqrt(pmax(
    rowSums((design %*% fixed$variance) * design), 0
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error,
    p_value = 2 * pnorm(-abs(log_estimate / standard_error))
  )
}

basis_matrix_cumulative <- function(temperatures, centre) {
  original_names <- colnames(cb_all)
  output <- matrix(
    0, length(temperatures), length(original_names),
    dimnames = list(as.character(temperatures), cb_names)
  )
  zero_vcov <- matrix(
    0, length(original_names), length(original_names),
    dimnames = list(original_names, original_names)
  )
  for (j in seq_along(original_names)) {
    coefficient <- rep(0, length(original_names))
    names(coefficient) <- original_names
    coefficient[[j]] <- 1
    prediction <- crosspred(
      cb_all, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = temperatures, cen = centre,
      cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

# Lag-specific estimates were already independently reconstructed and QA'd for
# the main Result 3 figure. Retain only the two prespecified contrasts.
lag_surfaces <- read.csv(
  surface_path, stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) |>
  filter(outcome %in% c(
    "Influenza A-only positivity", "Influenza B-only positivity"
  )) |>
  mutate(
    outcome = recode(
      outcome,
      `Influenza A-only positivity` = "Influenza A",
      `Influenza B-only positivity` = "Influenza B"
    )
  )

lag_minus5 <- lag_surfaces |>
  filter(abs(temperature + 5) < 1e-10) |>
  transmute(
    outcome,
    temperature_comparison = "-5 C vs 0 C",
    lag_window = paste0("Lag ", lag_week),
    lag_order = lag_week,
    altitude = altitude,
    effect_measure = "OR",
    estimate, low, high,
    standard_error,
    p_value = 2 * pnorm(-abs(log_estimate / standard_error))
  )
lag_0_vs5 <- lag_surfaces |>
  filter(abs(temperature - 5) < 1e-10) |>
  mutate(
    original_estimate = estimate,
    original_low = low,
    original_high = high
  ) |>
  transmute(
    outcome,
    temperature_comparison = "0 C vs 5 C",
    lag_window = paste0("Lag ", lag_week),
    lag_order = lag_week,
    altitude = altitude,
    effect_measure = "OR",
    estimate = 1 / original_estimate,
    low = 1 / original_high,
    high = 1 / original_low,
    standard_error,
    p_value = 2 * pnorm(-abs(log_estimate / standard_error))
  )

cumulative_basis <- basis_matrix_cumulative(
  c(-5, 0, 5), reference_temperature
)
cumulative_rows <- bind_rows(lapply(names(fits), function(outcome_name) {
  minus5_basis <- matrix(
    cumulative_basis["-5", ], nrow = 1,
    dimnames = list(NULL, cb_names)
  )
  plus5_basis <- matrix(
    cumulative_basis["5", ], nrow = 1,
    dimnames = list(NULL, cb_names)
  )
  bind_rows(
    bind_cols(
      tibble(
        outcome = outcome_name,
        temperature_comparison = "-5 C vs 0 C",
        lag_window = "Cumulative lag 0-2",
        lag_order = 3L,
        altitude = altitude,
        effect_measure = "OR"
      ),
      estimate_basis(fits[[outcome_name]], minus5_basis, sign = 1)
    ),
    bind_cols(
      tibble(
        outcome = outcome_name,
        temperature_comparison = "0 C vs 5 C",
        lag_window = "Cumulative lag 0-2",
        lag_order = 3L,
        altitude = altitude,
        effect_measure = "OR"
      ),
      estimate_basis(fits[[outcome_name]], plus5_basis, sign = -1)
    )
  )
})) |>
  select(
    outcome, temperature_comparison, lag_window, lag_order, altitude,
    effect_measure, estimate, low, high, standard_error, p_value
  )

lag_effects <- bind_rows(lag_minus5, lag_0_vs5, cumulative_rows) |>
  mutate(
    outcome = factor(outcome, levels = c("Influenza A", "Influenza B")),
    temperature_comparison = factor(
      temperature_comparison,
      levels = c("-5 C vs 0 C", "0 C vs 5 C")
    )
  ) |>
  arrange(outcome, temperature_comparison, lag_order) |>
  mutate(
    outcome = as.character(outcome),
    temperature_comparison = as.character(temperature_comparison)
  )

curves <- read.csv(
  file.path(phenology_dir, "standardized_subtype_seasonal_curves.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
peaks <- read.csv(
  file.path(phenology_dir, "subtype_peak_week_sensitivity.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
diagnostics <- read.csv(
  file.path(phenology_dir, "model_diagnostics.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

peak_estimates <- peaks |>
  left_join(
    curves |>
      select(
        outcome, harmonics, annual_week,
        peak_positivity = estimate,
        peak_positivity_low = low,
        peak_positivity_high = high
      ),
    by = c(
      "outcome", "harmonics", "peak_week" = "annual_week"
    )
  ) |>
  mutate(
    across(
      c(peak_positivity, peak_positivity_low, peak_positivity_high),
      ~ .x * 100
    )
  ) |>
  left_join(
    diagnostics |>
      select(outcome, harmonics, observations, convergence, pdHess, AIC),
    by = c("outcome", "harmonics")
  ) |>
  transmute(
    subtype = outcome_label,
    fourier_pairs = harmonics,
    specification = if_else(primary_specification, "Primary", "Sensitivity"),
    peak_week,
    peak_week_low,
    peak_week_high,
    peak_positivity,
    peak_positivity_low,
    peak_positivity_high,
    AIC,
    observations,
    convergence,
    pdHess
  ) |>
  arrange(subtype, fourier_pairs)

qa <- tibble(
  check = c(
    "Lag table contains two subtypes",
    "Lag table contains two temperature contrasts",
    "Each subtype-contrast has three lags and one cumulative row",
    "All effect estimates and intervals are finite",
    "Peak sensitivity covers 2-4 Fourier pairs",
    "All phenology models converged with positive-definite Hessians",
    "Primary peak values reproduce the main Result 3 source data"
  ),
  pass = c(
    n_distinct(lag_effects$outcome) == 2,
    n_distinct(lag_effects$temperature_comparison) == 2,
    all(lag_effects |>
          count(outcome, temperature_comparison) |>
          pull(n) == 4),
    all(is.finite(unlist(lag_effects[c("estimate", "low", "high", "p_value")]))),
    identical(sort(unique(peak_estimates$fourier_pairs)), 2:4),
    all(peak_estimates$convergence == 0 & peak_estimates$pdHess),
    all(abs(
      peak_estimates$peak_week[peak_estimates$specification == "Primary"] -
        c(51, 2.5)
    ) < 1e-10)
  )
)
if (any(!qa$pass)) stop("Supplementary Table S2 source data failed QA.")

write.csv(
  lag_effects,
  file.path(output_dir, "supplementary_table_s2_lag_effects.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  peak_estimates,
  file.path(output_dir, "supplementary_table_s2_peak_sensitivity.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  qa,
  file.path(output_dir, "supplementary_table_s2_source_QA.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

message("Supplementary Table S2 source data prepared in: ", output_dir)
