suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(glmmTMB)
})

# Confirmatory analysis contract
# 1) Test nonlinearity without selecting a temperature contrast from the result.
# 2) Test global altitude modification with an unpenalized low-df DLNM.
# 3) Estimate cold-side turning points only inside common P10-P90 support.
# 4) Treat a turning point as supported only when it is internal and its
#    uncertainty does not collapse onto a support boundary.

source(file.path("02_primary_models", "81_climate_context_models.R"))

analysis_id_confirm <- "confirmatory_anomaly_shape_20260831"
checkpoint_confirm <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", analysis_id_confirm
)
source_confirm <- file.path(
  "result0715", "source_data", paste0("v4_", analysis_id_confirm)
)
dir.create(checkpoint_confirm, recursive = TRUE, showWarnings = FALSE)
dir.create(source_confirm, recursive = TRUE, showWarnings = FALSE)

write_csv_utf8_confirm <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

control_confirm <- glmmTMBControl(
  optCtrl = list(iter.max = 6000, eval.max = 6000)
)

fit_or_load_confirm <- function(path, formula, data) {
  if (file.exists(path)) return(readRDS(path))
  fit <- glmmTMB(
    formula, family = betabinomial(link = "logit"),
    data = data, control = control_confirm
  )
  saveRDS(fit, path)
  fit
}

# Low-df lag structure: a linear polynomial across lags 0-2 (two lag basis
# columns). Exposure-response alternatives are linear (one exposure basis) and
# natural spline df=3. The spline interaction therefore uses six, rather than
# twelve, altitude-interaction coefficients and no priors.
cb_linear_grid <- crossbasis(
  weekly_grid$anomaly_for_basis,
  lag = c(0, 2),
  argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = 1),
  group = weekly_grid$county
)
cb_spline_grid <- crossbasis(
  weekly_grid$anomaly_for_basis,
  lag = c(0, 2),
  argvar = list(fun = "ns", df = 3),
  arglag = list(fun = "poly", degree = 1),
  group = weekly_grid$county
)

linear_names <- paste0("cf_lin_", seq_len(ncol(cb_linear_grid)))
spline_names <- paste0("cf_ns_", seq_len(ncol(cb_spline_grid)))
linear_int_names <- paste0("cf_lin_alt_", seq_len(ncol(cb_linear_grid)))
spline_int_names <- paste0("cf_ns_alt_", seq_len(ncol(cb_spline_grid)))

linear_df <- as.data.frame(cb_linear_grid[anomaly_model_dat$exposure_row_id, , drop = FALSE])
spline_df <- as.data.frame(cb_spline_grid[anomaly_model_dat$exposure_row_id, , drop = FALSE])
names(linear_df) <- linear_names
names(spline_df) <- spline_names

confirm_dat <- bind_cols(
  anomaly_model_dat,
  linear_df,
  spline_df
) |>
  filter(
    complete.cases(across(all_of(c(linear_names, spline_names))))
  )

for (i in seq_along(linear_names)) {
  confirm_dat[[linear_int_names[[i]]]] <-
    confirm_dat[[linear_names[[i]]]] * confirm_dat$altitude_1000_c
}
for (i in seq_along(spline_names)) {
  confirm_dat[[spline_int_names[[i]]]] <-
    confirm_dat[[spline_names[[i]]]] * confirm_dat$altitude_1000_c
}

random_terms <- c("(1 | city)", "(1 | county_nested)")
common_rhs <- c("altitude_1000_c", portable_names, random_terms)

formula_from_terms <- function(terms) {
  as.formula(paste(
    "cbind(positive, tests - positive) ~",
    paste(terms, collapse = " + ")
  ))
}

f_linear_no_int <- formula_from_terms(c(linear_names, common_rhs))
f_linear_alt <- formula_from_terms(c(
  linear_names, linear_int_names, common_rhs
))
f_spline_no_int <- formula_from_terms(c(spline_names, common_rhs))
f_spline_alt <- formula_from_terms(c(
  spline_names, spline_int_names, common_rhs
))

fit_linear_no_int <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "linear_no_altitude_interaction.rds"),
  f_linear_no_int, confirm_dat
)
fit_linear_alt <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "linear_altitude_interaction.rds"),
  f_linear_alt, confirm_dat
)
fit_spline_no_int <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "spline_no_altitude_interaction.rds"),
  f_spline_no_int, confirm_dat
)
fit_spline_alt <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "spline_altitude_interaction.rds"),
  f_spline_alt, confirm_dat
)

safe_lrt <- function(reduced, full, label) {
  tab <- anova(reduced, full)
  tibble(
    test = label,
    reduced_logLik = as.numeric(logLik(reduced)),
    full_logLik = as.numeric(logLik(full)),
    df_difference = tab$Df[[2]] - tab$Df[[1]],
    likelihood_ratio = tab$Chisq[[2]],
    p_value = tab$`Pr(>Chisq)`[[2]],
    reduced_AIC = AIC(reduced), full_AIC = AIC(full),
    delta_AIC_full_minus_reduced = AIC(full) - AIC(reduced),
    reduced_BIC = BIC(reduced), full_BIC = BIC(full),
    delta_BIC_full_minus_reduced = BIC(full) - BIC(reduced)
  )
}

confirm_tests <- bind_rows(
  safe_lrt(
    fit_linear_no_int, fit_spline_no_int,
    "Nonlinearity without altitude interaction"
  ),
  safe_lrt(
    fit_linear_alt, fit_spline_alt,
    "Nonlinearity allowing altitude interaction"
  ),
  safe_lrt(
    fit_spline_no_int, fit_spline_alt,
    "Global spline-by-altitude interaction"
  ),
  safe_lrt(
    fit_linear_no_int, fit_linear_alt,
    "Linear anomaly-by-altitude interaction"
  )
)

model_diagnostics_confirm <- bind_rows(
  tibble(model = "linear_no_int", fit = list(fit_linear_no_int)),
  tibble(model = "linear_alt", fit = list(fit_linear_alt)),
  tibble(model = "spline_no_int", fit = list(fit_spline_no_int)),
  tibble(model = "spline_alt", fit = list(fit_spline_alt))
) |>
  rowwise() |>
  transmute(
    model,
    n = nobs(fit),
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit), BIC = BIC(fit),
    convergence = fit$fit$convergence,
    pdHess = fit$sdr$pdHess
  ) |>
  ungroup()

extract_fixed_confirm <- function(fit) {
  list(beta = fixef(fit)$cond, vcov = as.matrix(vcov(fit)$cond))
}

cumulative_basis_confirm <- function(values, centre = 0) {
  values <- sort(unique(c(values, centre)))
  original_names <- colnames(cb_spline_grid)
  output <- matrix(
    NA_real_, length(values), length(original_names),
    dimnames = list(as.character(values), spline_names)
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
      cb_spline_grid, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = values, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

curve_design_confirm <- function(values, altitude_1000_c, fit = fit_spline_alt) {
  basis <- cumulative_basis_confirm(values, 0)
  fixed <- extract_fixed_confirm(fit)
  design <- matrix(
    0, nrow(basis), length(fixed$beta),
    dimnames = list(rownames(basis), names(fixed$beta))
  )
  design[, spline_names] <- basis
  design[, spline_int_names] <- basis * altitude_1000_c
  design
}

estimate_curve_confirm <- function(values, altitude_1000_c, fit = fit_spline_alt) {
  design <- curve_design_confirm(values, altitude_1000_c, fit)
  fixed <- extract_fixed_confirm(fit)
  log_est <- as.numeric(design %*% fixed$beta)
  se <- sqrt(pmax(rowSums((design %*% fixed$vcov) * design), 0))
  tibble(
    anomaly = as.numeric(rownames(design)),
    estimate = exp(log_est),
    low = exp(log_est - 1.96 * se),
    high = exp(log_est + 1.96 * se),
    log_estimate = log_est,
    standard_error = se
  )
}

common_p10_low_confirm <- max(anomaly_support$p10)
common_p90_high_confirm <- min(anomaly_support$p90)
curve_grid_confirm <- seq(
  ceiling(common_p10_low_confirm * 50) / 50,
  floor(common_p90_high_confirm * 50) / 50,
  by = 0.02
)
cold_grid_confirm <- curve_grid_confirm[curve_grid_confirm <= 0]

anchor_specs_confirm <- anchors |>
  select(altitude_anchor, altitude, altitude_1000_c)

curves_confirm <- bind_rows(lapply(seq_len(nrow(anchor_specs_confirm)), function(i) {
  spec <- anchor_specs_confirm[i, ]
  estimate_curve_confirm(curve_grid_confirm, spec$altitude_1000_c) |>
    mutate(
      altitude_anchor = spec$altitude_anchor,
      altitude = spec$altitude,
      altitude_1000_c = spec$altitude_1000_c
    )
}))

# Coefficient-draw bootstrap for turning-point uncertainty. This avoids treating
# a visual maximum from the fitted line as known. It is conditional on the
# fitted model; leave-one-year-out refits below assess temporal sensitivity.
set.seed(20260831)
n_draws <- 4000L
fixed_full <- extract_fixed_confirm(fit_spline_alt)
beta_draws <- MASS::mvrnorm(
  n = n_draws, mu = fixed_full$beta, Sigma = fixed_full$vcov
)
colnames(beta_draws) <- names(fixed_full$beta)

peak_by_anchor <- vector("list", nrow(anchor_specs_confirm))
peak_matrix <- matrix(
  NA_real_, nrow(anchor_specs_confirm), n_draws,
  dimnames = list(anchor_specs_confirm$altitude_anchor, NULL)
)

for (i in seq_len(nrow(anchor_specs_confirm))) {
  spec <- anchor_specs_confirm[i, ]
  design <- curve_design_confirm(cold_grid_confirm, spec$altitude_1000_c)
  eta <- design %*% t(beta_draws)
  peak_index <- max.col(t(eta), ties.method = "first")
  peak_value <- cold_grid_confirm[peak_index]
  peak_eta <- eta[cbind(peak_index, seq_len(n_draws))]
  peak_matrix[i, ] <- peak_value
  internal <- peak_index > 1L & peak_index < length(cold_grid_confirm)
  peak_by_anchor[[i]] <- tibble(
    altitude_anchor = spec$altitude_anchor,
    altitude = spec$altitude,
    median_peak_anomaly = median(peak_value),
    peak_anomaly_low = quantile(peak_value, 0.025),
    peak_anomaly_high = quantile(peak_value, 0.975),
    probability_internal_peak = mean(internal),
    median_peak_or = median(exp(peak_eta)),
    peak_or_low = quantile(exp(peak_eta), 0.025),
    peak_or_high = quantile(exp(peak_eta), 0.975),
    probability_peak_above_reference = mean(peak_eta > 0)
  )
}
turning_points_confirm <- bind_rows(peak_by_anchor)

# Full-common-support maxima. These complement the prespecified cold-side
# analysis above and prevent a positive-anomaly maximum at high altitude from
# being mislabeled as an absence of a response maximum.
full_peak_by_anchor <- vector("list", nrow(anchor_specs_confirm))
full_peak_matrix <- matrix(
  NA_real_, nrow(anchor_specs_confirm), n_draws,
  dimnames = list(anchor_specs_confirm$altitude_anchor, NULL)
)

for (i in seq_len(nrow(anchor_specs_confirm))) {
  spec <- anchor_specs_confirm[i, ]
  design <- curve_design_confirm(curve_grid_confirm, spec$altitude_1000_c)
  eta <- design %*% t(beta_draws)
  peak_index <- max.col(t(eta), ties.method = "first")
  peak_value <- curve_grid_confirm[peak_index]
  peak_eta <- eta[cbind(peak_index, seq_len(n_draws))]
  full_peak_matrix[i, ] <- peak_value
  internal <- peak_index > 1L & peak_index < length(curve_grid_confirm)
  full_peak_by_anchor[[i]] <- tibble(
    altitude_anchor = spec$altitude_anchor,
    altitude = spec$altitude,
    median_peak_anomaly = median(peak_value),
    peak_anomaly_low = quantile(peak_value, 0.025),
    peak_anomaly_high = quantile(peak_value, 0.975),
    probability_internal_peak = mean(internal),
    probability_peak_below_zero = mean(peak_value < 0),
    probability_peak_above_zero = mean(peak_value > 0),
    median_peak_or = median(exp(peak_eta)),
    peak_or_low = quantile(exp(peak_eta), 0.025),
    peak_or_high = quantile(exp(peak_eta), 0.975),
    probability_peak_above_reference = mean(peak_eta > 0)
  )
}
turning_points_full_confirm <- bind_rows(full_peak_by_anchor)

alt_km <- anchor_specs_confirm$altitude / 1000
alt_centered <- alt_km - mean(alt_km)
peak_shift_draws <- colSums(peak_matrix * alt_centered) /
  sum(alt_centered^2)
peak_shift_confirm <- tibble(
  estimand = "Change in cold-side peak anomaly per 1,000 m",
  median = median(peak_shift_draws),
  low = quantile(peak_shift_draws, 0.025),
  high = quantile(peak_shift_draws, 0.975),
  probability_shift_toward_colder_anomaly = mean(peak_shift_draws < 0)
)

full_peak_shift_draws <- colSums(full_peak_matrix * alt_centered) /
  sum(alt_centered^2)
peak_shift_full_confirm <- tibble(
  estimand = "Change in full-support peak anomaly per 1,000 m",
  median = median(full_peak_shift_draws),
  low = quantile(full_peak_shift_draws, 0.025),
  high = quantile(full_peak_shift_draws, 0.975),
  probability_shift_toward_colder_anomaly = mean(full_peak_shift_draws < 0),
  probability_shift_toward_warmer_anomaly = mean(full_peak_shift_draws > 0)
)

# Cold-warm asymmetry at +/-2 C under the low-df model.
contrast_design_confirm <- function(value, reference, altitude_1000_c) {
  basis <- cumulative_basis_confirm(c(value, reference), 0)
  diff <- basis[as.character(value), ] - basis[as.character(reference), ]
  fixed <- extract_fixed_confirm(fit_spline_alt)
  design <- rep(0, length(fixed$beta))
  names(design) <- names(fixed$beta)
  design[spline_names] <- diff
  design[spline_int_names] <- diff * altitude_1000_c
  design
}

estimate_design_confirm <- function(design) {
  fixed <- extract_fixed_confirm(fit_spline_alt)
  log_est <- sum(design * fixed$beta)
  se <- sqrt(as.numeric(t(design) %*% fixed$vcov %*% design))
  tibble(
    estimate = exp(log_est), low = exp(log_est - 1.96 * se),
    high = exp(log_est + 1.96 * se), log_estimate = log_est,
    standard_error = se
  )
}

signed_contrasts_confirm <- bind_rows(lapply(
  seq_len(nrow(anchor_specs_confirm)), function(i) {
    spec <- anchor_specs_confirm[i, ]
    bind_rows(
      estimate_design_confirm(
        contrast_design_confirm(-2, 0, spec$altitude_1000_c)
      ) |> mutate(contrast = "-2 C versus 0 C"),
      estimate_design_confirm(
        contrast_design_confirm(2, 0, spec$altitude_1000_c)
      ) |> mutate(contrast = "+2 C versus 0 C"),
      estimate_design_confirm(
        contrast_design_confirm(-2, 2, spec$altitude_1000_c)
      ) |> mutate(contrast = "-2 C versus +2 C")
    ) |>
      mutate(
        altitude_anchor = spec$altitude_anchor,
        altitude = spec$altitude
      )
  }
))

# Leave-one-year-out refits of the same confirmatory spline model. Each fold
# reports the likelihood-ratio test and the cold-side peak estimates.
loyo_results <- vector("list", 8L)
loyo_peaks <- vector("list", 8L)
years_confirm <- 2017:2024
for (j in seq_along(years_confirm)) {
  omitted <- years_confirm[[j]]
  fold_dat <- confirm_dat |> filter(iso_year != omitted)
  reduced_path <- file.path(
    checkpoint_confirm, sprintf("loyo_%d_spline_no_int.rds", omitted)
  )
  full_path <- file.path(
    checkpoint_confirm, sprintf("loyo_%d_spline_alt.rds", omitted)
  )
  reduced_fit <- fit_or_load_confirm(reduced_path, f_spline_no_int, fold_dat)
  full_fit <- fit_or_load_confirm(full_path, f_spline_alt, fold_dat)
  loyo_results[[j]] <- safe_lrt(
    reduced_fit, full_fit,
    sprintf("Global interaction; omitted %d", omitted)
  ) |>
    mutate(
      omitted_year = omitted,
      convergence_reduced = reduced_fit$fit$convergence,
      convergence_full = full_fit$fit$convergence,
      pdHess_reduced = reduced_fit$sdr$pdHess,
      pdHess_full = full_fit$sdr$pdHess
    )
  loyo_peaks[[j]] <- bind_rows(lapply(
    seq_len(nrow(anchor_specs_confirm)), function(i) {
      spec <- anchor_specs_confirm[i, ]
      curve <- estimate_curve_confirm(
        cold_grid_confirm, spec$altitude_1000_c, full_fit
      )
      top <- curve[which.max(curve$estimate), ]
      tibble(
        omitted_year = omitted,
        altitude_anchor = spec$altitude_anchor,
        altitude = spec$altitude,
        cold_side_peak_anomaly = top$anomaly,
        peak_or = top$estimate
      )
    }
  ))
}
loyo_interaction_confirm <- bind_rows(loyo_results)
loyo_peaks_confirm <- bind_rows(loyo_peaks)

# Province-week and county fixed-effects sensitivity. Local climatology remains
# in the model because it varies by county-week and is not absorbed by additive
# county and province-week effects.
# The altitude main effect is time-invariant within county and is therefore
# exactly absorbed by county fixed effects. Retain only anomaly-by-altitude
# interaction terms in the full fixed-effects model.
fixed_common <- c(climatology_names, "county", "week_id")
f_fe_no_int <- formula_from_terms(c(spline_names, fixed_common))
f_fe_alt <- formula_from_terms(c(
  spline_names, spline_int_names, fixed_common
))
fit_fe_no_int <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "twoway_fixed_effects_no_int_v2.rds"),
  f_fe_no_int, confirm_dat
)
fit_fe_alt <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "twoway_fixed_effects_alt_v2.rds"),
  f_fe_alt, confirm_dat
)
fixed_effects_test_confirm <- safe_lrt(
  fit_fe_no_int, fit_fe_alt,
  "Global interaction with county and province-week fixed effects"
) |>
  mutate(
    convergence_reduced = fit_fe_no_int$fit$convergence,
    convergence_full = fit_fe_alt$fit$convergence,
    pdHess_reduced = fit_fe_no_int$sdr$pdHess,
    pdHess_full = fit_fe_alt$sdr$pdHess
  )

# Absolute-humidity anomaly sensitivity on its complete-data subset.
weekly_grid_ah <- weekly_grid |>
  group_by(county, iso_week) |>
  mutate(
    ah_valid_here = is.finite(absolute_humidity_g_m3) & ah_valid_days >= 5,
    ah_total = sum(if_else(ah_valid_here, absolute_humidity_g_m3, 0), na.rm = TRUE),
    ah_n = sum(ah_valid_here, na.rm = TRUE),
    ah_other_n = ah_n - as.integer(ah_valid_here),
    loyo_ah_climatology = if_else(
      ah_other_n >= 1,
      (ah_total - if_else(ah_valid_here, absolute_humidity_g_m3, 0)) / ah_other_n,
      NA_real_
    ),
    ah_anomaly_for_basis = if_else(
      ah_valid_here & is.finite(loyo_ah_climatology),
      absolute_humidity_g_m3 - loyo_ah_climatology,
      NA_real_
    )
  ) |>
  ungroup()

cb_ah_grid <- crossbasis(
  weekly_grid_ah$ah_anomaly_for_basis,
  lag = c(0, 2),
  argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = 1),
  group = weekly_grid_ah$county
)
ah_names <- paste0("cf_ah_", seq_len(ncol(cb_ah_grid)))
ah_df <- as.data.frame(cb_ah_grid[confirm_dat$exposure_row_id, , drop = FALSE])
names(ah_df) <- ah_names
ah_clim <- weekly_grid_ah$loyo_ah_climatology[confirm_dat$exposure_row_id]
ah_clim_basis <- splines::ns(
  ah_clim, df = 3,
  Boundary.knots = range(ah_clim, finite = TRUE)
)
ah_clim_names <- paste0("cf_ah_clim_", seq_len(ncol(ah_clim_basis)))
colnames(ah_clim_basis) <- ah_clim_names

humidity_dat <- bind_cols(
  confirm_dat,
  ah_df,
  as.data.frame(ah_clim_basis)
) |>
  filter(complete.cases(across(all_of(c(ah_names, ah_clim_names)))))

humidity_common <- c(common_rhs, ah_names, ah_clim_names)
f_humidity_no_int <- formula_from_terms(c(spline_names, humidity_common))
f_humidity_alt <- formula_from_terms(c(
  spline_names, spline_int_names, humidity_common
))
fit_humidity_no_int <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "absolute_humidity_no_int.rds"),
  f_humidity_no_int, humidity_dat
)
fit_humidity_alt <- fit_or_load_confirm(
  file.path(checkpoint_confirm, "absolute_humidity_alt.rds"),
  f_humidity_alt, humidity_dat
)
humidity_test_confirm <- safe_lrt(
  fit_humidity_no_int, fit_humidity_alt,
  "Global interaction adjusted for absolute-humidity anomaly"
) |>
  mutate(
    n = nobs(fit_humidity_alt),
    convergence_reduced = fit_humidity_no_int$fit$convergence,
    convergence_full = fit_humidity_alt$fit$convergence,
    pdHess_reduced = fit_humidity_no_int$sdr$pdHess,
    pdHess_full = fit_humidity_alt$sdr$pdHess
  )

write_csv_utf8_confirm(
  confirm_tests,
  file.path(source_confirm, "confirmatory_model_tests.csv")
)
write_csv_utf8_confirm(
  model_diagnostics_confirm,
  file.path(source_confirm, "confirmatory_model_diagnostics.csv")
)
write_csv_utf8_confirm(
  curves_confirm,
  file.path(source_confirm, "confirmatory_anomaly_curves.csv")
)
write_csv_utf8_confirm(
  turning_points_confirm,
  file.path(source_confirm, "confirmatory_cold_side_turning_points.csv")
)
write_csv_utf8_confirm(
  turning_points_full_confirm,
  file.path(source_confirm, "confirmatory_full_support_turning_points.csv")
)
write_csv_utf8_confirm(
  peak_shift_confirm,
  file.path(source_confirm, "confirmatory_peak_shift.csv")
)
write_csv_utf8_confirm(
  peak_shift_full_confirm,
  file.path(source_confirm, "confirmatory_full_support_peak_shift.csv")
)
write_csv_utf8_confirm(
  signed_contrasts_confirm,
  file.path(source_confirm, "confirmatory_signed_contrasts.csv")
)
write_csv_utf8_confirm(
  loyo_interaction_confirm,
  file.path(source_confirm, "confirmatory_leave_one_year_out_interaction.csv")
)
write_csv_utf8_confirm(
  loyo_peaks_confirm,
  file.path(source_confirm, "confirmatory_leave_one_year_out_peaks.csv")
)
write_csv_utf8_confirm(
  fixed_effects_test_confirm,
  file.path(source_confirm, "confirmatory_twoway_fixed_effects.csv")
)
write_csv_utf8_confirm(
  humidity_test_confirm,
  file.path(source_confirm, "confirmatory_absolute_humidity_adjustment.csv")
)

saveRDS(
  list(
    confirm_dat = confirm_dat,
    cb_linear_grid = cb_linear_grid,
    cb_spline_grid = cb_spline_grid,
    linear_names = linear_names,
    spline_names = spline_names,
    linear_int_names = linear_int_names,
    spline_int_names = spline_int_names,
    fit_linear_no_int = fit_linear_no_int,
    fit_linear_alt = fit_linear_alt,
    fit_spline_no_int = fit_spline_no_int,
    fit_spline_alt = fit_spline_alt,
    confirm_tests = confirm_tests,
    curves_confirm = curves_confirm,
    turning_points_confirm = turning_points_confirm,
    turning_points_full_confirm = turning_points_full_confirm,
    peak_shift_confirm = peak_shift_confirm,
    peak_shift_full_confirm = peak_shift_full_confirm,
    signed_contrasts_confirm = signed_contrasts_confirm,
    loyo_interaction_confirm = loyo_interaction_confirm,
    loyo_peaks_confirm = loyo_peaks_confirm,
    fixed_effects_test_confirm = fixed_effects_test_confirm,
    humidity_test_confirm = humidity_test_confirm
  ),
  file.path(checkpoint_confirm, "confirmatory_anomaly_bundle.rds")
)

qa_confirm <- tibble(
  check = c(
    "All four primary models converged",
    "All four primary Hessians are positive definite",
    "Confirmatory model is unpenalized",
    "Turning points restricted to common P10-P90 support",
    "Full-support maxima restricted to common P10-P90 support",
    "Full-support peak-side probabilities are valid",
    "All leave-one-year-out models converged",
    "Two-way fixed-effects models converged",
    "Humidity-adjusted models converged"
  ),
  pass = c(
    all(model_diagnostics_confirm$convergence == 0),
    all(model_diagnostics_confirm$pdHess),
    is.null(fit_spline_alt$modelInfo$priors) ||
      nrow(fit_spline_alt$modelInfo$priors) == 0,
    all(turning_points_confirm$peak_anomaly_low >= common_p10_low_confirm) &&
      all(turning_points_confirm$peak_anomaly_high <= 0),
    all(turning_points_full_confirm$peak_anomaly_low >=
          min(curve_grid_confirm)) &&
      all(turning_points_full_confirm$peak_anomaly_high <=
            max(curve_grid_confirm)),
    all(turning_points_full_confirm$probability_peak_below_zero >= 0 &
          turning_points_full_confirm$probability_peak_below_zero <= 1),
    all(loyo_interaction_confirm$convergence_reduced == 0) &&
      all(loyo_interaction_confirm$convergence_full == 0),
    fixed_effects_test_confirm$convergence_reduced == 0 &&
      fixed_effects_test_confirm$convergence_full == 0,
    humidity_test_confirm$convergence_reduced == 0 &&
      humidity_test_confirm$convergence_full == 0
  )
)
write_csv_utf8_confirm(
  qa_confirm,
  file.path(source_confirm, "confirmatory_analysis_qa.csv")
)

message("Completed low-df unpenalized confirmatory anomaly-shape analysis.")
