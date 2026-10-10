suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(glmmTMB)
  library(geepack)
})

set.seed(20260901)

# Result 5 contract
# - Use only the updated 2-m air-temperature data.
# - Preserve the primary absolute-temperature estimand (-5 C versus 0 C).
# - Use temperature anomalies for local-context robustness analyses.
# - Keep sensitivity, falsification and transportability diagnostics distinct.

source(file.path("02_primary_models", "80_primary_models.R"))

result5_id <- "result5_robustness_20260901"
result5_checkpoint <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", result5_id
)
result5_source <- file.path(
  "result0715", "source_data", paste0("v4_", result5_id)
)
dir.create(result5_checkpoint, recursive = TRUE, showWarnings = FALSE)
dir.create(result5_source, recursive = TRUE, showWarnings = FALSE)

write_result5_csv <- function(x, filename) {
  write.csv(
    x, file.path(result5_source, filename), row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

primary_bundle <- readRDS(file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825", "v4_t2m_model_bundle.rds"
))

primary_dat <- primary_bundle$model_data
primary_cb <- primary_bundle$cb_all
primary_cb_names <- primary_bundle$cb_names
primary_int_names <- primary_bundle$z_int_names
primary_scale <- primary_bundle$interaction_scale
primary_anchors <- primary_bundle$altitude_anchors
primary_altitude_difference <-
  primary_anchors$altitude_1000_c[primary_anchors$altitude_anchor == "P90"] -
  primary_anchors$altitude_1000_c[primary_anchors$altitude_anchor == "P10"]

primary_basis_cumulative <- function(values, centre = 0) {
  values <- sort(unique(c(values, centre)))
  original_names <- colnames(primary_cb)
  output <- matrix(
    NA_real_, length(values), length(original_names),
    dimnames = list(as.character(values), primary_cb_names)
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
      primary_cb, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = values, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

extract_glmm_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, variance = as.matrix(vcov(fit)$cond))
}

extract_gee_fixed <- function(fit) {
  beta <- coef(fit)
  variance <- fit$geese$vbeta
  dimnames(variance) <- list(names(beta), names(beta))
  list(beta = beta, variance = variance)
}

estimate_primary_contrast_generic <- function(
    fixed, colder, warmer, altitude_1000_c,
    interaction_only = FALSE) {
  basis <- primary_basis_cumulative(c(colder, warmer), 0)
  difference <- basis[as.character(colder), ] - basis[as.character(warmer), ]
  design <- rep(0, length(fixed$beta))
  names(design) <- names(fixed$beta)
  if (!interaction_only) design[primary_cb_names] <- difference
  design[primary_int_names] <-
    altitude_1000_c * difference / primary_scale
  log_estimate <- sum(design * fixed$beta)
  standard_error <- sqrt(as.numeric(
    t(design) %*% fixed$variance %*% design
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

primary_absolute <- estimate_primary_contrast_generic(
  extract_glmm_fixed(primary_bundle$fits$positive_full),
  -5, 0, 1, TRUE
) |>
  mutate(
    analysis = "Primary full laboratory panel",
    outcome = "Laboratory positivity",
    estimand = "Ratio of ORs per 1,000 m for -5 C versus 0 C",
    observations = nobs(primary_bundle$fits$positive_full)
  )

excluded_absolute <- estimate_primary_contrast_generic(
  extract_glmm_fixed(primary_bundle$fits$positive_excluded),
  -5, 0, 1, TRUE
) |>
  mutate(
    analysis = "Spring-2023 weeks excluded",
    outcome = "Laboratory positivity",
    estimand = "Ratio of ORs per 1,000 m for -5 C versus 0 C",
    observations = nobs(primary_bundle$fits$positive_excluded)
  )

# Monitoring-geography sensitivity: retain the same province-wide P10 and P90
# altitude anchors, but restrict the laboratory outcome to the 18 counties with
# ILI surveillance. This isolates geography from a change in the estimand.
ili_counties <- primary_dat |>
  filter(ili_available %in% TRUE) |>
  distinct(county) |>
  pull(county) |>
  as.character()

restricted_dat <- primary_dat |>
  filter(as.character(county) %in% ili_counties) |>
  droplevels()

restricted_path <- file.path(
  result5_checkpoint, "laboratory_ili_footprint_hierarchical.rds"
)
if (file.exists(restricted_path)) {
  fit_restricted <- readRDS(restricted_path)
} else {
  fit_restricted <- glmmTMB(
    formula(primary_bundle$fits$positive_full),
    family = betabinomial(link = "logit"), data = restricted_dat,
    priors = primary_bundle$fits$positive_full$modelInfo$priors,
    control = glmmTMBControl(optCtrl = list(iter.max = 5000, eval.max = 5000))
  )
  saveRDS(fit_restricted, restricted_path)
}

restricted_absolute <- estimate_primary_contrast_generic(
  extract_glmm_fixed(fit_restricted),
  -5, 0, 1, TRUE
) |>
  mutate(
    analysis = "Laboratory panel restricted to the ILI footprint",
    outcome = "Laboratory positivity",
    estimand = "Ratio of ORs per 1,000 m for -5 C versus 0 C",
    observations = nobs(fit_restricted)
  )

# Alternative population-averaged estimator for the full laboratory panel.
gee_dat <- primary_dat |>
  arrange(county, week_index) |>
  mutate(positivity = positive / tests, county = factor(county))
gee_terms <- c(
  primary_cb_names, "altitude_1000_c", primary_int_names,
  grep("^season_(sin|cos)_", names(gee_dat), value = TRUE),
  grep("^portable_trend_", names(gee_dat), value = TRUE),
  "pandemic_period"
)
gee_formula <- as.formula(paste(
  "positivity ~", paste(gee_terms, collapse = " + ")
))
gee_path <- file.path(result5_checkpoint, "laboratory_full_gee_ar1.rds")
if (file.exists(gee_path)) {
  fit_lab_gee <- readRDS(gee_path)
} else {
  fit_lab_gee <- geeglm(
    gee_formula, family = binomial(link = "logit"), data = gee_dat,
    weights = tests, id = county, waves = week_index,
    corstr = "ar1", std.err = "san.se"
  )
  saveRDS(fit_lab_gee, gee_path)
}

gee_absolute <- estimate_primary_contrast_generic(
  extract_gee_fixed(fit_lab_gee),
  -5, 0, 1, TRUE
) |>
  mutate(
    analysis = "Population-averaged GEE-AR1",
    outcome = "Laboratory positivity",
    estimand = "Ratio of ORs per 1,000 m for -5 C versus 0 C",
    observations = nrow(gee_dat)
  )

absolute_temperature_results <- bind_rows(
  primary_absolute, restricted_absolute, excluded_absolute, gee_absolute
) |>
  select(
    analysis, outcome, estimand, estimate, low, high,
    log_estimate, standard_error, observations
  )

# Endpoint triangulation. The ILI models use the identical exposure basis,
# lag window and temporal adjustment but a different outcome and surveillance
# footprint. They are not treated as numerical validation of positivity.
ili_cb_df <- as.data.frame(cb_all)
names(ili_cb_df) <- cb_names
ili_dat <- bind_cols(
  weekly_grid |> mutate(exposure_row_id = row_number()), ili_cb_df
) |>
  filter(
    ili_available %in% TRUE, visits >= 1, ILI <= visits,
    complete.cases(across(all_of(cb_names))),
    !is.na(altitude_1000_c), !is.na(city)
  ) |>
  mutate(
    city = factor(city), county = factor(county),
    county_nested = interaction(city, county, drop = TRUE)
  ) |>
  arrange(county, week_index)

for (harmonic in 1:3) {
  ili_dat[[paste0("season_sin_", harmonic)]] <-
    sin(2 * pi * harmonic * ili_dat$week_index / annual_period_weeks)
  ili_dat[[paste0("season_cos_", harmonic)]] <-
    cos(2 * pi * harmonic * ili_dat$week_index / annual_period_weeks)
}
ili_trend_basis <- splines::ns(
  ili_dat$week_index, df = 4,
  Boundary.knots = range(primary_dat$week_index)
)
colnames(ili_trend_basis) <- trend_names
ili_dat <- bind_cols(ili_dat, as.data.frame(ili_trend_basis)) |>
  mutate(pandemic_period = as.integer(iso_year %in% 2020:2022))

for (i in seq_along(cb_names)) {
  ili_dat[[raw_int_names[[i]]]] <-
    ili_dat[[cb_names[[i]]]] * ili_dat$altitude_1000_c
  ili_dat[[z_int_names[[i]]]] <-
    (ili_dat[[raw_int_names[[i]]]] - interaction_center[[i]]) /
    interaction_scale[[i]]
}

ili_formula_hierarchical <- as.formula(paste(
  "cbind(ILI, visits - ILI) ~",
  paste(base_rhs, collapse = " + ")
))
ili_hier_path <- file.path(result5_checkpoint, "ili_hierarchical.rds")
if (file.exists(ili_hier_path)) {
  fit_ili_hier <- readRDS(ili_hier_path)
} else {
  fit_ili_hier <- glmmTMB(
    ili_formula_hierarchical,
    family = betabinomial(link = "logit"), data = ili_dat,
    priors = make_priors(), control = control
  )
  saveRDS(fit_ili_hier, ili_hier_path)
}

ili_median_altitude <- ili_dat |>
  distinct(county, altitude, altitude_1000_c) |>
  summarise(
    altitude = median(altitude),
    altitude_1000_c = (altitude - altitude_mean) / 1000
  )

matched_laboratory_contrast <- estimate_primary_contrast_generic(
  extract_glmm_fixed(fit_restricted), -5, 0,
  ili_median_altitude$altitude_1000_c, FALSE
) |>
  mutate(
    analysis = "Laboratory model in the matched ILI footprint",
    outcome = "Laboratory positivity",
    estimand = sprintf(
      "OR for -5 C versus 0 C at %.0f m", ili_median_altitude$altitude
    ),
    observations = nobs(fit_restricted)
  )

ili_hierarchical_contrast <- estimate_primary_contrast_generic(
  extract_glmm_fixed(fit_ili_hier), -5, 0,
  ili_median_altitude$altitude_1000_c, FALSE
) |>
  mutate(
    analysis = "Hierarchical beta-binomial model",
    outcome = "ILI consultation rate",
    estimand = sprintf(
      "OR for -5 C versus 0 C at %.0f m", ili_median_altitude$altitude
    ),
    observations = nobs(fit_ili_hier)
  )

ili_gee_dat <- ili_dat |>
  mutate(ili_rate = ILI / visits, county = factor(county))
ili_gee_formula <- as.formula(paste(
  "ili_rate ~",
  paste(c(
    cb_names, "altitude_1000_c", z_int_names,
    season_names, trend_names, "pandemic_period"
  ), collapse = " + ")
))
ili_gee_path <- file.path(result5_checkpoint, "ili_gee_ar1.rds")
if (file.exists(ili_gee_path)) {
  fit_ili_gee <- readRDS(ili_gee_path)
} else {
  fit_ili_gee <- geeglm(
    ili_gee_formula, family = binomial(link = "logit"),
    data = ili_gee_dat, weights = visits, id = county,
    waves = week_index, corstr = "ar1", std.err = "san.se"
  )
  saveRDS(fit_ili_gee, ili_gee_path)
}

ili_gee_contrast <- estimate_primary_contrast_generic(
  extract_gee_fixed(fit_ili_gee), -5, 0,
  ili_median_altitude$altitude_1000_c, FALSE
) |>
  mutate(
    analysis = "Population-averaged GEE-AR1",
    outcome = "ILI consultation rate",
    estimand = sprintf(
      "OR for -5 C versus 0 C at %.0f m", ili_median_altitude$altitude
    ),
    observations = nrow(ili_gee_dat)
  )

endpoint_results <- bind_rows(
  matched_laboratory_contrast,
  ili_hierarchical_contrast,
  ili_gee_contrast
) |>
  mutate(
    reportable = is.finite(estimate) & is.finite(low) & is.finite(high)
  ) |>
  select(
    analysis, outcome, estimand, estimate, low, high,
    log_estimate, standard_error, observations, reportable
  )

# Load the finalized confirmatory anomaly analysis, then add pressure,
# future-exposure and geographic-transportability checks.
source(file.path("02_primary_models", "82_confirmatory_anomaly_shape.R"))
source(file.path("02_primary_models", "83_confirmatory_twoway_gee.R"))

result5_checkpoint <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", result5_id
)
result5_source <- file.path(
  "result0715", "source_data", paste0("v4_", result5_id)
)

# Pressure sensitivity uses short-term pressure anomalies and the local pressure
# climatology, avoiding direct adjustment for raw pressure and altitude together.
weekly_grid_pressure <- weekly_grid |>
  group_by(county, iso_week) |>
  mutate(
    pressure_valid_here =
      is.finite(pressure_hpa) & pressure_valid_days >= 5,
    pressure_total = sum(
      if_else(pressure_valid_here, pressure_hpa, 0), na.rm = TRUE
    ),
    pressure_n = sum(pressure_valid_here, na.rm = TRUE),
    pressure_other_n = pressure_n - as.integer(pressure_valid_here),
    loyo_pressure_climatology = if_else(
      pressure_other_n >= 1,
      (pressure_total - if_else(pressure_valid_here, pressure_hpa, 0)) /
        pressure_other_n,
      NA_real_
    ),
    pressure_anomaly_for_basis = if_else(
      pressure_valid_here & is.finite(loyo_pressure_climatology),
      pressure_hpa - loyo_pressure_climatology,
      NA_real_
    )
  ) |>
  ungroup()

cb_pressure_grid <- crossbasis(
  weekly_grid_pressure$pressure_anomaly_for_basis,
  lag = c(0, 2), argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = 1),
  group = weekly_grid_pressure$county
)
pressure_names <- paste0("r5_pressure_", seq_len(ncol(cb_pressure_grid)))
pressure_df <- as.data.frame(
  cb_pressure_grid[confirm_dat$exposure_row_id, , drop = FALSE]
)
names(pressure_df) <- pressure_names
pressure_clim <- weekly_grid_pressure$loyo_pressure_climatology[
  confirm_dat$exposure_row_id
]
pressure_clim_basis <- splines::ns(
  pressure_clim, df = 3,
  Boundary.knots = range(pressure_clim, finite = TRUE)
)
pressure_clim_names <- paste0(
  "r5_pressure_clim_", seq_len(ncol(pressure_clim_basis))
)
colnames(pressure_clim_basis) <- pressure_clim_names
pressure_dat <- bind_cols(
  confirm_dat, pressure_df, as.data.frame(pressure_clim_basis)
) |>
  filter(complete.cases(across(all_of(c(
    pressure_names, pressure_clim_names
  )))))

pressure_common <- c(common_rhs, pressure_names, pressure_clim_names)
f_pressure_no_int <- formula_from_terms(c(spline_names, pressure_common))
f_pressure_alt <- formula_from_terms(c(
  spline_names, spline_int_names, pressure_common
))
fit_pressure_no_int <- fit_or_load_confirm(
  file.path(result5_checkpoint, "pressure_no_int.rds"),
  f_pressure_no_int, pressure_dat
)
fit_pressure_alt <- fit_or_load_confirm(
  file.path(result5_checkpoint, "pressure_alt.rds"),
  f_pressure_alt, pressure_dat
)
pressure_test <- safe_lrt(
  fit_pressure_no_int, fit_pressure_alt,
  "Global interaction adjusted for pressure anomaly"
) |>
  mutate(
    n = nobs(fit_pressure_alt),
    convergence_reduced = fit_pressure_no_int$fit$convergence,
    convergence_full = fit_pressure_alt$fit$convergence,
    pdHess_reduced = fit_pressure_no_int$sdr$pdHess,
    pdHess_full = fit_pressure_alt$sdr$pdHess
  )

# Negative-control exposure: lead the anomaly by three weeks, then apply lags
# 0-2. The resulting basis contains only future weeks t+1 through t+3.
weekly_grid_future <- weekly_grid |>
  group_by(county) |>
  arrange(week_start, .by_group = TRUE) |>
  mutate(future_anomaly_for_basis = lead(anomaly_for_basis, 3)) |>
  ungroup()

cb_future_grid <- crossbasis(
  weekly_grid_future$future_anomaly_for_basis,
  lag = c(0, 2), argvar = list(fun = "ns", df = 3),
  arglag = list(fun = "poly", degree = 1),
  group = weekly_grid_future$county
)
future_names <- paste0("r5_future_ns_", seq_len(ncol(cb_future_grid)))
future_int_names <- paste0("r5_future_alt_", seq_len(ncol(cb_future_grid)))
future_df <- as.data.frame(
  cb_future_grid[confirm_dat$exposure_row_id, , drop = FALSE]
)
names(future_df) <- future_names
future_dat <- bind_cols(confirm_dat, future_df) |>
  filter(complete.cases(across(all_of(future_names))))
for (i in seq_along(future_names)) {
  future_dat[[future_int_names[[i]]]] <-
    future_dat[[future_names[[i]]]] * future_dat$altitude_1000_c
}
f_future_no_int <- formula_from_terms(c(future_names, common_rhs))
f_future_alt <- formula_from_terms(c(
  future_names, future_int_names, common_rhs
))
fit_future_no_int <- fit_or_load_confirm(
  file.path(result5_checkpoint, "future_anomaly_no_int.rds"),
  f_future_no_int, future_dat
)
fit_future_alt <- fit_or_load_confirm(
  file.path(result5_checkpoint, "future_anomaly_alt.rds"),
  f_future_alt, future_dat
)
future_test <- safe_lrt(
  fit_future_no_int, fit_future_alt,
  "Future-anomaly spline-by-altitude interaction"
) |>
  mutate(
    n = nobs(fit_future_alt),
    convergence_reduced = fit_future_no_int$fit$convergence,
    convergence_full = fit_future_alt$fit$convergence,
    pdHess_reduced = fit_future_no_int$sdr$pdHess,
    pdHess_full = fit_future_alt$sdr$pdHess
  )

# Leave-one-prefecture-out predictive validation. Random effects are omitted
# when predicting held-out prefectures, so the diagnostic evaluates transport
# of fixed temperature-response structure rather than fitted local intercepts.
binomial_log_loss <- function(successes, totals, probability) {
  probability <- pmin(pmax(probability, 1e-8), 1 - 1e-8)
  -sum(
    successes * log(probability) +
      (totals - successes) * log(1 - probability)
  ) / sum(totals)
}

cities <- sort(unique(as.character(confirm_dat$city)))
lopo_rows <- vector("list", length(cities))
for (i in seq_along(cities)) {
  omitted_city <- cities[[i]]
  train <- confirm_dat |>
    filter(as.character(city) != omitted_city) |>
    droplevels()
  test <- confirm_dat |>
    filter(as.character(city) == omitted_city)
  safe_city <- sprintf("city_%02d", i)
  reduced_path <- file.path(
    result5_checkpoint, paste0("lopo_", safe_city, "_no_int.rds")
  )
  full_path <- file.path(
    result5_checkpoint, paste0("lopo_", safe_city, "_alt.rds")
  )
  reduced_fit <- fit_or_load_confirm(reduced_path, f_spline_no_int, train)
  full_fit <- fit_or_load_confirm(full_path, f_spline_alt, train)
  reduced_probability <- predict(
    reduced_fit, newdata = test, type = "response",
    re.form = NA, allow.new.levels = TRUE
  )
  full_probability <- predict(
    full_fit, newdata = test, type = "response",
    re.form = NA, allow.new.levels = TRUE
  )
  lopo_rows[[i]] <- tibble(
    omitted_prefecture = omitted_city,
    county_weeks = nrow(test), tests = sum(test$tests),
    reduced_log_loss = binomial_log_loss(
      test$positive, test$tests, reduced_probability
    ),
    interaction_log_loss = binomial_log_loss(
      test$positive, test$tests, full_probability
    ),
    interaction_minus_reduced = interaction_log_loss - reduced_log_loss,
    interaction_favoured = interaction_log_loss < reduced_log_loss,
    convergence_reduced = reduced_fit$fit$convergence,
    convergence_full = full_fit$fit$convergence,
    pdHess_reduced = reduced_fit$sdr$pdHess,
    pdHess_full = full_fit$sdr$pdHess
  )
}
lopo_results <- bind_rows(lopo_rows)

# Compact summaries for the main-text table and detailed supplementary sheets.
primary_global <- confirm_tests |>
  filter(test == "Global spline-by-altitude interaction")
humidity_global <- humidity_test_confirm
gee_fixed_global <- read.csv(file.path(
  "result0715", "source_data", "v4_confirmatory_anomaly_shape_20260831",
  "confirmatory_twoway_fixed_effects_gee.csv"
))
loyo_summary <- loyo_interaction_confirm |>
  summarise(
    folds = n(), p_lt_005 = sum(p_value < .05),
    aic_favours = sum(delta_AIC_full_minus_reduced < 0),
    bic_favours = sum(delta_BIC_full_minus_reduced < 0)
  )
lopo_summary <- lopo_results |>
  summarise(
    folds = n(),
    interaction_favours = sum(interaction_favoured),
    median_log_loss_difference = median(interaction_minus_reduced)
  )

phase_results <- read.csv(file.path(
  "result0715", "source_data", "v4_result4_climate_context_phase_20260826",
  "phase_specific_temperature_altitude_results.csv"
))
phase_tests <- read.csv(file.path(
  "result0715", "source_data", "v4_result4_climate_context_phase_20260826",
  "phase_heterogeneity_tests.csv"
))

main_table <- bind_rows(
  absolute_temperature_results |>
    transmute(
      domain = "Absolute-temperature robustness",
      analysis, estimand,
      result = sprintf("%.2f (%.2f-%.2f)", estimate, low, high),
      p_value = NA_real_, observations
    ),
  endpoint_results |>
    filter(reportable) |>
    transmute(
      domain = "Endpoint triangulation", analysis, estimand,
      result = sprintf("%.2f (%.2f-%.2f)", estimate, low, high),
      p_value = NA_real_, observations
    ),
  tibble(
    domain = "Local-anomaly model",
    analysis = c(
      "Primary hierarchical model",
      "Adjusted for absolute-humidity anomaly",
      "Adjusted for pressure anomaly",
      "County and province-week fixed-effects GEE"
    ),
    estimand = c(
      "Global spline-by-altitude interaction",
      "Global spline-by-altitude interaction",
      "Global spline-by-altitude interaction",
      "Joint Wald test of spline-by-altitude terms"
    ),
    result = c(
      sprintf(
        "Delta AIC %.1f; Delta BIC %.1f",
        primary_global$delta_AIC_full_minus_reduced,
        primary_global$delta_BIC_full_minus_reduced
      ),
      sprintf(
        "Delta AIC %.1f; Delta BIC %.1f",
        humidity_global$delta_AIC_full_minus_reduced,
        humidity_global$delta_BIC_full_minus_reduced
      ),
      sprintf(
        "Delta AIC %.1f; Delta BIC %.1f",
        pressure_test$delta_AIC_full_minus_reduced,
        pressure_test$delta_BIC_full_minus_reduced
      ),
      sprintf("Wald chi-square %.1f (%d df)",
              gee_fixed_global$statistic, gee_fixed_global$df)
    ),
    p_value = c(
      primary_global$p_value, humidity_global$p_value,
      pressure_test$p_value, gee_fixed_global$p_value
    ),
    observations = c(
      nobs(fit_spline_alt), humidity_global$n,
      pressure_test$n, gee_fixed_global$observations
    )
  ),
  tibble(
    domain = c(
      "Falsification", "Temporal repeatability",
      "Geographic transportability"
    ),
    analysis = c(
      "Future temperature anomaly (weeks t+1 to t+3)",
      "Leave-one-year-out refitting",
      "Leave-one-prefecture-out prediction"
    ),
    estimand = c(
      "Global spline-by-altitude interaction",
      "Interaction model favoured",
      "Lower held-out weighted log loss"
    ),
    result = c(
      sprintf(
        "Delta AIC %.1f; Delta BIC %.1f",
        future_test$delta_AIC_full_minus_reduced,
        future_test$delta_BIC_full_minus_reduced
      ),
      sprintf(
        "AIC %d/%d; BIC %d/%d folds",
        loyo_summary$aic_favours, loyo_summary$folds,
        loyo_summary$bic_favours, loyo_summary$folds
      ),
      sprintf(
        "%d/%d prefectures",
        lopo_summary$interaction_favours, lopo_summary$folds
      )
    ),
    p_value = c(future_test$p_value, NA_real_, NA_real_),
    observations = c(
      future_test$n, nobs(fit_spline_alt), nobs(fit_spline_alt)
    )
  )
)

model_qa <- bind_rows(
  tibble(
    model = c(
      "Primary absolute-temperature hierarchical",
      "ILI-footprint laboratory hierarchical",
      "Spring-2023-excluded hierarchical",
      "ILI hierarchical",
      "Pressure-adjusted reduced",
      "Pressure-adjusted interaction",
      "Future-anomaly reduced",
      "Future-anomaly interaction"
    ),
    convergence = c(
      primary_bundle$fits$positive_full$fit$convergence,
      fit_restricted$fit$convergence,
      primary_bundle$fits$positive_excluded$fit$convergence,
      fit_ili_hier$fit$convergence,
      fit_pressure_no_int$fit$convergence,
      fit_pressure_alt$fit$convergence,
      fit_future_no_int$fit$convergence,
      fit_future_alt$fit$convergence
    ),
    pdHess = c(
      primary_bundle$fits$positive_full$sdr$pdHess,
      fit_restricted$sdr$pdHess,
      primary_bundle$fits$positive_excluded$sdr$pdHess,
      fit_ili_hier$sdr$pdHess,
      fit_pressure_no_int$sdr$pdHess,
      fit_pressure_alt$sdr$pdHess,
      fit_future_no_int$sdr$pdHess,
      fit_future_alt$sdr$pdHess
    )
  ),
  lopo_results |>
    transmute(
      model = paste0("LOPO: ", omitted_prefecture),
      convergence = pmax(convergence_reduced, convergence_full),
      pdHess = pdHess_reduced & pdHess_full
    )
)

write_result5_csv(
  absolute_temperature_results,
  "result5_absolute_temperature_robustness.csv"
)
write_result5_csv(endpoint_results, "result5_endpoint_triangulation.csv")
write_result5_csv(pressure_test, "result5_pressure_anomaly_adjustment.csv")
write_result5_csv(future_test, "result5_future_anomaly_falsification.csv")
write_result5_csv(lopo_results, "result5_leave_one_prefecture_out.csv")
write_result5_csv(loyo_interaction_confirm, "result5_leave_one_year_out.csv")
write_result5_csv(phase_results, "result5_phase_specific_results.csv")
write_result5_csv(phase_tests, "result5_phase_heterogeneity_tests.csv")
write_result5_csv(main_table, "result5_main_table_source.csv")
write_result5_csv(model_qa, "result5_model_qa.csv")

if (any(model_qa$convergence != 0) || any(!model_qa$pdHess)) {
  warning("At least one Result 5 hierarchical fit failed convergence/Hessian QA.")
}

message("Completed updated 2-m air-temperature Result 5 analyses.")
