suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(glmmTMB)
})

set.seed(20260826)

analysis_id <- "result4_climate_context_phase_20260826"
model_root <- file.path("model_outputs", "reanalysis_v4_t2m_20260825")
panel_path <- file.path(model_root, "county_week_panel_v4_t2m_2016_2025.rds")
primary_bundle_path <- file.path(
  model_root, "primary_result2_result4_models_20260825",
  "v4_t2m_model_bundle.rds"
)
checkpoint_dir <- file.path(model_root, analysis_id)
source_dir <- file.path("result0715", "source_data", paste0("v4_", analysis_id))
qa_dir <- file.path(
  "result0715", "figure_qa", "v4_t2m_20260825", analysis_id
)
for (path in c(checkpoint_dir, source_dir, qa_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

panel <- readRDS(panel_path)
primary <- readRDS(primary_bundle_path)
anchors <- primary$altitude_anchors
altitude_mean <- primary$altitude_mean

# -------------------------------------------------------------------------
# 1. Empirical county-level air-temperature gradients across altitude.
# Each weekly slope compares counties under the same synoptic week. Monthly and
# seasonal summaries therefore describe the observed county-level gradient, not
# a universal free-atmosphere lapse rate.
# -------------------------------------------------------------------------
met <- panel |>
  filter(
    iso_year >= 2016, iso_year <= 2024,
    met_valid_5day_main %in% TRUE,
    is.finite(temperature), is.finite(altitude)
  ) |>
  mutate(
    altitude_km = altitude / 1000,
    month = month(week_start),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% 3:5 ~ "Spring",
      month %in% 6:8 ~ "Summer",
      TRUE ~ "Autumn"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  )

weekly_gradient <- bind_rows(lapply(
  split(met, met$week_start),
  function(dat) {
    if (n_distinct(dat$county) < 30L) return(NULL)
    fit <- lm(temperature ~ altitude_km, data = dat)
    coefficient <- summary(fit)$coefficients["altitude_km", ]
    tibble(
      week_start = unique(dat$week_start)[[1]],
      iso_year = unique(dat$iso_year)[[1]],
      iso_week = unique(dat$iso_week)[[1]],
      month = unique(dat$month)[[1]],
      season = as.character(unique(dat$season)[[1]]),
      counties = n_distinct(dat$county),
      slope_c_per_1000m = unname(coefficient[["Estimate"]]),
      slope_standard_error = unname(coefficient[["Std. Error"]]),
      r_squared = summary(fit)$r.squared
    )
  }
))

bootstrap_mean_ci <- function(x, replicates = 3000L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(mean = mean(x), low = NA, high = NA))
  means <- replicate(replicates, mean(sample(x, length(x), replace = TRUE)))
  c(
    mean = mean(x),
    low = unname(quantile(means, 0.025)),
    high = unname(quantile(means, 0.975))
  )
}

gradient_monthly <- bind_rows(lapply(sort(unique(weekly_gradient$month)), function(m) {
  z <- weekly_gradient |> filter(month == m)
  ci <- bootstrap_mean_ci(z$slope_c_per_1000m)
  tibble(
    period_type = "Month", period = month.abb[[m]], period_order = m,
    weeks = nrow(z), county_weeks = sum(z$counties),
    estimate = ci[["mean"]], low = ci[["low"]], high = ci[["high"]],
    median = median(z$slope_c_per_1000m),
    iqr_low = quantile(z$slope_c_per_1000m, .25),
    iqr_high = quantile(z$slope_c_per_1000m, .75)
  )
}))

season_levels <- c("Winter", "Spring", "Summer", "Autumn")
gradient_seasonal <- bind_rows(lapply(seq_along(season_levels), function(i) {
  label <- season_levels[[i]]
  z <- weekly_gradient |> filter(season == label)
  ci <- bootstrap_mean_ci(z$slope_c_per_1000m)
  tibble(
    period_type = "Season", period = label, period_order = i,
    weeks = nrow(z), county_weeks = sum(z$counties),
    estimate = ci[["mean"]], low = ci[["low"]], high = ci[["high"]],
    median = median(z$slope_c_per_1000m),
    iqr_low = quantile(z$slope_c_per_1000m, .25),
    iqr_high = quantile(z$slope_c_per_1000m, .75)
  )
}))

annual_ci <- bootstrap_mean_ci(weekly_gradient$slope_c_per_1000m)
gradient_annual <- tibble(
  period_type = "Annual", period = "Annual", period_order = 0L,
  weeks = nrow(weekly_gradient), county_weeks = sum(weekly_gradient$counties),
  estimate = annual_ci[["mean"]], low = annual_ci[["low"]],
  high = annual_ci[["high"]], median = median(weekly_gradient$slope_c_per_1000m),
  iqr_low = quantile(weekly_gradient$slope_c_per_1000m, .25),
  iqr_high = quantile(weekly_gradient$slope_c_per_1000m, .75)
)
gradient_summary <- bind_rows(gradient_annual, gradient_seasonal, gradient_monthly)

# -------------------------------------------------------------------------
# 2. Leave-one-year-out county/week climatology and temperature anomaly.
# The expected temperature for each county-week is estimated from the same ISO
# week in all other years from 2016-2024, preventing self-inclusion.
# -------------------------------------------------------------------------
climatology_rows <- lapply(2016:2024, function(target_year) {
  reference <- met |>
    filter(iso_year != target_year) |>
    group_by(county, iso_week) |>
    summarise(
      loyo_climatology = mean(temperature),
      climatology_years = n_distinct(iso_year),
      climatology_observations = n(),
      .groups = "drop"
    )
  met |>
    filter(iso_year == target_year) |>
    select(county, week_start, iso_year, iso_week, temperature) |>
    left_join(reference, by = c("county", "iso_week")) |>
    mutate(
      temperature_anomaly_loyo = temperature - loyo_climatology
    )
})
anomaly_weather <- bind_rows(climatology_rows)

week_min <- as.Date("2017-01-02")
week_max <- as.Date("2024-12-23")
county_attributes <- panel |>
  filter(!is.na(altitude)) |>
  distinct(county, city, altitude, altitude_1000_c)

weekly_grid <- tidyr::expand_grid(
  county = sort(unique(county_attributes$county)),
  week_start = seq(week_min, week_max, by = "7 days")
) |>
  mutate(
    iso_year = isoyear(week_start), iso_week = isoweek(week_start),
    week_id = sprintf("%04d-W%02d", iso_year, iso_week)
  ) |>
  left_join(
    panel |> filter(week_start >= week_min, week_start <= week_max),
    by = c("county", "week_start", "iso_year", "iso_week", "week_id"),
    suffix = c("", "_panel")
  ) |>
  left_join(
    anomaly_weather |>
      select(
        county, week_start, loyo_climatology,
        climatology_years, climatology_observations,
        temperature_anomaly_loyo
      ),
    by = c("county", "week_start")
  ) |>
  left_join(
    county_attributes |>
      select(
        county, city_attribute = city, altitude_attribute = altitude,
        altitude_1000_c_attribute = altitude_1000_c
      ),
    by = "county"
  ) |>
  mutate(
    city = coalesce(city, city_attribute),
    altitude = coalesce(altitude, altitude_attribute),
    altitude_1000_c = coalesce(altitude_1000_c, altitude_1000_c_attribute),
    anomaly_for_basis = if_else(
      met_valid_5day_main %in% TRUE,
      temperature_anomaly_loyo, NA_real_
    ),
    week_index = as.integer((week_start - week_min) / 7) + 1L
  ) |>
  select(
    -ends_with("_panel"), -city_attribute, -altitude_attribute,
    -altitude_1000_c_attribute
  ) |>
  arrange(county, week_start)

cb_anomaly <- crossbasis(
  weekly_grid$anomaly_for_basis,
  lag = c(0, 2),
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "integer"),
  group = weekly_grid$county
)
cb_anomaly_df <- as.data.frame(cb_anomaly)
cb_anomaly_names <- paste0("cb_anom_", seq_len(ncol(cb_anomaly_df)))
names(cb_anomaly_df) <- cb_anomaly_names

anomaly_model_dat <- bind_cols(
  weekly_grid |> mutate(exposure_row_id = row_number()),
  cb_anomaly_df
) |>
  filter(
    lab_available %in% TRUE, tests >= 1, positive <= tests,
    complete.cases(across(all_of(cb_anomaly_names))),
    is.finite(loyo_climatology), !is.na(altitude_1000_c), !is.na(city)
  ) |>
  mutate(
    city = factor(city), county = factor(county),
    county_nested = interaction(city, county, drop = TRUE)
  ) |>
  arrange(county, week_index)

climatology_basis <- splines::ns(
  anomaly_model_dat$loyo_climatology, df = 4,
  Boundary.knots = range(anomaly_model_dat$loyo_climatology)
)
climatology_names <- paste0("clim_", seq_len(ncol(climatology_basis)))
colnames(climatology_basis) <- climatology_names
anomaly_model_dat <- bind_cols(
  anomaly_model_dat, as.data.frame(climatology_basis)
)

annual_period_weeks <- 365.2425 / 7
season_names <- character()
for (harmonic in 1:3) {
  sn <- paste0("season_sin_", harmonic)
  cs <- paste0("season_cos_", harmonic)
  anomaly_model_dat[[sn]] <- sin(
    2 * pi * harmonic * anomaly_model_dat$week_index / annual_period_weeks
  )
  anomaly_model_dat[[cs]] <- cos(
    2 * pi * harmonic * anomaly_model_dat$week_index / annual_period_weeks
  )
  season_names <- c(season_names, sn, cs)
}
trend_basis <- splines::ns(
  anomaly_model_dat$week_index, df = 4,
  Boundary.knots = range(anomaly_model_dat$week_index)
)
trend_names <- paste0("portable_trend_", seq_len(ncol(trend_basis)))
colnames(trend_basis) <- trend_names
anomaly_model_dat <- bind_cols(
  anomaly_model_dat, as.data.frame(trend_basis)
) |>
  mutate(
    pandemic_period = as.integer(iso_year %in% 2020:2022),
    spring2023_event = as.integer(
      week_start <= as.Date("2023-05-31") &
        week_start + 6 >= as.Date("2023-02-01")
    )
  )

raw_anomaly_int_names <- paste0("raw_anom_int_", seq_along(cb_anomaly_names))
z_anomaly_int_names <- paste0("z_anom_int_", seq_along(cb_anomaly_names))
for (i in seq_along(cb_anomaly_names)) {
  anomaly_model_dat[[raw_anomaly_int_names[[i]]]] <-
    anomaly_model_dat[[cb_anomaly_names[[i]]]] *
    anomaly_model_dat$altitude_1000_c
}
anomaly_interaction_center <- vapply(
  anomaly_model_dat[raw_anomaly_int_names], mean, numeric(1), na.rm = TRUE
)
anomaly_interaction_scale <- vapply(
  anomaly_model_dat[raw_anomaly_int_names], sd, numeric(1), na.rm = TRUE
)
if (any(!is.finite(anomaly_interaction_scale) | anomaly_interaction_scale <= 0)) {
  stop("Invalid anomaly-interaction scaling.")
}
for (i in seq_along(raw_anomaly_int_names)) {
  anomaly_model_dat[[z_anomaly_int_names[[i]]]] <-
    (anomaly_model_dat[[raw_anomaly_int_names[[i]]]] -
       anomaly_interaction_center[[i]]) /
    anomaly_interaction_scale[[i]]
}

anomaly_support <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  a <- anchors$altitude[[i]]
  z <- anomaly_model_dat |> filter(abs(altitude - a) <= 250)
  tibble(
    altitude_anchor = anchors$altitude_anchor[[i]],
    altitude = a, n = nrow(z),
    p01 = quantile(z$temperature_anomaly_loyo, .01),
    p05 = quantile(z$temperature_anomaly_loyo, .05),
    p10 = quantile(z$temperature_anomaly_loyo, .10),
    median = median(z$temperature_anomaly_loyo),
    p90 = quantile(z$temperature_anomaly_loyo, .90),
    p95 = quantile(z$temperature_anomaly_loyo, .95),
    p99 = quantile(z$temperature_anomaly_loyo, .99)
  )
}))
common_anomaly_p10_lower <- max(anomaly_support$p10)
common_anomaly_p10_upper <- min(anomaly_support$p90)
common_anomaly_p05_lower <- max(anomaly_support$p05)
common_anomaly_p05_upper <- min(anomaly_support$p95)
if (!(0 >= common_anomaly_p10_lower && 0 <= common_anomaly_p10_upper)) {
  stop("Zero anomaly is outside common P10-P90 support.")
}

candidate_cold_anomalies <- c(-3, -2, -1)
valid_candidates <- candidate_cold_anomalies[
  candidate_cold_anomalies >= common_anomaly_p10_lower
]
if (length(valid_candidates) == 0L) {
  stop("No whole-degree cold-anomaly contrast lies in common P10-P90 support.")
}
primary_cold_anomaly <- min(valid_candidates)
primary_anomaly_width <- abs(primary_cold_anomaly)

anomaly_basis_cumulative <- function(values, centre = 0) {
  values <- sort(unique(c(values, centre)))
  original_names <- colnames(cb_anomaly)
  output <- matrix(
    NA_real_, length(values), length(original_names),
    dimnames = list(as.character(values), cb_anomaly_names)
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
      cb_anomaly, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = values, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

altitude_span_km <- diff(anchors$altitude_1000_c[c(1, 3)])
warmer_grid <- seq(
  ceiling((common_anomaly_p10_lower + primary_anomaly_width) * 4) / 4,
  floor(common_anomaly_p10_upper * 4) / 4,
  by = 0.25
)
anomaly_contrast_grid <- tibble(
  colder = warmer_grid - primary_anomaly_width,
  warmer = warmer_grid
)
calibration_values <- sort(unique(c(
  anomaly_contrast_grid$colder,
  anomaly_contrast_grid$warmer, 0
)))
calibration_basis <- anomaly_basis_cumulative(calibration_values, 0)
get_anomaly_basis <- function(value) {
  calibration_basis[
    match(as.character(value), rownames(calibration_basis)), , drop = TRUE
  ]
}
contrast_norms <- vapply(seq_len(nrow(anomaly_contrast_grid)), function(i) {
  basis_difference <-
    get_anomaly_basis(anomaly_contrast_grid$colder[[i]]) -
    get_anomaly_basis(anomaly_contrast_grid$warmer[[i]])
  z_difference <-
    altitude_span_km * basis_difference / anomaly_interaction_scale
  sqrt(sum(z_difference^2))
}, numeric(1))
prior_or_bound <- 3
anomaly_prior_sigma <-
  log(prior_or_bound) / (1.96 * max(contrast_norms))

portable_names <- c(
  season_names, trend_names, "pandemic_period", climatology_names
)
anomaly_base_rhs <- c(
  cb_anomaly_names, "altitude_1000_c", portable_names,
  "(1 | city)", "(1 | county_nested)"
)
anomaly_full_rhs <- c(
  cb_anomaly_names, "altitude_1000_c", z_anomaly_int_names,
  portable_names, "(1 | city)", "(1 | county_nested)"
)
anomaly_reduced_formula <- as.formula(paste(
  "cbind(positive, tests - positive) ~",
  paste(anomaly_base_rhs, collapse = " + ")
))
anomaly_full_formula <- as.formula(paste(
  "cbind(positive, tests - positive) ~",
  paste(anomaly_full_rhs, collapse = " + ")
))
anomaly_event_formula <- as.formula(paste(
  "cbind(positive, tests - positive) ~",
  paste(c(anomaly_full_rhs, "spring2023_event"), collapse = " + ")
))
anomaly_event_reduced_formula <- as.formula(paste(
  "cbind(positive, tests - positive) ~",
  paste(c(anomaly_base_rhs, "spring2023_event"), collapse = " + ")
))
anomaly_priors <- data.frame(
  prior = rep(
    sprintf("normal(0,%.12f)", anomaly_prior_sigma),
    length(z_anomaly_int_names)
  ),
  class = "fixef", coef = z_anomaly_int_names,
  stringsAsFactors = FALSE
)
control <- glmmTMBControl(optCtrl = list(iter.max = 4000, eval.max = 4000))

fit_or_load <- function(path, formula, priors = NULL, data) {
  if (file.exists(path)) return(readRDS(path))
  fit <- if (is.null(priors)) {
    glmmTMB(
      formula, family = betabinomial(link = "logit"),
      data = data, control = control
    )
  } else {
    glmmTMB(
      formula, family = betabinomial(link = "logit"),
      data = data, priors = priors, control = control
    )
  }
  saveRDS(fit, path)
  fit
}

anomaly_reduced_fit <- fit_or_load(
  file.path(checkpoint_dir, "temperature_anomaly_reduced.rds"),
  anomaly_reduced_formula, NULL, anomaly_model_dat
)
anomaly_full_fit <- fit_or_load(
  file.path(checkpoint_dir, "temperature_anomaly_full.rds"),
  anomaly_full_formula, anomaly_priors, anomaly_model_dat
)
anomaly_event_fit <- fit_or_load(
  file.path(checkpoint_dir, "temperature_anomaly_event_adjusted_full.rds"),
  anomaly_event_formula, anomaly_priors, anomaly_model_dat
)
anomaly_event_reduced_fit <- fit_or_load(
  file.path(checkpoint_dir, "temperature_anomaly_event_adjusted_reduced.rds"),
  anomaly_event_reduced_formula, NULL, anomaly_model_dat
)

extract_fixed <- function(fit) {
  list(
    beta = fixef(fit)$cond,
    vcov = as.matrix(vcov(fit)$cond)
  )
}

joint_wald <- function(fit, terms) {
  fixed <- extract_fixed(fit)
  keep <- intersect(terms, names(fixed$beta))
  beta <- fixed$beta[keep]
  variance <- fixed$vcov[keep, keep, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% solve(variance, beta))
  tibble(
    statistic = statistic, df = length(keep),
    p_value = pchisq(statistic, df = length(keep), lower.tail = FALSE)
  )
}
anomaly_global_interaction <- joint_wald(
  anomaly_full_fit, z_anomaly_int_names
) |>
  mutate(
    test = "Joint Wald test of anomaly-temperature by altitude interaction",
    reduced_AIC = AIC(anomaly_reduced_fit),
    full_AIC = AIC(anomaly_full_fit),
    delta_AIC_full_minus_reduced = AIC(anomaly_full_fit) - AIC(anomaly_reduced_fit),
    reduced_BIC = BIC(anomaly_reduced_fit),
    full_BIC = BIC(anomaly_full_fit),
    delta_BIC_full_minus_reduced = BIC(anomaly_full_fit) - BIC(anomaly_reduced_fit)
  )
anomaly_event_global_interaction <- joint_wald(
  anomaly_event_fit, z_anomaly_int_names
) |>
  mutate(
    test = paste(
      "Joint Wald test after adjustment for the February-May 2023 event"
    ),
    reduced_AIC = AIC(anomaly_event_reduced_fit),
    full_AIC = AIC(anomaly_event_fit),
    delta_AIC_full_minus_reduced =
      AIC(anomaly_event_fit) - AIC(anomaly_event_reduced_fit),
    reduced_BIC = BIC(anomaly_event_reduced_fit),
    full_BIC = BIC(anomaly_event_fit),
    delta_BIC_full_minus_reduced =
      BIC(anomaly_event_fit) - BIC(anomaly_event_reduced_fit)
  )

estimate_anomaly_basis <- function(
    basis, altitude_1000_c, fit = anomaly_full_fit) {
  fixed <- extract_fixed(fit)
  design <- matrix(
    0, nrow(basis), length(fixed$beta),
    dimnames = list(NULL, names(fixed$beta))
  )
  design[, cb_anomaly_names] <- basis
  design[, z_anomaly_int_names] <-
    sweep(basis, 2, anomaly_interaction_scale, "/") * altitude_1000_c
  log_estimate <- as.numeric(design %*% fixed$beta)
  standard_error <- sqrt(pmax(
    rowSums((design %*% fixed$vcov) * design), 0
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

estimate_anomaly_contrast <- function(
    colder, warmer, altitude_1000_c, interaction_only = FALSE,
    fit = anomaly_full_fit) {
  basis <- anomaly_basis_cumulative(c(colder, warmer), 0)
  basis_difference <-
    basis[as.character(colder), ] - basis[as.character(warmer), ]
  fixed <- extract_fixed(fit)
  design <- rep(0, length(fixed$beta))
  names(design) <- names(fixed$beta)
  if (!interaction_only) design[cb_anomaly_names] <- basis_difference
  design[z_anomaly_int_names] <-
    altitude_1000_c * basis_difference / anomaly_interaction_scale
  log_estimate <- sum(design * fixed$beta)
  standard_error <- sqrt(as.numeric(
    t(design) %*% fixed$vcov %*% design
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

anomaly_curves <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  spec <- anchors[i, ] |>
    left_join(anomaly_support, by = c("altitude_anchor", "altitude"))
  values <- sort(unique(c(
    seq(
      floor(spec$p05[[1]] * 4) / 4,
      ceiling(spec$p95[[1]] * 4) / 4,
      by = 0.25
    ),
    0
  )))
  basis <- anomaly_basis_cumulative(values, 0)
  bind_cols(
    tibble(
      altitude_anchor = spec$altitude_anchor[[1]],
      altitude = spec$altitude[[1]],
      altitude_1000_c = spec$altitude_1000_c[[1]],
      anomaly = as.numeric(rownames(basis)),
      reference_anomaly = 0,
      support_p05 = spec$p05[[1]], support_p95 = spec$p95[[1]]
    ),
    estimate_anomaly_basis(basis, spec$altitude_1000_c[[1]])
  )
}))
anomaly_event_curves <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  spec <- anchors[i, ] |>
    left_join(anomaly_support, by = c("altitude_anchor", "altitude"))
  values <- sort(unique(c(
    seq(
      floor(spec$p05[[1]] * 4) / 4,
      ceiling(spec$p95[[1]] * 4) / 4,
      by = 0.25
    ),
    0
  )))
  basis <- anomaly_basis_cumulative(values, 0)
  bind_cols(
    tibble(
      altitude_anchor = spec$altitude_anchor[[1]],
      altitude = spec$altitude[[1]],
      altitude_1000_c = spec$altitude_1000_c[[1]],
      anomaly = as.numeric(rownames(basis)),
      reference_anomaly = 0,
      support_p05 = spec$p05[[1]], support_p95 = spec$p95[[1]]
    ),
    estimate_anomaly_basis(
      basis, spec$altitude_1000_c[[1]], anomaly_event_fit
    )
  )
}))

anomaly_altitude_values <- anomaly_model_dat$altitude
anomaly_altitude_grid <- seq(
  floor(min(anomaly_altitude_values) / 25) * 25,
  ceiling(max(anomaly_altitude_values) / 25) * 25,
  by = 25
)
anomaly_fixed_curve <- bind_rows(lapply(anomaly_altitude_grid, function(a) {
  bind_cols(
    tibble(
      colder_anomaly = primary_cold_anomaly,
      warmer_anomaly = 0,
      altitude = a
    ),
    estimate_anomaly_contrast(
      primary_cold_anomaly, 0,
      (a - altitude_mean) / 1000,
      interaction_only = FALSE
    )
  )
}))
anomaly_event_fixed_curve <- bind_rows(lapply(
  anomaly_altitude_grid, function(a) {
    bind_cols(
      tibble(
        colder_anomaly = primary_cold_anomaly,
        warmer_anomaly = 0,
        altitude = a
      ),
      estimate_anomaly_contrast(
        primary_cold_anomaly, 0,
        (a - altitude_mean) / 1000,
        FALSE, anomaly_event_fit
      )
    )
  }
))
anomaly_anchor_contrasts <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  bind_cols(
    anchors[i, ],
    tibble(
      colder_anomaly = primary_cold_anomaly,
      warmer_anomaly = 0,
      effect_measure = "Odds ratio"
    ),
    estimate_anomaly_contrast(
      primary_cold_anomaly, 0, anchors$altitude_1000_c[[i]], FALSE
    )
  )
}))
anomaly_per1000m <- bind_cols(
  tibble(
    altitude_anchor = "Per 1,000 m",
    altitude = NA_real_, altitude_1000_c = 1,
    colder_anomaly = primary_cold_anomaly,
    warmer_anomaly = 0,
    effect_measure = "Ratio of odds ratios per 1,000 m"
  ),
  estimate_anomaly_contrast(
    primary_cold_anomaly, 0, 1, interaction_only = TRUE
  )
)
anomaly_key_contrasts <- bind_rows(
  anomaly_anchor_contrasts, anomaly_per1000m
)
anomaly_event_anchor_contrasts <- bind_rows(lapply(
  seq_len(nrow(anchors)), function(i) {
    bind_cols(
      anchors[i, ],
      tibble(
        colder_anomaly = primary_cold_anomaly,
        warmer_anomaly = 0,
        effect_measure = "Odds ratio"
      ),
      estimate_anomaly_contrast(
        primary_cold_anomaly, 0, anchors$altitude_1000_c[[i]],
        FALSE, anomaly_event_fit
      )
    )
  }
))
anomaly_event_per1000m <- bind_cols(
  tibble(
    altitude_anchor = "Per 1,000 m",
    altitude = NA_real_, altitude_1000_c = 1,
    colder_anomaly = primary_cold_anomaly,
    warmer_anomaly = 0,
    effect_measure = "Ratio of odds ratios per 1,000 m"
  ),
  estimate_anomaly_contrast(
    primary_cold_anomaly, 0, 1, TRUE, anomaly_event_fit
  )
)
anomaly_event_key_contrasts <- bind_rows(
  anomaly_event_anchor_contrasts, anomaly_event_per1000m
)

# -------------------------------------------------------------------------
# 3. Joint altitude-temperature environmental contrasts from the unchanged
# full-sample absolute-temperature models. These contrasts simultaneously move
# altitude by +1,000 m and temperature by the locally estimated gradient.
# -------------------------------------------------------------------------
main_cb <- primary$cb_all
main_cb_names <- primary$cb_names
main_z_names <- primary$z_int_names
main_interaction_scale <- primary$interaction_scale
main_reference <- primary$reference_temperature

main_basis_cumulative <- function(values, centre = main_reference) {
  values <- sort(unique(c(values, centre)))
  original_names <- colnames(main_cb)
  output <- matrix(
    NA_real_, length(values), length(original_names),
    dimnames = list(as.character(values), main_cb_names)
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
      main_cb, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = values, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

estimate_joint_environment <- function(
    fit, low_temperature, high_temperature,
    low_altitude, high_altitude) {
  basis <- main_basis_cumulative(
    c(low_temperature, high_temperature), main_reference
  )
  b_low <- basis[as.character(low_temperature), ]
  b_high <- basis[as.character(high_temperature), ]
  a_low <- (low_altitude - altitude_mean) / 1000
  a_high <- (high_altitude - altitude_mean) / 1000
  fixed <- extract_fixed(fit)
  design <- rep(0, length(fixed$beta))
  names(design) <- names(fixed$beta)
  design[main_cb_names] <- b_high - b_low
  design[["altitude_1000_c"]] <- a_high - a_low
  design[main_z_names] <-
    (b_high * a_high - b_low * a_low) / main_interaction_scale
  log_estimate <- sum(design * fixed$beta)
  standard_error <- sqrt(as.numeric(
    t(design) %*% fixed$vcov %*% design
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

gradient_for_joint <- bind_rows(
  gradient_annual |> select(period, estimate),
  gradient_seasonal |> select(period, estimate)
) |>
  rename(temperature_gradient_c_per_1000m = estimate)

typical_temperature <- bind_rows(lapply(gradient_for_joint$period, function(label) {
  z <- met |> filter(abs(altitude - anchors$altitude[[2]]) <= 250)
  if (label != "Annual") z <- z |> filter(season == label)
  tibble(
    period = label,
    low_altitude_typical_temperature = median(z$temperature, na.rm = TRUE)
  )
}))

joint_specs <- gradient_for_joint |>
  left_join(typical_temperature, by = "period") |>
  mutate(
    low_altitude = anchors$altitude[[2]],
    high_altitude = anchors$altitude[[2]] + 1000,
    high_altitude_expected_temperature =
      low_altitude_typical_temperature + temperature_gradient_c_per_1000m
  )

joint_fits <- list(
  `Overall positivity` = primary$fits$positive_full,
  `Influenza A-only positivity` = primary$fits$fluA_only_full,
  `Influenza B-only positivity` = primary$fits$fluB_only_full
)
joint_environment <- bind_rows(lapply(seq_len(nrow(joint_specs)), function(i) {
  bind_rows(lapply(names(joint_fits), function(outcome_name) {
    bind_cols(
      tibble(
        outcome = outcome_name,
        period = joint_specs$period[[i]],
        low_altitude = joint_specs$low_altitude[[i]],
        high_altitude = joint_specs$high_altitude[[i]],
        low_altitude_temperature =
          joint_specs$low_altitude_typical_temperature[[i]],
        high_altitude_temperature =
          joint_specs$high_altitude_expected_temperature[[i]],
        temperature_gradient_c_per_1000m =
          joint_specs$temperature_gradient_c_per_1000m[[i]]
      ),
      estimate_joint_environment(
        joint_fits[[outcome_name]],
        joint_specs$low_altitude_typical_temperature[[i]],
        joint_specs$high_altitude_expected_temperature[[i]],
        joint_specs$low_altitude[[i]], joint_specs$high_altitude[[i]]
      )
    )
  }))
})) |>
  mutate(
    low_temperature_supported =
      low_altitude_temperature >= primary$local_support$p05[[2]] &
      low_altitude_temperature <= primary$local_support$p95[[2]],
    high_temperature_supported =
      high_altitude_temperature >= primary$local_support$p05[[3]] &
      high_altitude_temperature <= primary$local_support$p95[[3]],
    both_temperatures_supported =
      low_temperature_supported & high_temperature_supported
  )

# -------------------------------------------------------------------------
# 4. Externally defined epidemic-period analysis. Separate period-specific fits
# use the same absolute-temperature basis, interaction scaling, prior and
# prespecified -5 versus 0 C contrast. The analysis is explicitly exploratory.
# -------------------------------------------------------------------------
phase_model_dat <- primary$model_data |>
  mutate(
    phase = case_when(
      calendar_year %in% 2017:2019 ~ "Pre-pandemic (2017-2019)",
      calendar_year %in% 2020:2022 ~ "Pandemic disruption (2020-2022)",
      calendar_year %in% 2023:2024 ~ "Post-transition (2023-2024)",
      TRUE ~ NA_character_
    ),
    phase = factor(
      phase,
      levels = c(
        "Pre-pandemic (2017-2019)",
        "Pandemic disruption (2020-2022)",
        "Post-transition (2023-2024)"
      )
    ),
    calendar_year_factor = factor(calendar_year)
  )

phase_support <- bind_rows(lapply(levels(phase_model_dat$phase), function(label) {
  bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
    a <- anchors$altitude[[i]]
    z <- phase_model_dat |>
      filter(phase == label, abs(altitude - a) <= 250)
    tibble(
      phase = label,
      altitude_anchor = anchors$altitude_anchor[[i]],
      altitude = a, n = nrow(z),
      p10 = quantile(z$temperature, .10),
      p90 = quantile(z$temperature, .90),
      contrast_supported = -5 >= quantile(z$temperature, .10) &
        0 <= quantile(z$temperature, .90)
    )
  }))
}))

phase_rhs <- c(
  main_cb_names, "altitude_1000_c", main_z_names,
  paste0("season_sin_", 1:3), paste0("season_cos_", 1:3),
  "calendar_year_factor", "(1 | city)", "(1 | county_nested)"
)
phase_formula <- as.formula(paste(
  "cbind(positive, tests - positive) ~",
  paste(phase_rhs, collapse = " + ")
))
phase_priors <- data.frame(
  prior = rep(
    sprintf("normal(0,%.12f)", primary$prior_sigma),
    length(main_z_names)
  ),
  class = "fixef", coef = main_z_names,
  stringsAsFactors = FALSE
)

safe_phase_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}
phase_fits <- list()
for (label in levels(phase_model_dat$phase)) {
  dat <- phase_model_dat |>
    filter(phase == label) |>
    droplevels()
  path <- file.path(
    checkpoint_dir, paste0("phase_", safe_phase_name(label), ".rds")
  )
  phase_fits[[label]] <- fit_or_load(
    path, phase_formula, phase_priors, dat
  )
}

estimate_main_contrast <- function(
    fit, colder, warmer, altitude_1000_c, interaction_only = FALSE) {
  basis <- main_basis_cumulative(c(colder, warmer), main_reference)
  basis_difference <-
    basis[as.character(colder), ] - basis[as.character(warmer), ]
  fixed <- extract_fixed(fit)
  design <- rep(0, length(fixed$beta))
  names(design) <- names(fixed$beta)
  if (!interaction_only) design[main_cb_names] <- basis_difference
  design[main_z_names] <-
    altitude_1000_c * basis_difference / main_interaction_scale
  log_estimate <- sum(design * fixed$beta)
  standard_error <- sqrt(as.numeric(
    t(design) %*% fixed$vcov %*% design
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

phase_results <- bind_rows(lapply(names(phase_fits), function(label) {
  fit <- phase_fits[[label]]
  global <- joint_wald(fit, main_z_names)
  bind_rows(
    bind_cols(
      tibble(
        phase = label,
        estimand = "-5 vs 0 C at 2,980 m",
        effect_measure = "Odds ratio",
        global_interaction_p = global$p_value,
        observations = nobs(fit),
        convergence = fit$fit$convergence,
        pdHess = fit$sdr$pdHess
      ),
      estimate_main_contrast(
        fit, -5, 0, anchors$altitude_1000_c[[2]], FALSE
      )
    ),
    bind_cols(
      tibble(
        phase = label,
        estimand = "Altitude modification of -5 vs 0 C",
        effect_measure = "Ratio of odds ratios per 1,000 m",
        global_interaction_p = global$p_value,
        observations = nobs(fit),
        convergence = fit$fit$convergence,
        pdHess = fit$sdr$pdHess
      ),
      estimate_main_contrast(fit, -5, 0, 1, TRUE)
    )
  )
}))

phase_candidate_results <- bind_rows(lapply(c(-3, -2), function(colder) {
  bind_rows(lapply(names(phase_fits), function(label) {
    fit <- phase_fits[[label]]
    global <- joint_wald(fit, main_z_names)
    bind_rows(
      bind_cols(
        tibble(
          phase = label,
          colder_temperature = colder,
          warmer_temperature = 0,
          estimand = sprintf("%s vs 0 C at 2,980 m", colder),
          effect_measure = "Odds ratio",
          global_interaction_p = global$p_value,
          observations = nobs(fit),
          convergence = fit$fit$convergence,
          pdHess = fit$sdr$pdHess
        ),
        estimate_main_contrast(
          fit, colder, 0, anchors$altitude_1000_c[[2]], FALSE
        )
      ),
      bind_cols(
        tibble(
          phase = label,
          colder_temperature = colder,
          warmer_temperature = 0,
          estimand = sprintf(
            "Altitude modification of %s vs 0 C", colder
          ),
          effect_measure = "Ratio of odds ratios per 1,000 m",
          global_interaction_p = global$p_value,
          observations = nobs(fit),
          convergence = fit$fit$convergence,
          pdHess = fit$sdr$pdHess
        ),
        estimate_main_contrast(fit, colder, 0, 1, TRUE)
      )
    )
  }))
}))

heterogeneity_test <- function(dat) {
  weights <- 1 / dat$standard_error^2
  pooled <- sum(weights * dat$log_estimate) / sum(weights)
  q <- sum(weights * (dat$log_estimate - pooled)^2)
  tibble(
    estimand = unique(dat$estimand),
    Q = q, df = nrow(dat) - 1,
    p_value = pchisq(q, df = nrow(dat) - 1, lower.tail = FALSE),
    fixed_effect_pooled_estimate = exp(pooled)
  )
}
phase_heterogeneity <- bind_rows(lapply(
  split(phase_results, phase_results$estimand),
  heterogeneity_test
))
phase_candidate_heterogeneity <- bind_rows(lapply(
  split(
    phase_candidate_results,
    interaction(
      phase_candidate_results$colder_temperature,
      phase_candidate_results$effect_measure,
      drop = TRUE
    )
  ),
  heterogeneity_test
))

# -------------------------------------------------------------------------
# Output, diagnostics and QA.
# -------------------------------------------------------------------------
anomaly_diagnostics <- tibble(
  model = c(
    "Anomaly reduced", "Anomaly interaction",
    "Anomaly reduced plus spring-2023 event",
    "Anomaly interaction plus spring-2023 event"
  ),
  observations = c(
    nobs(anomaly_reduced_fit), nobs(anomaly_full_fit),
    nobs(anomaly_event_reduced_fit), nobs(anomaly_event_fit)
  ),
  convergence = c(
    anomaly_reduced_fit$fit$convergence,
    anomaly_full_fit$fit$convergence,
    anomaly_event_reduced_fit$fit$convergence,
    anomaly_event_fit$fit$convergence
  ),
  pdHess = c(
    anomaly_reduced_fit$sdr$pdHess,
    anomaly_full_fit$sdr$pdHess,
    anomaly_event_reduced_fit$sdr$pdHess,
    anomaly_event_fit$sdr$pdHess
  ),
  AIC = c(
    AIC(anomaly_reduced_fit), AIC(anomaly_full_fit),
    AIC(anomaly_event_reduced_fit), AIC(anomaly_event_fit)
  ),
  BIC = c(
    BIC(anomaly_reduced_fit), BIC(anomaly_full_fit),
    BIC(anomaly_event_reduced_fit), BIC(anomaly_event_fit)
  )
)

qa <- tibble(
  check = c(
    "Meteorology spans 2016-2024 and all 45 counties",
    "Weekly gradient regressions retain at least 30 counties",
    "LOYO climatology excludes the target year",
    "Anomaly model retains at least 95% of primary model observations",
    "Zero anomaly lies within common P10-P90 support",
    "Primary cold anomaly was chosen only from exposure support",
    "Anomaly models converged with positive-definite Hessians",
    "All joint contrasts use supported P05-P95 temperatures",
    "All phase models converged with positive-definite Hessians",
    "-5 versus 0 C is supported in every phase/altitude anchor"
  ),
  pass = c(
    min(met$iso_year) == 2016 && max(met$iso_year) == 2024 &&
      n_distinct(met$county) == 45,
    min(weekly_gradient$counties) >= 30,
    min(anomaly_weather$climatology_years, na.rm = TRUE) >= 1,
    nobs(anomaly_full_fit) >= .95 * nobs(primary$fits$positive_full),
    0 >= common_anomaly_p10_lower && 0 <= common_anomaly_p10_upper,
    primary_cold_anomaly %in% candidate_cold_anomalies &&
      primary_cold_anomaly >= common_anomaly_p10_lower,
    all(anomaly_diagnostics$convergence == 0) &&
      all(anomaly_diagnostics$pdHess),
    all(joint_environment$both_temperatures_supported),
    all(phase_results$convergence == 0) && all(phase_results$pdHess),
    all(phase_support$contrast_supported)
  )
)

write_csv_utf8(weekly_gradient, file.path(source_dir, "weekly_temperature_altitude_gradients.csv"))
write_csv_utf8(gradient_summary, file.path(source_dir, "temperature_altitude_gradient_summary.csv"))
write_csv_utf8(anomaly_support, file.path(source_dir, "temperature_anomaly_support.csv"))
write_csv_utf8(anomaly_curves, file.path(source_dir, "temperature_anomaly_curves.csv"))
write_csv_utf8(anomaly_event_curves, file.path(source_dir, "temperature_anomaly_event_adjusted_curves.csv"))
write_csv_utf8(anomaly_fixed_curve, file.path(source_dir, "temperature_anomaly_fixed_contrast_across_altitude.csv"))
write_csv_utf8(anomaly_event_fixed_curve, file.path(source_dir, "temperature_anomaly_event_adjusted_fixed_contrast_across_altitude.csv"))
write_csv_utf8(anomaly_key_contrasts, file.path(source_dir, "temperature_anomaly_key_contrasts.csv"))
write_csv_utf8(anomaly_global_interaction, file.path(source_dir, "temperature_anomaly_global_interaction.csv"))
write_csv_utf8(anomaly_event_key_contrasts, file.path(source_dir, "temperature_anomaly_event_adjusted_key_contrasts.csv"))
write_csv_utf8(anomaly_event_global_interaction, file.path(source_dir, "temperature_anomaly_event_adjusted_global_interaction.csv"))
write_csv_utf8(anomaly_diagnostics, file.path(source_dir, "temperature_anomaly_model_diagnostics.csv"))
write_csv_utf8(joint_environment, file.path(source_dir, "joint_altitude_temperature_environmental_contrasts.csv"))
write_csv_utf8(phase_support, file.path(source_dir, "phase_temperature_support.csv"))
write_csv_utf8(phase_results, file.path(source_dir, "phase_specific_temperature_altitude_results.csv"))
write_csv_utf8(phase_heterogeneity, file.path(source_dir, "phase_heterogeneity_tests.csv"))
write_csv_utf8(phase_candidate_results, file.path(source_dir, "phase_candidate_contrast_results.csv"))
write_csv_utf8(phase_candidate_heterogeneity, file.path(source_dir, "phase_candidate_heterogeneity_tests.csv"))
write_csv_utf8(qa, file.path(qa_dir, "result4_climate_context_phase_qa.csv"))

saveRDS(list(
  weekly_gradient = weekly_gradient,
  gradient_summary = gradient_summary,
  anomaly_model_data = anomaly_model_dat,
  cb_anomaly = cb_anomaly,
  cb_anomaly_names = cb_anomaly_names,
  z_anomaly_int_names = z_anomaly_int_names,
  anomaly_interaction_scale = anomaly_interaction_scale,
  anomaly_support = anomaly_support,
  primary_cold_anomaly = primary_cold_anomaly,
  anomaly_prior_sigma = anomaly_prior_sigma,
  anomaly_reduced_fit = anomaly_reduced_fit,
  anomaly_full_fit = anomaly_full_fit,
  anomaly_event_fit = anomaly_event_fit,
  anomaly_event_reduced_fit = anomaly_event_reduced_fit,
  anomaly_curves = anomaly_curves,
  anomaly_event_curves = anomaly_event_curves,
  anomaly_fixed_curve = anomaly_fixed_curve,
  anomaly_event_fixed_curve = anomaly_event_fixed_curve,
  anomaly_key_contrasts = anomaly_key_contrasts,
  anomaly_event_key_contrasts = anomaly_event_key_contrasts,
  anomaly_global_interaction = anomaly_global_interaction,
  anomaly_event_global_interaction = anomaly_event_global_interaction,
  joint_environment = joint_environment,
  phase_fits = phase_fits,
  phase_results = phase_results,
  phase_heterogeneity = phase_heterogeneity,
  phase_candidate_results = phase_candidate_results,
  phase_candidate_heterogeneity = phase_candidate_heterogeneity
), file.path(checkpoint_dir, "result4_climate_context_phase_bundle.rds"))

if (!all(qa$pass)) {
  warning("At least one Result 4 QA item requires review; outputs were retained.")
}

message("Completed Result 4 climate-context and epidemic-period analyses.")
