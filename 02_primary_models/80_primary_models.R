suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(glmmTMB)
})

set.seed(20260825)

version_id <- "v4_t2m_20260825"
analysis_id <- "primary_result2_result4_models_20260825"
model_root <- file.path("model_outputs", paste0("reanalysis_", version_id))
panel_path <- file.path(model_root, "county_week_panel_v4_t2m_2016_2025.rds")
checkpoint_dir <- file.path(model_root, analysis_id)
analysis_dir <- file.path("05_frequentist_dlnm", paste0("reanalysis_", version_id), analysis_id)
source_dir <- file.path("result0715", "source_data", paste0("v4_", analysis_id))
for (z in c(checkpoint_dir, analysis_dir, source_dir)) {
  dir.create(z, recursive = TRUE, showWarnings = FALSE)
}
if (!file.exists(panel_path)) stop("Missing v4 T2m panel.")

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

panel <- readRDS(panel_path)
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
    county_attributes |>
      select(county, city_attribute = city, altitude_attribute = altitude,
             altitude_1000_c_attribute = altitude_1000_c),
    by = "county"
  ) |>
  mutate(
    city = coalesce(city, city_attribute),
    altitude = coalesce(altitude, altitude_attribute),
    altitude_1000_c = coalesce(altitude_1000_c, altitude_1000_c_attribute),
    temperature_for_basis = if_else(met_valid_5day_main %in% TRUE,
                                    temperature, NA_real_),
    week_index = as.integer((week_start - week_min) / 7) + 1L
  ) |>
  select(-ends_with("_panel"), -city_attribute, -altitude_attribute,
         -altitude_1000_c_attribute) |>
  arrange(county, week_start)

cb_all <- crossbasis(
  weekly_grid$temperature_for_basis,
  lag = c(0, 2),
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "integer"),
  group = weekly_grid$county
)
cb_df <- as.data.frame(cb_all)
cb_names <- paste0("cb_", seq_len(ncol(cb_df)))
names(cb_df) <- cb_names

model_dat <- bind_cols(
  weekly_grid |> mutate(exposure_row_id = row_number()), cb_df
) |>
  filter(
    lab_available %in% TRUE, tests >= 1, positive <= tests,
    complete.cases(across(all_of(cb_names))),
    !is.na(altitude_1000_c), !is.na(city)
  ) |>
  mutate(
    city = factor(city), county = factor(county),
    county_nested = interaction(city, county, drop = TRUE)
  ) |>
  arrange(county, week_index)

annual_period_weeks <- 365.2425 / 7
season_names <- character()
for (harmonic in 1:3) {
  sn <- paste0("season_sin_", harmonic)
  cs <- paste0("season_cos_", harmonic)
  model_dat[[sn]] <- sin(2 * pi * harmonic * model_dat$week_index / annual_period_weeks)
  model_dat[[cs]] <- cos(2 * pi * harmonic * model_dat$week_index / annual_period_weeks)
  season_names <- c(season_names, sn, cs)
}
trend_basis <- splines::ns(
  model_dat$week_index, df = 4,
  Boundary.knots = range(model_dat$week_index)
)
trend_names <- paste0("portable_trend_", seq_len(ncol(trend_basis)))
colnames(trend_basis) <- trend_names
model_dat <- bind_cols(model_dat, as.data.frame(trend_basis)) |>
  mutate(
    pandemic_period = as.integer(iso_year %in% 2020:2022),
    spring2023_event = as.integer(
      week_start <= as.Date("2023-05-31") &
        week_start + 6 >= as.Date("2023-02-01")
    )
  )

raw_int_names <- paste0("raw_int_", seq_along(cb_names))
z_int_names <- paste0("z_int_", seq_along(cb_names))
for (i in seq_along(cb_names)) {
  model_dat[[raw_int_names[[i]]]] <-
    model_dat[[cb_names[[i]]]] * model_dat$altitude_1000_c
}
interaction_center <- vapply(model_dat[raw_int_names], mean, numeric(1), na.rm = TRUE)
interaction_scale <- vapply(model_dat[raw_int_names], sd, numeric(1), na.rm = TRUE)
if (any(!is.finite(interaction_scale) | interaction_scale <= 0)) {
  stop("Invalid interaction scale.")
}
for (i in seq_along(raw_int_names)) {
  model_dat[[z_int_names[[i]]]] <-
    (model_dat[[raw_int_names[[i]]]] - interaction_center[[i]]) /
    interaction_scale[[i]]
}

altitude_values <- county_attributes |> distinct(county, altitude) |> pull(altitude)
altitude_mean <- mean(altitude_values)
altitude_anchor_values <- as.numeric(quantile(altitude_values, c(.1, .5, .9)))
altitude_anchors <- tibble(
  altitude_anchor = c("P10", "P50", "P90"),
  altitude = altitude_anchor_values,
  altitude_1000_c = (altitude_anchor_values - altitude_mean) / 1000
)
altitude_difference <- altitude_anchors$altitude_1000_c[3] -
  altitude_anchors$altitude_1000_c[1]

local_support <- bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
  a <- altitude_anchors$altitude[[i]]
  z <- model_dat |> filter(abs(altitude - a) <= 250)
  tibble(
    altitude_anchor = altitude_anchors$altitude_anchor[[i]],
    altitude = a, n = nrow(z),
    p01 = quantile(z$temperature, .01), p05 = quantile(z$temperature, .05),
    p10 = quantile(z$temperature, .10), median = median(z$temperature),
    p90 = quantile(z$temperature, .90), p95 = quantile(z$temperature, .95),
    p99 = quantile(z$temperature, .99)
  )
}))

reference_temperature <- 0
common_p10_lower <- max(local_support$p10)
common_p10_upper <- min(local_support$p90)
if (!(reference_temperature >= common_p10_lower &&
      reference_temperature <= common_p10_upper)) {
  stop("0 C is not within common local P10-P90 support.")
}

basis_matrix_cumulative <- function(temperatures, reference_temperature = 0) {
  temperatures <- sort(unique(c(temperatures, reference_temperature)))
  original_names <- colnames(cb_all)
  output <- matrix(
    NA_real_, nrow = length(temperatures), ncol = length(original_names),
    dimnames = list(as.character(temperatures), original_names)
  )
  zero_vcov <- matrix(
    0, nrow = length(original_names), ncol = length(original_names),
    dimnames = list(original_names, original_names)
  )
  for (j in seq_along(original_names)) {
    coefficient <- rep(0, length(original_names))
    names(coefficient) <- original_names
    coefficient[[j]] <- 1
    prediction <- crosspred(
      cb_all, coef = coefficient, vcov = zero_vcov, model.link = "logit",
      at = temperatures, cen = reference_temperature, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  colnames(output) <- cb_names
  output
}

# Calibrate the interaction prior using all 5 C contrasts within the shared
# P10-P90 support of the three altitude slices.
warmer_grid <- seq(
  ceiling((common_p10_lower + 5) * 4) / 4,
  floor(common_p10_upper * 4) / 4,
  by = 0.25
)
contrast_grid <- tibble(
  colder_temperature = warmer_grid - 5,
  warmer_temperature = warmer_grid
)
calibration_temperatures <- sort(unique(c(
  contrast_grid$colder_temperature,
  contrast_grid$warmer_temperature,
  reference_temperature
)))
calibration_basis <- basis_matrix_cumulative(
  calibration_temperatures, reference_temperature
)
get_basis <- function(temp) {
  calibration_basis[match(as.character(temp), rownames(calibration_basis)), , drop = TRUE]
}
contrast_norms <- vapply(seq_len(nrow(contrast_grid)), function(i) {
  basis_difference <- get_basis(contrast_grid$colder_temperature[[i]]) -
    get_basis(contrast_grid$warmer_temperature[[i]])
  z_difference <- altitude_difference * basis_difference / interaction_scale
  sqrt(sum(z_difference^2))
}, numeric(1))
prior_or_bound <- 3
prior_sigma <- log(prior_or_bound) / (1.96 * max(contrast_norms))

portable_names <- c(season_names, trend_names, "pandemic_period")
base_rhs <- c(
  cb_names, "altitude_1000_c", z_int_names, portable_names,
  "(1 | city)", "(1 | county_nested)"
)
make_formula <- function(numerator, event_adjusted = FALSE) {
  rhs <- if (event_adjusted) c(base_rhs, "spring2023_event") else base_rhs
  as.formula(paste(
    sprintf("cbind(%s, tests - %s) ~", numerator, numerator),
    paste(rhs, collapse = " + ")
  ))
}
make_priors <- function() {
  data.frame(
    prior = rep(sprintf("normal(0,%.12f)", prior_sigma), length(z_int_names)),
    class = rep("fixef", length(z_int_names)), coef = z_int_names,
    stringsAsFactors = FALSE
  )
}
control <- glmmTMBControl(optCtrl = list(iter.max = 4000, eval.max = 4000))

fit_or_load <- function(numerator, variant = "full") {
  path <- file.path(checkpoint_dir, paste0(numerator, "__", variant, ".rds"))
  if (file.exists(path)) {
    message("Loading checkpoint: ", basename(path))
    return(readRDS(path))
  }
  fit_data <- model_dat
  event_adjusted <- FALSE
  if (variant == "excluded") {
    fit_data <- model_dat |> filter(spring2023_event == 0) |> droplevels()
  } else if (variant == "event_adjusted") {
    event_adjusted <- TRUE
  }
  message("Fitting ", numerator, " / ", variant, " (n=", nrow(fit_data), ")")
  fit <- glmmTMB(
    make_formula(numerator, event_adjusted),
    family = betabinomial(link = "logit"), data = fit_data,
    priors = make_priors(), control = control
  )
  saveRDS(fit, path)
  fit
}

fits <- list(
  positive_full = fit_or_load("positive", "full"),
  fluA_only_full = fit_or_load("fluA_only", "full"),
  fluB_only_full = fit_or_load("fluB_only", "full"),
  positive_event_adjusted = fit_or_load("positive", "event_adjusted"),
  fluA_only_event_adjusted = fit_or_load("fluA_only", "event_adjusted"),
  positive_excluded = fit_or_load("positive", "excluded"),
  fluA_only_excluded = fit_or_load("fluA_only", "excluded"),
  fluB_only_excluded = fit_or_load("fluB_only", "excluded")
)

diagnostics <- bind_rows(lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  tibble(
    model = nm, observations = nobs(fit), convergence = fit$fit$convergence,
    pdHess = fit$sdr$pdHess, logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit), BIC = BIC(fit)
  )
}))
if (any(diagnostics$convergence != 0) || any(!diagnostics$pdHess)) {
  stop("At least one v4 model failed convergence/Hessian QA.")
}

extract_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, vcov = as.matrix(vcov(fit)$cond))
}
estimate_basis <- function(fit, basis, altitude_1000_c) {
  fixed <- extract_fixed(fit)
  design <- matrix(0, nrow(basis), length(fixed$beta),
                   dimnames = list(NULL, names(fixed$beta)))
  design[, cb_names] <- basis
  design[, z_int_names] <- sweep(basis, 2, interaction_scale, "/") * altitude_1000_c
  log_estimate <- as.numeric(design %*% fixed$beta)
  se <- sqrt(pmax(rowSums((design %*% fixed$vcov) * design), 0))
  tibble(
    estimate = exp(log_estimate), low = exp(log_estimate - 1.96 * se),
    high = exp(log_estimate + 1.96 * se), log_estimate = log_estimate,
    standard_error = se
  )
}
estimate_contrast <- function(fit, colder, warmer, altitude_1000_c,
                              interaction_only = FALSE) {
  b <- basis_matrix_cumulative(c(colder, warmer), reference_temperature)
  bd <- b[as.character(colder), ] - b[as.character(warmer), ]
  fixed <- extract_fixed(fit)
  design <- rep(0, length(fixed$beta)); names(design) <- names(fixed$beta)
  if (!interaction_only) design[cb_names] <- bd
  design[z_int_names] <- altitude_1000_c * bd / interaction_scale
  log_estimate <- sum(design * fixed$beta)
  se <- sqrt(as.numeric(t(design) %*% fixed$vcov %*% design))
  tibble(
    estimate = exp(log_estimate), low = exp(log_estimate - 1.96 * se),
    high = exp(log_estimate + 1.96 * se), log_estimate = log_estimate,
    standard_error = se
  )
}

# Result 2: overall curve at three altitude anchors.
result2_curves <- bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
  spec <- altitude_anchors[i, ] |> left_join(local_support, by = c("altitude_anchor", "altitude"))
  temperatures <- sort(unique(c(
    seq(floor(spec$p05[[1]] * 4) / 4,
        ceiling(spec$p95[[1]] * 4) / 4, by = 0.25),
    reference_temperature
  )))
  basis <- basis_matrix_cumulative(temperatures, reference_temperature)
  bind_cols(
    tibble(
      altitude_anchor = spec$altitude_anchor[[1]], altitude = spec$altitude[[1]],
      altitude_1000_c = spec$altitude_1000_c[[1]],
      reference_temperature = reference_temperature,
      support_p05 = spec$p05[[1]], support_p95 = spec$p95[[1]],
      temperature = as.numeric(rownames(basis))
    ),
    estimate_basis(fits$positive_full, basis, spec$altitude_1000_c[[1]])
  )
}))

moving_modification <- bind_rows(lapply(seq_len(nrow(contrast_grid)), function(i) {
  bind_cols(
    contrast_grid[i, ],
    estimate_contrast(
      fits$positive_full,
      contrast_grid$colder_temperature[[i]], contrast_grid$warmer_temperature[[i]],
      altitude_difference, interaction_only = TRUE
    )
  )
}))

fixed_contrasts <- tibble(
  contrast = c("-5 vs 0 C", "0 vs 5 C"),
  colder_temperature = c(-5, 0), warmer_temperature = c(0, 5)
)
altitude_grid <- seq(min(altitude_values), max(altitude_values), by = 25)
fixed_altitude_curves <- bind_rows(lapply(seq_len(nrow(fixed_contrasts)), function(i) {
  bind_rows(lapply(altitude_grid, function(a) {
    bind_cols(
      fixed_contrasts[i, ], tibble(altitude = a),
      estimate_contrast(
        fits$positive_full,
        fixed_contrasts$colder_temperature[[i]],
        fixed_contrasts$warmer_temperature[[i]],
        (a - altitude_mean) / 1000
      )
    )
  }))
}))

primary_cold_contrast <- bind_rows(
  bind_cols(
    tibble(altitude_anchor = "P10", altitude = altitude_anchors$altitude[[1]]),
    estimate_contrast(fits$positive_full, -5, 0,
                      altitude_anchors$altitude_1000_c[[1]])
  ),
  bind_cols(
    tibble(altitude_anchor = "P50", altitude = altitude_anchors$altitude[[2]]),
    estimate_contrast(fits$positive_full, -5, 0,
                      altitude_anchors$altitude_1000_c[[2]])
  ),
  bind_cols(
    tibble(altitude_anchor = "P90", altitude = altitude_anchors$altitude[[3]]),
    estimate_contrast(fits$positive_full, -5, 0,
                      altitude_anchors$altitude_1000_c[[3]])
  ),
  bind_cols(
    tibble(altitude_anchor = "P90_vs_P10", altitude = NA_real_),
    estimate_contrast(fits$positive_full, -5, 0, altitude_difference, TRUE)
  )
)

# Result 4: excluded curves for overall, A-only and B-only.
excluded_specs <- list(
  list(outcome = "Overall positivity", fit = fits$positive_excluded),
  list(outcome = "Influenza A-only positivity", fit = fits$fluA_only_excluded),
  list(outcome = "Influenza B-only positivity", fit = fits$fluB_only_excluded)
)
result4_excluded_curves <- bind_rows(lapply(excluded_specs, function(outcome_spec) {
  bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
    spec <- altitude_anchors[i, ] |> left_join(local_support, by = c("altitude_anchor", "altitude"))
    temperatures <- sort(unique(c(
      seq(floor(spec$p05[[1]] * 4) / 4,
          ceiling(spec$p95[[1]] * 4) / 4, by = 0.25),
      reference_temperature
    )))
    basis <- basis_matrix_cumulative(temperatures, reference_temperature)
    bind_cols(
      tibble(
        outcome = outcome_spec$outcome,
        altitude_anchor = spec$altitude_anchor[[1]], altitude = spec$altitude[[1]],
        altitude_1000_c = spec$altitude_1000_c[[1]],
        reference_temperature = reference_temperature,
        support_p05 = spec$p05[[1]], support_p95 = spec$p95[[1]],
        temperature = as.numeric(rownames(basis))
      ),
      estimate_basis(outcome_spec$fit, basis, spec$altitude_1000_c[[1]])
    )
  }))
}))

# Directly comparable all-years curves using the identical reference,
# local-support grids and axis-ready source structure as the excluded figure.
full_specs <- list(
  list(outcome = "Overall positivity", fit = fits$positive_full),
  list(outcome = "Influenza A-only positivity", fit = fits$fluA_only_full),
  list(outcome = "Influenza B-only positivity", fit = fits$fluB_only_full)
)
result4_all_years_curves <- bind_rows(lapply(full_specs, function(outcome_spec) {
  bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
    spec <- altitude_anchors[i, ] |> left_join(local_support, by = c("altitude_anchor", "altitude"))
    temperatures <- sort(unique(c(
      seq(floor(spec$p05[[1]] * 4) / 4,
          ceiling(spec$p95[[1]] * 4) / 4, by = 0.25),
      reference_temperature
    )))
    basis <- basis_matrix_cumulative(temperatures, reference_temperature)
    bind_cols(
      tibble(
        outcome = outcome_spec$outcome,
        altitude_anchor = spec$altitude_anchor[[1]], altitude = spec$altitude[[1]],
        altitude_1000_c = spec$altitude_1000_c[[1]],
        reference_temperature = reference_temperature,
        support_p05 = spec$p05[[1]], support_p95 = spec$p95[[1]],
        temperature = as.numeric(rownames(basis))
      ),
      estimate_basis(outcome_spec$fit, basis, spec$altitude_1000_c[[1]])
    )
  }))
}))

# Result 4 event evidence: baseline, event-adjusted and excluded curves at P10.
event_curve_specs <- list(
  list(outcome = "Overall positivity", method = "All years", fit = fits$positive_full),
  list(outcome = "Overall positivity", method = "Event-adjusted", fit = fits$positive_event_adjusted),
  list(outcome = "Overall positivity", method = "Excluding spring 2023", fit = fits$positive_excluded),
  list(outcome = "Influenza A-only positivity", method = "All years", fit = fits$fluA_only_full),
  list(outcome = "Influenza A-only positivity", method = "Event-adjusted", fit = fits$fluA_only_event_adjusted),
  list(outcome = "Influenza A-only positivity", method = "Excluding spring 2023", fit = fits$fluA_only_excluded)
)
low_spec <- altitude_anchors[1, ] |> left_join(local_support, by = c("altitude_anchor", "altitude"))
event_temperatures <- sort(unique(c(
  seq(floor(low_spec$p05[[1]] * 4) / 4,
      ceiling(low_spec$p95[[1]] * 4) / 4, by = 0.25),
  reference_temperature
)))
event_basis <- basis_matrix_cumulative(event_temperatures, reference_temperature)
result4_event_curves <- bind_rows(lapply(event_curve_specs, function(z) {
  bind_cols(
    tibble(
      outcome = z$outcome, method = z$method,
      altitude_anchor = "P10", altitude = low_spec$altitude[[1]],
      reference_temperature = reference_temperature,
      temperature = as.numeric(rownames(event_basis))
    ),
    estimate_basis(z$fit, event_basis, low_spec$altitude_1000_c[[1]])
  )
}))

event_indicators <- bind_rows(lapply(c("positive", "fluA_only"), function(numerator) {
  baseline <- fits[[paste0(numerator, "_full")]]
  adjusted <- fits[[paste0(numerator, "_event_adjusted")]]
  fixed <- extract_fixed(adjusted)
  term <- grep("^spring2023_event", names(fixed$beta), value = TRUE)[1]
  b <- fixed$beta[[term]]; se <- sqrt(fixed$vcov[term, term])
  ll0 <- logLik(baseline); ll1 <- logLik(adjusted)
  lr <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  df <- attr(ll1, "df") - attr(ll0, "df")
  tibble(
    outcome = ifelse(numerator == "positive", "Overall positivity",
                     "Influenza A-only positivity"),
    estimate = exp(b), low = exp(b - 1.96 * se), high = exp(b + 1.96 * se),
    baseline_AIC = AIC(baseline), event_adjusted_AIC = AIC(adjusted),
    delta_AIC = AIC(adjusted) - AIC(baseline), LR = lr, df = df,
    p_value = pchisq(lr, df = df, lower.tail = FALSE)
  )
}))

exclusion_summary <- model_dat |>
  summarise(
    full_county_weeks = n(), excluded_county_weeks = sum(spring2023_event == 1),
    retained_county_weeks = sum(spring2023_event == 0),
    full_tests = sum(tests), excluded_tests = sum(tests[spring2023_event == 1]),
    retained_tests = sum(tests[spring2023_event == 0]),
    full_A_only = sum(fluA_only), excluded_A_only = sum(fluA_only[spring2023_event == 1]),
    retained_A_only = sum(fluA_only[spring2023_event == 0])
  )

prior_calibration <- tibble(
  prior_or_bound = prior_or_bound, prior_sigma = prior_sigma,
  maximum_contrast_design_norm = max(contrast_norms),
  common_p10_lower = common_p10_lower, common_p10_upper = common_p10_upper,
  reference_temperature = reference_temperature
)

saveRDS(list(
  model_data = model_dat, cb_all = cb_all, cb_names = cb_names,
  z_int_names = z_int_names, raw_int_names = raw_int_names,
  interaction_center = interaction_center, interaction_scale = interaction_scale,
  altitude_anchors = altitude_anchors, altitude_mean = altitude_mean,
  local_support = local_support, prior_sigma = prior_sigma,
  reference_temperature = reference_temperature, fits = fits
), file.path(checkpoint_dir, "v4_t2m_model_bundle.rds"))

write_csv_utf8(diagnostics, file.path(source_dir, "model_diagnostics.csv"))
write_csv_utf8(local_support, file.path(source_dir, "local_temperature_support.csv"))
write_csv_utf8(prior_calibration, file.path(source_dir, "prior_calibration.csv"))
write_csv_utf8(result2_curves, file.path(source_dir, "result2_anchor_curves.csv"))
write_csv_utf8(moving_modification, file.path(source_dir, "result2_moving_altitude_modification.csv"))
write_csv_utf8(fixed_altitude_curves, file.path(source_dir, "result2_fixed_contrasts_across_altitude.csv"))
write_csv_utf8(primary_cold_contrast, file.path(source_dir, "result2_primary_cold_contrast.csv"))
write_csv_utf8(result4_excluded_curves, file.path(source_dir, "result4_excluded_curves.csv"))
write_csv_utf8(result4_all_years_curves, file.path(source_dir, "result4_all_years_curves.csv"))
write_csv_utf8(result4_event_curves, file.path(source_dir, "result4_event_curves.csv"))
write_csv_utf8(event_indicators, file.path(source_dir, "result4_event_indicator_estimates.csv"))
write_csv_utf8(exclusion_summary, file.path(source_dir, "result4_exclusion_summary.csv"))

message("Completed v4 T2m Result 2 / Result 4 model set.")
