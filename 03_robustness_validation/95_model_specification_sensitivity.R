suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(glmmTMB)
  library(readr)
})

set.seed(20260902)

# Load the updated-air-temperature data preparation and primary model bundle.
# Existing primary fits are loaded from checkpoints by the sourced script.
source(file.path("02_primary_models", "80_primary_models.R"))

analysis_id <- "result5_model_specification_sensitivity_20260902"
checkpoint_dir_s6 <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", analysis_id
)
output_dir_s6 <- file.path("outputs", analysis_id)
source_dir_s6 <- file.path(output_dir_s6, "source_data")
dir.create(checkpoint_dir_s6, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir_s6, recursive = TRUE, showWarnings = FALSE)

primary_bundle_s6 <- readRDS(file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825", "v4_t2m_model_bundle.rds"
))

specifications <- tibble::tribble(
  ~specification_id, ~analysis_domain, ~specification_label, ~lag_max, ~temperature_df, ~knot_scheme, ~trend_control,
  "primary", "Primary", "Primary specification", 2L, 4L, "Default quantile knots", "Natural spline, 4 df",
  "lag01", "Lag window", "Cumulative lag 0-1 weeks", 1L, 4L, "Default quantile knots", "Natural spline, 4 df",
  "lag03", "Lag window", "Cumulative lag 0-3 weeks", 3L, 4L, "Default quantile knots", "Natural spline, 4 df",
  "temp_df3", "Temperature basis", "Temperature spline, 3 df", 2L, 3L, "Default quantile knots", "Natural spline, 4 df",
  "temp_df5", "Temperature basis", "Temperature spline, 5 df", 2L, 5L, "Default quantile knots", "Natural spline, 4 df",
  "tail_knots", "Knot placement", "Temperature knots at P10, P50 and P90", 2L, 4L, "P10/P50/P90 knots", "Natural spline, 4 df",
  "trend_df8", "Temporal control", "Long-term trend spline, 8 df", 2L, 4L, "Default quantile knots", "Natural spline, 8 df",
  "year_fixed", "Temporal control", "Calendar-year fixed effects", 2L, 4L, "Default quantile knots", "Calendar-year fixed effects"
)

control_s6 <- glmmTMBControl(optCtrl = list(iter.max = 5000, eval.max = 5000))
annual_period_weeks_s6 <- 365.2425 / 7

make_cumulative_basis <- function(cb, cb_names, temperatures, lag_max, centre = 0) {
  temperatures <- sort(unique(c(temperatures, centre)))
  original_names <- colnames(cb)
  zero_vcov <- matrix(
    0, nrow = length(original_names), ncol = length(original_names),
    dimnames = list(original_names, original_names)
  )
  out <- matrix(
    NA_real_, nrow = length(temperatures), ncol = length(original_names),
    dimnames = list(as.character(temperatures), cb_names)
  )
  for (j in seq_along(original_names)) {
    coefficient <- rep(0, length(original_names))
    names(coefficient) <- original_names
    coefficient[[j]] <- 1
    prediction <- crosspred(
      cb, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = temperatures, cen = centre, cumul = TRUE
    )
    out[, j] <- prediction$cumfit[, paste0("lag", lag_max)]
  }
  out
}

joint_wald <- function(fit, coefficient_names) {
  beta <- fixef(fit)$cond[coefficient_names]
  covariance <- vcov(fit)$cond[coefficient_names, coefficient_names, drop = FALSE]
  statistic <- as.numeric(crossprod(beta, solve(covariance, beta)))
  tibble(
    global_wald = statistic,
    global_df = length(coefficient_names),
    global_p = pchisq(statistic, df = length(coefficient_names), lower.tail = FALSE)
  )
}

altitude_modification_contrast <- function(
    fit, cb, cb_names, z_names, interaction_scale, lag_max,
    colder = -5, warmer = 0) {
  basis <- make_cumulative_basis(
    cb, cb_names, c(colder, warmer), lag_max = lag_max, centre = 0
  )
  basis_difference <- basis[as.character(colder), ] - basis[as.character(warmer), ]
  design <- basis_difference / interaction_scale
  beta <- fixef(fit)$cond[z_names]
  covariance <- vcov(fit)$cond[z_names, z_names, drop = FALSE]
  log_estimate <- sum(design * beta)
  standard_error <- sqrt(as.numeric(t(design) %*% covariance %*% design))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

prepare_and_fit_specification <- function(spec_row) {
  specification_id <- spec_row$specification_id[[1]]
  lag_max <- spec_row$lag_max[[1]]
  temperature_df <- spec_row$temperature_df[[1]]
  knot_scheme <- spec_row$knot_scheme[[1]]
  trend_control <- spec_row$trend_control[[1]]

  valid_temperature <- weekly_grid$temperature_for_basis[
    is.finite(weekly_grid$temperature_for_basis)
  ]
  argvar <- if (knot_scheme == "P10/P50/P90 knots") {
    list(
      fun = "ns",
      knots = as.numeric(quantile(valid_temperature, c(.10, .50, .90), na.rm = TRUE)),
      Boundary.knots = range(valid_temperature, na.rm = TRUE)
    )
  } else {
    list(fun = "ns", df = temperature_df)
  }

  cb <- crossbasis(
    weekly_grid$temperature_for_basis,
    lag = c(0, lag_max),
    argvar = argvar,
    arglag = list(fun = "integer"),
    group = weekly_grid$county
  )
  cb_df <- as.data.frame(cb)
  cb_names <- paste0("s6_cb_", seq_len(ncol(cb_df)))
  names(cb_df) <- cb_names

  dat <- bind_cols(
    weekly_grid |> mutate(exposure_row_id_s6 = row_number()), cb_df
  ) |>
    filter(
      lab_available %in% TRUE, tests >= 1, positive <= tests,
      complete.cases(across(all_of(cb_names))),
      !is.na(altitude_1000_c), !is.na(city)
    ) |>
    mutate(
      city = factor(city), county = factor(county),
      county_nested = interaction(city, county, drop = TRUE),
      pandemic_period = as.integer(iso_year %in% 2020:2022),
      calendar_year = factor(iso_year)
    ) |>
    arrange(county, week_index)

  season_names <- character()
  for (harmonic in 1:3) {
    sn <- paste0("s6_season_sin_", harmonic)
    cs <- paste0("s6_season_cos_", harmonic)
    dat[[sn]] <- sin(2 * pi * harmonic * dat$week_index / annual_period_weeks_s6)
    dat[[cs]] <- cos(2 * pi * harmonic * dat$week_index / annual_period_weeks_s6)
    season_names <- c(season_names, sn, cs)
  }

  temporal_names <- character()
  if (trend_control == "Calendar-year fixed effects") {
    temporal_names <- "calendar_year"
  } else {
    trend_df <- if (trend_control == "Natural spline, 8 df") 8L else 4L
    trend_basis <- splines::ns(
      dat$week_index, df = trend_df,
      Boundary.knots = range(dat$week_index)
    )
    trend_names <- paste0("s6_trend_", seq_len(ncol(trend_basis)))
    colnames(trend_basis) <- trend_names
    dat <- bind_cols(dat, as.data.frame(trend_basis))
    temporal_names <- c(trend_names, "pandemic_period")
  }

  raw_names <- paste0("s6_raw_int_", seq_along(cb_names))
  z_names <- paste0("s6_z_int_", seq_along(cb_names))
  for (i in seq_along(cb_names)) {
    dat[[raw_names[[i]]]] <- dat[[cb_names[[i]]]] * dat$altitude_1000_c
  }
  interaction_center <- vapply(dat[raw_names], mean, numeric(1), na.rm = TRUE)
  interaction_scale <- vapply(dat[raw_names], sd, numeric(1), na.rm = TRUE)
  if (any(!is.finite(interaction_scale) | interaction_scale <= 0)) {
    stop("Invalid interaction scaling for ", specification_id)
  }
  for (i in seq_along(raw_names)) {
    dat[[z_names[[i]]]] <-
      (dat[[raw_names[[i]]]] - interaction_center[[i]]) / interaction_scale[[i]]
  }

  local_support <- bind_rows(lapply(seq_len(nrow(primary_bundle_s6$altitude_anchors)), function(i) {
    altitude_value <- primary_bundle_s6$altitude_anchors$altitude[[i]]
    local_data <- dat |> filter(abs(altitude - altitude_value) <= 250)
    tibble(
      altitude_anchor = primary_bundle_s6$altitude_anchors$altitude_anchor[[i]],
      p10 = quantile(local_data$temperature, .10, na.rm = TRUE),
      p90 = quantile(local_data$temperature, .90, na.rm = TRUE)
    )
  }))
  common_lower <- max(local_support$p10)
  common_upper <- min(local_support$p90)
  contrast_supported <- -5 >= common_lower && 0 <= common_upper
  if (!contrast_supported) {
    stop("The prespecified -5 C versus 0 C contrast lacks common support for ", specification_id)
  }

  warmer_grid <- seq(
    ceiling((common_lower + 5) * 4) / 4,
    floor(common_upper * 4) / 4,
    by = 0.25
  )
  calibration_temperatures <- sort(unique(c(warmer_grid - 5, warmer_grid, 0)))
  calibration_basis <- make_cumulative_basis(
    cb, cb_names, calibration_temperatures, lag_max = lag_max, centre = 0
  )
  get_basis <- function(value) {
    calibration_basis[match(as.character(value), rownames(calibration_basis)), , drop = TRUE]
  }
  altitude_difference <-
    primary_bundle_s6$altitude_anchors$altitude_1000_c[
      primary_bundle_s6$altitude_anchors$altitude_anchor == "P90"
    ] -
    primary_bundle_s6$altitude_anchors$altitude_1000_c[
      primary_bundle_s6$altitude_anchors$altitude_anchor == "P10"
    ]
  contrast_norms <- vapply(seq_along(warmer_grid), function(i) {
    basis_difference <- get_basis(warmer_grid[[i]] - 5) - get_basis(warmer_grid[[i]])
    z_difference <- altitude_difference * basis_difference / interaction_scale
    sqrt(sum(z_difference^2))
  }, numeric(1))
  prior_sigma <- log(3) / (1.96 * max(contrast_norms))

  rhs <- c(
    cb_names, "altitude_1000_c", z_names, season_names, temporal_names,
    "(1 | city)", "(1 | county_nested)"
  )
  formula <- as.formula(paste(
    "cbind(positive, tests - positive) ~", paste(rhs, collapse = " + ")
  ))
  priors <- data.frame(
    prior = rep(sprintf("normal(0,%.12f)", prior_sigma), length(z_names)),
    class = rep("fixef", length(z_names)),
    coef = z_names,
    stringsAsFactors = FALSE
  )

  checkpoint <- file.path(checkpoint_dir_s6, paste0(specification_id, ".rds"))
  if (file.exists(checkpoint)) {
    message("Loading S6 checkpoint: ", specification_id)
    fit <- readRDS(checkpoint)
  } else {
    message("Fitting S6 specification: ", specification_id, " (n=", nrow(dat), ")")
    fit <- glmmTMB(
      formula, family = betabinomial(link = "logit"), data = dat,
      priors = priors, control = control_s6
    )
    saveRDS(fit, checkpoint)
  }

  wald <- joint_wald(fit, z_names)
  contrast <- altitude_modification_contrast(
    fit, cb, cb_names, z_names, interaction_scale, lag_max,
    colder = -5, warmer = 0
  )
  bind_cols(
    spec_row,
    tibble(
      observations = nobs(fit),
      tests = sum(dat$tests),
      interaction_terms = length(z_names),
      common_support_lower = common_lower,
      common_support_upper = common_upper,
      contrast_supported = contrast_supported,
      prior_sigma = prior_sigma,
      convergence = fit$fit$convergence,
      pdHess = fit$sdr$pdHess,
      logLik = as.numeric(logLik(fit)),
      AIC = AIC(fit),
      BIC = BIC(fit)
    ),
    wald,
    contrast
  )
}

primary_result <- local({
  fit <- primary_bundle_s6$fits$positive_full
  wald <- joint_wald(fit, primary_bundle_s6$z_int_names)
  contrast <- altitude_modification_contrast(
    fit = fit,
    cb = primary_bundle_s6$cb_all,
    cb_names = primary_bundle_s6$cb_names,
    z_names = primary_bundle_s6$z_int_names,
    interaction_scale = primary_bundle_s6$interaction_scale,
    lag_max = 2,
    colder = -5,
    warmer = 0
  )
  bind_cols(
    specifications |> filter(specification_id == "primary"),
    tibble(
      observations = nobs(fit),
      tests = sum(primary_bundle_s6$model_data$tests),
      interaction_terms = length(primary_bundle_s6$z_int_names),
      common_support_lower = max(primary_bundle_s6$local_support$p10),
      common_support_upper = min(primary_bundle_s6$local_support$p90),
      contrast_supported = TRUE,
      prior_sigma = primary_bundle_s6$prior_sigma,
      convergence = fit$fit$convergence,
      pdHess = fit$sdr$pdHess,
      logLik = as.numeric(logLik(fit)),
      AIC = AIC(fit),
      BIC = BIC(fit)
    ),
    wald,
    contrast
  )
})

alternative_results <- bind_rows(lapply(
  split(specifications |> filter(specification_id != "primary"),
        specifications$specification_id[specifications$specification_id != "primary"]),
  prepare_and_fit_specification
))

results <- bind_rows(primary_result, alternative_results) |>
  mutate(
    specification_id = factor(specification_id, levels = specifications$specification_id)
  ) |>
  arrange(specification_id) |>
  mutate(
    specification_id = as.character(specification_id),
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", estimate, low, high),
    global_p_display = case_when(
      global_p < .001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", global_p)
    ),
    direction_consistent = estimate >= 1,
    qa_pass = convergence == 0 & pdHess & contrast_supported &
      is.finite(estimate) & is.finite(low) & is.finite(high) & low > 0
  )

write_csv_via_temp <- function(data, destination) {
  temporary_file <- tempfile(fileext = ".csv")
  readr::write_csv(data, temporary_file)
  if (!file.copy(temporary_file, destination, overwrite = TRUE)) {
    stop("Unable to copy CSV output to: ", destination)
  }
}

write_csv_via_temp(
  results,
  file.path(source_dir_s6, "table_s6_model_specification_sensitivity.csv")
)
write_csv_via_temp(
  results |>
    select(
      specification_id, specification_label, observations, interaction_terms,
      convergence, pdHess, contrast_supported, qa_pass
    ),
  file.path(source_dir_s6, "table_s6_model_qa.csv")
)

if (!all(results$qa_pass)) {
  stop("At least one Supplementary Table S6 specification failed QA.")
}

message("Completed Supplementary Table S6 model-specification sensitivity analyses.")
