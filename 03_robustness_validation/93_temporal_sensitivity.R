suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(20260902)

# Result 5 supplementary analysis contract
# 1. The primary period definition follows the China-policy transition:
#    2017-2019, 2020-2022, and 2023-2024.
# 2. The WHO PHEIC boundary is an alternative definition, not a replacement.
# 3. Spring 2023 is handled in three prespecified ways: retained, adjusted by
#    an externally anchored February-May indicator, and excluded.
# 4. All comparisons use the same -5 C versus 0 C altitude-modification
#    estimand and the updated 2-m air-temperature exposure.
# 5. The future-temperature negative control conditions on the observed lag
#    0-2 exposure and uses non-overlapping future weeks t+3 to t+5.

source(file.path("03_robustness_validation", "92_robustness.R"))

analysis_id <- "result5_temporal_supplement_20260902"
checkpoint_dir <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", analysis_id
)
output_dir <- file.path("outputs", analysis_id)
source_dir <- file.path(output_dir, "source_data")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

write_csv_utf8 <- function(x, filename) {
  write.csv(
    x, file.path(source_dir, filename), row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

fit_or_load_r5 <- function(path, formula, data, priors = NULL) {
  if (file.exists(path)) return(readRDS(path))
  fit <- if (is.null(priors)) {
    glmmTMB(
      formula, family = betabinomial(link = "logit"), data = data,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 6000, eval.max = 6000)
      )
    )
  } else {
    glmmTMB(
      formula, family = betabinomial(link = "logit"), data = data,
      priors = priors,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 6000, eval.max = 6000)
      )
    )
  }
  saveRDS(fit, path)
  fit
}

joint_wald_terms_r5 <- function(fit, terms) {
  fixed <- fixef(fit)$cond
  variance <- as.matrix(vcov(fit)$cond)
  keep <- intersect(terms, names(fixed))
  coefficient <- fixed[keep]
  covariance <- variance[keep, keep, drop = FALSE]
  statistic <- as.numeric(
    t(coefficient) %*% solve(covariance, coefficient)
  )
  tibble(
    statistic = statistic,
    df = length(keep),
    p_value = pchisq(statistic, df = length(keep), lower.tail = FALSE)
  )
}

model_qa_row <- function(label, fit) {
  tibble(
    model = label,
    observations = nobs(fit),
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit),
    convergence = fit$fit$convergence,
    pdHess = fit$sdr$pdHess
  )
}

# -------------------------------------------------------------------------
# A. Primary China-policy phase definition and data support.
# -------------------------------------------------------------------------
phase_support_summary <- phase_model_dat |>
  filter(!is.na(phase)) |>
  group_by(phase) |>
  summarise(
    start_week = min(week_start),
    end_week = max(week_start),
    surveillance_weeks = n_distinct(week_start),
    county_weeks = n(),
    tests = sum(tests),
    positives = sum(positive),
    positivity = positives / tests,
    temperature_p05 = quantile(temperature, 0.05, na.rm = TRUE),
    temperature_median = median(temperature, na.rm = TRUE),
    temperature_p95 = quantile(temperature, 0.95, na.rm = TRUE),
    event_county_weeks = sum(spring2023_event == 1),
    .groups = "drop"
  )

phase_altitude_results <- phase_results |>
  filter(estimand == "Altitude modification of -5 vs 0 C") |>
  mutate(
    phase_definition = "China-policy calendar phases",
    contrast = "-5 \u00B0C versus 0 \u00B0C",
    metric = "Ratio of odds ratios per 1,000 m"
  )

phase_effect_results <- phase_results |>
  filter(estimand == "-5 vs 0 C at 2,980 m") |>
  mutate(
    phase_definition = "China-policy calendar phases",
    contrast = "-5 \u00B0C versus 0 \u00B0C",
    metric = "Odds ratio at 2,980 m"
  )

phase_heterogeneity_clean <- phase_tests |>
  mutate(phase_definition = "China-policy calendar phases")

# -------------------------------------------------------------------------
# B. Spring-2023 event handling in the full sample and post-transition phase.
# -------------------------------------------------------------------------
extract_altitude_modification <- function(fit, estimator = c("primary", "phase")) {
  estimator <- match.arg(estimator)
  if (estimator == "primary") {
    estimate_primary_contrast_generic(
      extract_glmm_fixed(fit), -5, 0, 1, TRUE
    )
  } else {
    estimate_main_contrast(fit, -5, 0, 1, TRUE)
  }
}

full_event_fits <- list(
  "Peak retained" = primary_bundle$fits$positive_full,
  "Event indicator adjusted" = primary_bundle$fits$positive_event_adjusted,
  "February-May 2023 excluded" = primary_bundle$fits$positive_excluded
)

full_event_results <- bind_rows(lapply(names(full_event_fits), function(label) {
  fit <- full_event_fits[[label]]
  bind_cols(
    tibble(scope = "Full 2017-2024 sample", event_handling = label),
    extract_altitude_modification(fit, "primary"),
    joint_wald_terms_r5(fit, primary_int_names) |>
      rename(global_wald = statistic, global_df = df, global_p = p_value),
    tibble(observations = nobs(fit))
  )
}))

full_event_A_fits <- list(
  "Peak retained" = primary_bundle$fits$fluA_only_full,
  "Event indicator adjusted" = primary_bundle$fits$fluA_only_event_adjusted,
  "February-May 2023 excluded" = primary_bundle$fits$fluA_only_excluded
)

full_event_A_results <- bind_rows(lapply(names(full_event_A_fits), function(label) {
  fit <- full_event_A_fits[[label]]
  bind_cols(
    tibble(
      outcome = "Influenza A-only positivity",
      scope = "Full 2017-2024 sample", event_handling = label
    ),
    extract_altitude_modification(fit, "primary"),
    joint_wald_terms_r5(fit, primary_int_names) |>
      rename(global_wald = statistic, global_df = df, global_p = p_value),
    tibble(observations = nobs(fit))
  )
}))

post_label <- "Post-transition (2023-2024)"
post_dat <- phase_model_dat |>
  filter(phase == post_label) |>
  droplevels()

post_event_formula <- update(
  phase_formula, . ~ . + spring2023_event
)
post_event_fit <- fit_or_load_r5(
  file.path(checkpoint_dir, "post_transition_event_adjusted.rds"),
  post_event_formula, post_dat, phase_priors
)
post_excluded_dat <- post_dat |>
  filter(spring2023_event == 0) |>
  droplevels()
post_excluded_fit <- fit_or_load_r5(
  file.path(checkpoint_dir, "post_transition_spring2023_excluded.rds"),
  phase_formula, post_excluded_dat, phase_priors
)

post_event_fits <- list(
  "Peak retained" = phase_fits[[post_label]],
  "Event indicator adjusted" = post_event_fit,
  "February-May 2023 excluded" = post_excluded_fit
)

post_event_results <- bind_rows(lapply(names(post_event_fits), function(label) {
  fit <- post_event_fits[[label]]
  bind_cols(
    tibble(scope = "Post-transition 2023-2024", event_handling = label),
    extract_altitude_modification(fit, "phase"),
    joint_wald_terms_r5(fit, main_z_names) |>
      rename(global_wald = statistic, global_df = df, global_p = p_value),
    tibble(observations = nobs(fit))
  )
}))

event_sensitivity_results <- bind_rows(
  full_event_results, post_event_results
) |>
  mutate(
    outcome = "Overall laboratory positivity",
    contrast = "-5 \u00B0C versus 0 \u00B0C",
    metric = "Ratio of odds ratios per 1,000 m"
  ) |>
  select(
    outcome, scope, event_handling, contrast, metric,
    estimate, low, high, global_wald, global_df, global_p,
    observations, log_estimate, standard_error
  )

full_event_temperature_effects <- bind_rows(lapply(names(full_event_fits), function(label) {
  fit <- full_event_fits[[label]]
  bind_cols(
    tibble(
      outcome = "Overall laboratory positivity",
      scope = "Full 2017-2024 sample", event_handling = label,
      contrast = "-5 \u00B0C versus 0 \u00B0C",
      metric = "Odds ratio at 2,980 m"
    ),
    estimate_primary_contrast_generic(
      extract_glmm_fixed(fit), -5, 0,
      primary_anchors$altitude_1000_c[
        primary_anchors$altitude_anchor == "P50"
      ],
      FALSE
    ),
    tibble(observations = nobs(fit))
  )
}))

post_event_temperature_effects <- bind_rows(lapply(names(post_event_fits), function(label) {
  fit <- post_event_fits[[label]]
  bind_cols(
    tibble(
      outcome = "Overall laboratory positivity",
      scope = "Post-transition 2023-2024", event_handling = label,
      contrast = "-5 \u00B0C versus 0 \u00B0C",
      metric = "Odds ratio at 2,980 m"
    ),
    estimate_main_contrast(
      fit, -5, 0,
      anchors$altitude_1000_c[anchors$altitude_anchor == "P50"],
      FALSE
    ),
    tibble(observations = nobs(fit))
  )
}))

event_temperature_effect_results <- bind_rows(
  full_event_temperature_effects, post_event_temperature_effects
)

# -------------------------------------------------------------------------
# C. Alternative WHO PHEIC phase boundary.
# -------------------------------------------------------------------------
who_phase_dat <- phase_model_dat |>
  mutate(
    who_phase = case_when(
      week_start < as.Date("2020-02-03") ~ "Pre-PHEIC",
      week_start <= as.Date("2023-05-01") ~ "WHO emergency phase",
      TRUE ~ "Post-PHEIC"
    ),
    who_phase = factor(
      who_phase,
      levels = c("Pre-PHEIC", "WHO emergency phase", "Post-PHEIC")
    ),
    calendar_year_factor = factor(calendar_year)
  )

who_phase_fits <- list()
for (label in levels(who_phase_dat$who_phase)) {
  dat <- who_phase_dat |>
    filter(who_phase == label) |>
    droplevels()
  safe_label <- gsub("[^A-Za-z0-9]+", "_", label)
  who_phase_fits[[label]] <- fit_or_load_r5(
    file.path(checkpoint_dir, paste0("who_phase_", safe_label, ".rds")),
    phase_formula, dat, phase_priors
  )
}

who_phase_results <- bind_rows(lapply(names(who_phase_fits), function(label) {
  fit <- who_phase_fits[[label]]
  bind_rows(
    bind_cols(
      tibble(
        phase_definition = "WHO PHEIC boundary",
        phase = label,
        estimand = "-5 vs 0 \u00B0C at 2,980 m",
        effect_measure = "Odds ratio"
      ),
      estimate_main_contrast(
        fit, -5, 0,
        anchors$altitude_1000_c[anchors$altitude_anchor == "P50"],
        FALSE
      )
    ),
    bind_cols(
      tibble(
        phase_definition = "WHO PHEIC boundary",
        phase = label,
        estimand = "Altitude modification of -5 vs 0 \u00B0C",
        effect_measure = "Ratio of odds ratios per 1,000 m"
      ),
      estimate_main_contrast(fit, -5, 0, 1, TRUE)
    )
  ) |>
    mutate(
      observations = nobs(fit),
      global_interaction_p = joint_wald_terms_r5(
        fit, main_z_names
      )$p_value
    )
}))

heterogeneity_q <- function(dat) {
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

who_phase_heterogeneity <- bind_rows(lapply(
  split(who_phase_results, who_phase_results$estimand),
  heterogeneity_q
)) |>
  mutate(phase_definition = "WHO PHEIC boundary")

who_phase_support <- who_phase_dat |>
  group_by(who_phase) |>
  summarise(
    start_week = min(week_start),
    end_week = max(week_start),
    surveillance_weeks = n_distinct(week_start),
    county_weeks = n(),
    tests = sum(tests),
    positives = sum(positive),
    positivity = positives / tests,
    temperature_p05 = quantile(temperature, .05, na.rm = TRUE),
    temperature_median = median(temperature, na.rm = TRUE),
    temperature_p95 = quantile(temperature, .95, na.rm = TRUE),
    event_county_weeks = sum(spring2023_event == 1),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# D. Conditional future-temperature negative control.
# -------------------------------------------------------------------------
weekly_grid_future_conditional <- weekly_grid |>
  arrange(county, week_index) |>
  group_by(county) |>
  mutate(future_anomaly_t3_t5 = lead(anomaly_for_basis, 5)) |>
  ungroup()

cb_future_conditional <- crossbasis(
  weekly_grid_future_conditional$future_anomaly_t3_t5,
  lag = c(0, 2),
  argvar = list(fun = "ns", df = 3),
  arglag = list(fun = "poly", degree = 1),
  group = weekly_grid_future_conditional$county
)
future_cond_names <- paste0(
  "r5_future_cond_ns_", seq_len(ncol(cb_future_conditional))
)
future_cond_int_names <- paste0(
  "r5_future_cond_alt_", seq_len(ncol(cb_future_conditional))
)
future_cond_df <- as.data.frame(
  cb_future_conditional[confirm_dat$exposure_row_id, , drop = FALSE]
)
names(future_cond_df) <- future_cond_names
future_cond_dat <- bind_cols(confirm_dat, future_cond_df) |>
  filter(complete.cases(across(all_of(future_cond_names))))
for (i in seq_along(future_cond_names)) {
  future_cond_dat[[future_cond_int_names[[i]]]] <-
    future_cond_dat[[future_cond_names[[i]]]] *
    future_cond_dat$altitude_1000_c
}

conditional_base_formula <- formula_from_terms(c(
  spline_names, spline_int_names, common_rhs
))
conditional_future_main_formula <- formula_from_terms(c(
  spline_names, spline_int_names,
  future_cond_names, common_rhs
))
conditional_future_full_formula <- formula_from_terms(c(
  spline_names, spline_int_names,
  future_cond_names, future_cond_int_names, common_rhs
))

fit_conditional_base <- fit_or_load_r5(
  file.path(checkpoint_dir, "conditional_future_base.rds"),
  conditional_base_formula, future_cond_dat
)
fit_conditional_future_main <- fit_or_load_r5(
  file.path(checkpoint_dir, "conditional_future_main.rds"),
  conditional_future_main_formula, future_cond_dat
)
fit_conditional_future_full <- fit_or_load_r5(
  file.path(checkpoint_dir, "conditional_future_interaction.rds"),
  conditional_future_full_formula, future_cond_dat
)

conditional_future_tests <- bind_rows(
  safe_lrt(
    fit_conditional_base, fit_conditional_future_main,
    "Incremental future-temperature main effect, conditional on lag 0-2"
  ),
  safe_lrt(
    fit_conditional_future_main, fit_conditional_future_full,
    "Incremental future-temperature-by-altitude interaction, conditional on lag 0-2"
  )
) |>
  mutate(
    future_window = "Outcome week t evaluated against future weeks t+3 to t+5",
    observations = nobs(fit_conditional_future_full)
  )

future_correlation <- weekly_grid_future_conditional |>
  filter(
    is.finite(anomaly_for_basis), is.finite(future_anomaly_t3_t5)
  ) |>
  summarise(
    comparison = "Current anomaly versus five-week lead",
    correlation = cor(anomaly_for_basis, future_anomaly_t3_t5),
    observations = n()
  )

# -------------------------------------------------------------------------
# E. Model QA and tabular source files.
# -------------------------------------------------------------------------
new_model_qa <- bind_rows(
  model_qa_row("Post-transition event adjusted", post_event_fit),
  model_qa_row("Post-transition spring 2023 excluded", post_excluded_fit),
  bind_rows(lapply(names(who_phase_fits), function(label) {
    model_qa_row(paste("WHO phase", label), who_phase_fits[[label]])
  })),
  model_qa_row("Conditional future base", fit_conditional_base),
  model_qa_row("Conditional future main", fit_conditional_future_main),
  model_qa_row("Conditional future interaction", fit_conditional_future_full)
)

write_csv_utf8(phase_support_summary, "table_s4_china_phase_support.csv")
write_csv_utf8(phase_altitude_results, "table_s4_china_phase_altitude_results.csv")
write_csv_utf8(phase_effect_results, "table_s4_china_phase_effect_results.csv")
write_csv_utf8(phase_heterogeneity_clean, "table_s4_china_phase_heterogeneity.csv")
write_csv_utf8(event_sensitivity_results, "table_s4_spring2023_event_sensitivity.csv")
write_csv_utf8(event_temperature_effect_results, "table_s4_spring2023_event_temperature_effects.csv")
write_csv_utf8(full_event_A_results, "table_s4_influenzaA_event_sensitivity.csv")
write_csv_utf8(who_phase_support, "table_s4_who_phase_support.csv")
write_csv_utf8(who_phase_results, "table_s4_who_phase_results.csv")
write_csv_utf8(who_phase_heterogeneity, "table_s4_who_phase_heterogeneity.csv")
write_csv_utf8(conditional_future_tests, "table_s5_conditional_future_negative_control.csv")
write_csv_utf8(future_correlation, "table_s5_future_temperature_correlation.csv")
write_csv_utf8(new_model_qa, "result5_supplement_model_qa.csv")

# -------------------------------------------------------------------------
# F. Supplementary Figure S4.
# -------------------------------------------------------------------------
palette <- c(
  blue = "#2F6FAE",
  blue_light = "#92B6D5",
  orange = "#C58A35",
  red = "#C85A4A",
  grey_dark = "#4D5660",
  grey_mid = "#929CA5",
  grey_light = "#DDE3E8"
)

theme_r5 <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#30343B"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#30343B"),
      axis.title = element_text(size = 7.2, colour = "#20252B"),
      axis.text = element_text(size = 6.6, colour = "#30343B"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6.2),
      plot.title = element_text(size = 8.2, face = "bold", colour = "#18212B"),
      plot.subtitle = element_text(size = 6.5, colour = "#657482"),
      plot.tag = element_text(size = 9, face = "bold", colour = "#111820"),
      plot.margin = margin(4, 6, 4, 4),
      panel.grid = element_blank()
    )
}

timeline <- tibble::tribble(
  ~track, ~label, ~start, ~end, ~segment,
  "Primary phases", "Pre-pandemic", as.Date("2017-01-01"), as.Date("2019-12-31"), "Pre-pandemic",
  "Primary phases", "Pandemic disruption", as.Date("2020-01-01"), as.Date("2022-12-31"), "Pandemic disruption",
  "Primary phases", "Post-transition", as.Date("2023-01-01"), as.Date("2024-12-31"), "Post-transition",
  "Event window", "Spring 2023", as.Date("2023-02-01"), as.Date("2023-05-31"), "Spring 2023 event",
  "WHO boundary", "PHEIC", as.Date("2020-01-30"), as.Date("2023-05-05"), "WHO PHEIC"
) |>
  mutate(
    track = factor(
      track,
      levels = c("WHO boundary", "Event window", "Primary phases")
    ),
    midpoint = start + (end - start) / 2
  )

timeline_colors <- c(
  "Pre-pandemic" = palette[["grey_mid"]],
  "Pandemic disruption" = palette[["orange"]],
  "Post-transition" = palette[["blue"]],
  "Spring 2023 event" = palette[["red"]],
  "WHO PHEIC" = palette[["grey_dark"]]
)

p_a <- ggplot(timeline, aes(y = track)) +
  geom_segment(
    aes(x = start, xend = end, yend = track, colour = segment),
    linewidth = 7.2, lineend = "butt"
  ) +
  geom_text(
    data = timeline |> filter(track == "Primary phases"),
    aes(x = midpoint, label = label),
    colour = "white", size = 2.05, fontface = "bold"
  ) +
  geom_text(
    data = timeline |> filter(track != "Primary phases"),
    aes(x = midpoint, label = label),
    colour = "#25303A", size = 2.05, vjust = -1.25
  ) +
  scale_colour_manual(values = timeline_colors, guide = "none") +
  scale_x_date(
    date_breaks = "1 year", date_labels = "%Y",
    limits = as.Date(c("2016-11-01", "2025-02-01")),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Externally anchored temporal definitions",
    subtitle = "The spring-2023 epidemic is analysed as an event, not reassigned post hoc",
    x = NULL, y = NULL
  ) +
  theme_r5() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 6.4),
    plot.margin = margin(4, 8, 7, 4)
  )

phase_plot_data <- phase_altitude_results |>
  mutate(
    display = recode(
      phase,
      "Pre-pandemic (2017-2019)" = "Pre-pandemic\n2017-2019",
      "Pandemic disruption (2020-2022)" = "Pandemic disruption\n2020-2022",
      "Post-transition (2023-2024)" = "Post-transition\n2023-2024"
    ),
    display = factor(
      display,
      levels = c(
        "Post-transition\n2023-2024",
        "Pandemic disruption\n2020-2022",
        "Pre-pandemic\n2017-2019"
      )
    ),
    highlight = if_else(
      phase == "Post-transition (2023-2024)", "Post-transition", "Other phases"
    ),
    estimate_label = sprintf("%.2f (%.2f-%.2f)", estimate, low, high)
  )

phase_het_p <- phase_heterogeneity_clean |>
  filter(estimand == "Altitude modification of -5 vs 0 C") |>
  pull(p_value)

p_b <- ggplot(phase_plot_data, aes(x = estimate, y = display)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, colour = palette[["grey_mid"]]) +
  geom_errorbarh(
    aes(xmin = low, xmax = high, colour = highlight),
    height = 0, linewidth = 0.75
  ) +
  geom_point(
    aes(colour = highlight), shape = 21, fill = "white",
    stroke = 0.8, size = 2.6
  ) +
  geom_text(
    aes(x = estimate, label = estimate_label),
    hjust = 0.5, vjust = -1.15, size = 2.15, colour = "#30343B"
  ) +
  scale_colour_manual(values = c(
    "Post-transition" = palette[["blue"]],
    "Other phases" = palette[["grey_dark"]]
  ), guide = "none") +
  scale_x_log10(
    limits = c(0.55, 2.70), breaks = c(0.6, 0.8, 1, 1.5, 2.0, 2.5),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    title = "Modification across epidemic periods",
    subtitle = sprintf("Cochran Q heterogeneity P = %.3f", phase_het_p),
    x = "Ratio of odds ratios per 1,000 m", y = NULL
  ) +
  theme_r5() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

event_plot_data <- event_sensitivity_results |>
  mutate(
    event_handling = factor(
      event_handling,
      levels = c(
        "February-May 2023 excluded",
        "Event indicator adjusted",
        "Peak retained"
      )
    ),
    scope_short = recode(
      scope,
      "Full 2017-2024 sample" = "Full sample",
      "Post-transition 2023-2024" = "Post-transition"
    ),
    scope_short = factor(
      scope_short, levels = c("Post-transition", "Full sample")
    ),
    estimate_label = sprintf("%.2f (%.2f-%.2f)", estimate, low, high)
  )

p_c <- ggplot(
  event_plot_data,
  aes(x = estimate, y = event_handling, colour = scope_short)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, colour = palette[["grey_mid"]]) +
  geom_errorbarh(
    aes(xmin = low, xmax = high),
    position = position_dodge(width = 0.34), height = 0, linewidth = 0.7
  ) +
  geom_point(
    position = position_dodge(width = 0.34),
    shape = 21, fill = "white", stroke = 0.8, size = 2.35
  ) +
  geom_text(
    aes(x = estimate, label = estimate_label),
    hjust = 0.5, vjust = -1.15, size = 2.05, colour = "#30343B",
    show.legend = FALSE
  ) +
  facet_grid(scope_short ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(
    "Full sample" = palette[["grey_dark"]],
    "Post-transition" = palette[["blue"]]
  ), guide = "none") +
  scale_x_log10(
    limits = c(0.75, 2.40), breaks = c(0.8, 1, 1.25, 1.5, 2),
    labels = label_number(accuracy = 0.05),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Sensitivity to the spring-2023 epidemic",
    subtitle = "Same prespecified contrast under three event-handling strategies",
    x = "Ratio of odds ratios per 1,000 m", y = NULL
  ) +
  theme_r5() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 6.4, face = "bold", angle = 0),
    panel.spacing.y = unit(4, "pt")
  )

figure_s4 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(0.68, 1.42)) +
  plot_annotation(
    title = "Temporal robustness of altitude modification",
    subtitle = paste0(
      "Overall laboratory positivity; 2-m air temperature; cumulative lag 0-2; ",
      "contrast: -5 \u00B0C versus 0 \u00B0C"
    ),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(
        family = "Arial", size = 10.2, face = "bold", colour = "#14202B"
      ),
      plot.subtitle = element_text(
        family = "Arial", size = 7, colour = "#657482"
      ),
      plot.tag = element_text(
        family = "Arial", size = 9, face = "bold", colour = "#111820"
      ),
      plot.margin = margin(5, 6, 4, 5)
    )
  )

figure_base <- file.path(
  output_dir, "Supplementary_Figure_S4"
)
width_in <- 183 / 25.4
height_in <- 132 / 25.4

svglite::svglite(paste0(figure_base, ".svg"), width = width_in, height = height_in)
print(figure_s4)
dev.off()

grDevices::cairo_pdf(
  paste0(figure_base, ".pdf"), width = width_in, height = height_in,
  family = "Arial"
)
print(figure_s4)
dev.off()

temporary_tiff <- file.path(tempdir(), "Supplementary_Figure_S4.tiff")
ragg::agg_tiff(
  temporary_tiff, width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure_s4)
dev.off()
if (!file.copy(temporary_tiff, paste0(figure_base, ".tiff"), overwrite = TRUE)) {
  stop("Unable to copy the TIFF export to the final output folder.")
}

temporary_png <- file.path(tempdir(), "Supplementary_Figure_S4.png")
ragg::agg_png(
  temporary_png, width = width_in, height = height_in,
  units = "in", res = 300
)
print(figure_s4)
dev.off()
if (!file.copy(temporary_png, paste0(figure_base, ".png"), overwrite = TRUE)) {
  stop("Unable to copy the PNG export to the final output folder.")
}

legend_text <- paste(
  "Supplementary Figure S4. Temporal robustness of altitude modification.",
  "(a) Externally anchored temporal definitions used for sensitivity analyses.",
  "The primary three-period classification follows calendar years surrounding",
  "the China-policy transition, whereas the WHO public-health-emergency boundary",
  "is evaluated as an alternative definition. The February-May 2023 influenza A",
  "epidemic is treated as a prespecified event window. (b) Period-specific ratios",
  "of odds ratios per 1,000-m increase in altitude for the cumulative lag 0-2",
  "contrast of -5 \u00B0C versus 0 \u00B0C. Points are estimates and horizontal bars are",
  "95% confidence intervals; the dashed line denotes no altitude modification.",
  "The displayed P value is the inverse-variance Cochran Q test for heterogeneity",
  "across periods. (c) The same altitude-modification estimand after retaining,",
  "adjusting for, or excluding the externally defined spring-2023 event window,",
  "shown for the full sample and the post-transition period. Estimates are from",
  "hierarchical beta-binomial distributed-lag models using updated weekly mean",
  "2-m air temperature. Period-specific findings are exploratory because the",
  "post-transition stratum contains two calendar years.",
  sep = " "
)
writeLines(
  strwrap(legend_text, width = 110),
  file.path(output_dir, "Supplementary_Figure_S4_legend.txt"),
  useBytes = TRUE
)

figure_qa <- tibble(
  check = c(
    "All newly fitted models converged",
    "All newly fitted models have positive-definite Hessians",
    "All plotted confidence intervals are finite and positive",
    "Phase contrast has observational support",
    "Figure source rows reconcile with plotted rows",
    "Conditional future negative control uses non-overlapping weeks"
  ),
  pass = c(
    all(new_model_qa$convergence == 0),
    all(new_model_qa$pdHess),
    all(is.finite(event_plot_data$low)) && all(event_plot_data$low > 0) &&
      all(is.finite(phase_plot_data$low)) && all(phase_plot_data$low > 0),
    all(phase_support$contrast_supported),
    nrow(event_plot_data) == 6 && nrow(phase_plot_data) == 3,
    TRUE
  )
)
write_csv_utf8(figure_qa, "figure_s4_qa.csv")

if (!all(figure_qa$pass)) {
  warning("At least one Figure S4 QA check failed; outputs retained for audit.")
}
if (any(new_model_qa$convergence != 0) || any(!new_model_qa$pdHess)) {
  warning("At least one new Result 5 model requires convergence review.")
}

message("Completed Result 5 temporal supplement, conditional negative control, and Figure S4.")
