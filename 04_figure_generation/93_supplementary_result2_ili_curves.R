suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(glmmTMB)
  library(dlnm)
  library(svglite)
  library(ragg)
  library(scales)
})

set.seed(20260902)

# Figure contract ----------------------------------------------------------
# Core conclusion: within the identical 18-county and county-week footprint,
# laboratory positivity and ILI consultation rates show outcome-specific
# nonlinear temperature-response curves, and their altitude pattern should be
# judged only within the observed ILI altitude support.
# Archetype: quantitative grid (outcomes by supported altitude anchors).
# Backend: R only.
# Outputs: editable SVG/PDF, 600-dpi TIFF, PNG preview, source CSV, model RDS,
# legend, figure contract, and QA CSV.

analysis_id <- "result2_supplementary_ili_aligned_20260902"
bundle_path <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825", "v4_t2m_model_bundle.rds"
)
model_dir <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", analysis_id
)
figure_dir <- file.path(
  "result0715", "supplementary_figures", "v4_t2m_20260825", analysis_id
)
source_dir <- file.path(
  "result0715", "source_data", paste0("v4_", analysis_id)
)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
primary_dat <- bundle$model_data
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
z_int_names <- bundle$z_int_names
interaction_scale <- bundle$interaction_scale
altitude_mean <- bundle$altitude_mean
reference_temperature <- 0
degree_c <- paste0(intToUtf8(176), "C")

# Exact endpoint alignment -------------------------------------------------
# primary_dat is already the complete laboratory model sample. Requiring a
# valid ILI numerator and denominator therefore creates the exact intersection
# of the laboratory and ILI county-weeks.
matched_dat <- primary_dat |>
  filter(
    ili_available %in% TRUE,
    is.finite(ILI), is.finite(visits),
    visits >= 1, ILI >= 0, ILI <= visits
  ) |>
  droplevels()

county_altitudes <- matched_dat |>
  distinct(county, altitude) |>
  arrange(altitude)

if (nrow(county_altitudes) != 18L) {
  stop("The exact matched footprint does not contain the expected 18 counties.")
}

anchor_values <- as.numeric(
  quantile(county_altitudes$altitude, c(.10, .50, .90), names = FALSE)
)
anchor_display <- floor(anchor_values + .5)
altitude_anchors <- tibble(
  altitude_anchor = factor(c("P10", "P50", "P90"),
                           levels = c("P10", "P50", "P90")),
  altitude = anchor_values,
  altitude_1000_c = (anchor_values - altitude_mean) / 1000,
  altitude_label = factor(
    sprintf("%s: %s m", c("P10", "P50", "P90"),
            format(anchor_display, big.mark = ",", scientific = FALSE)),
    levels = sprintf("%s: %s m", c("P10", "P50", "P90"),
                     format(anchor_display, big.mark = ",", scientific = FALSE))
  )
)

# Local temperature support follows the Result 2 convention: county-weeks
# within +/-250 m of each continuous-altitude display anchor.
local_support <- bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
  a <- altitude_anchors$altitude[[i]]
  z <- matched_dat |> filter(abs(altitude - a) <= 250)
  if (nrow(z) < 100) stop("Insufficient local support around an ILI altitude anchor.")
  tibble(
    altitude_anchor = altitude_anchors$altitude_anchor[[i]],
    altitude = a,
    local_county_weeks = nrow(z),
    local_counties = n_distinct(z$county),
    p01 = as.numeric(quantile(z$temperature, .01)),
    p05 = as.numeric(quantile(z$temperature, .05)),
    p10 = as.numeric(quantile(z$temperature, .10)),
    median = median(z$temperature),
    p90 = as.numeric(quantile(z$temperature, .90)),
    p95 = as.numeric(quantile(z$temperature, .95)),
    p99 = as.numeric(quantile(z$temperature, .99))
  )
}))

if (any(reference_temperature < local_support$p10 |
        reference_temperature > local_support$p90)) {
  stop("The 0 C reference is outside local P10-P90 support at an ILI anchor.")
}

# Comparable hierarchical beta-binomial models ---------------------------
season_names <- grep("^season_(sin|cos)_", names(matched_dat), value = TRUE)
trend_names <- grep("^portable_trend_", names(matched_dat), value = TRUE)
rhs_full <- c(
  cb_names, "altitude_1000_c", z_int_names,
  season_names, trend_names, "pandemic_period",
  "(1 | city)", "(1 | county_nested)"
)
formula_lab <- as.formula(paste(
  "cbind(positive, tests - positive) ~", paste(rhs_full, collapse = " + ")
))
formula_ili <- as.formula(paste(
  "cbind(ILI, visits - ILI) ~", paste(rhs_full, collapse = " + ")
))
prior_spec <- bundle$fits$positive_full$modelInfo$priors
control <- glmmTMBControl(optCtrl = list(iter.max = 5000, eval.max = 5000))

fit_or_load <- function(path, formula) {
  if (file.exists(path)) return(readRDS(path))
  fit <- glmmTMB(
    formula, family = betabinomial(link = "logit"), data = matched_dat,
    priors = prior_spec, control = control
  )
  saveRDS(fit, path)
  fit
}

fit_lab <- fit_or_load(
  file.path(model_dir, "matched_18county_laboratory_hierarchical.rds"),
  formula_lab
)
fit_ili <- fit_or_load(
  file.path(model_dir, "matched_18county_ili_hierarchical.rds"),
  formula_ili
)

extract_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, variance = as.matrix(vcov(fit)$cond))
}

joint_wald <- function(fit) {
  fixed <- extract_fixed(fit)
  beta_int <- fixed$beta[z_int_names]
  variance_int <- fixed$variance[z_int_names, z_int_names, drop = FALSE]
  statistic <- as.numeric(t(beta_int) %*% solve(variance_int, beta_int))
  tibble(
    chisq = statistic, df = length(beta_int),
    p_value = pchisq(statistic, df = length(beta_int), lower.tail = FALSE)
  )
}

wald_lab <- joint_wald(fit_lab) |> mutate(outcome = "Laboratory positivity")
wald_ili <- joint_wald(fit_ili) |> mutate(outcome = "ILI consultation rate")

basis_matrix_cumulative <- function(temperatures, centre = 0) {
  temperatures <- sort(unique(c(temperatures, centre)))
  original_names <- colnames(cb_all)
  output <- matrix(
    NA_real_, nrow = length(temperatures), ncol = length(original_names),
    dimnames = list(as.character(temperatures), cb_names)
  )
  zero_variance <- matrix(
    0, nrow = length(original_names), ncol = length(original_names),
    dimnames = list(original_names, original_names)
  )
  for (j in seq_along(original_names)) {
    coefficient <- rep(0, length(original_names))
    names(coefficient) <- original_names
    coefficient[[j]] <- 1
    prediction <- crosspred(
      cb_all, coef = coefficient, vcov = zero_variance,
      model.link = "logit", at = temperatures, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

estimate_curve <- function(fit, basis, altitude_1000_c) {
  fixed <- extract_fixed(fit)
  design <- matrix(
    0, nrow = nrow(basis), ncol = length(fixed$beta),
    dimnames = list(NULL, names(fixed$beta))
  )
  design[, cb_names] <- basis
  design[, z_int_names] <-
    sweep(basis, 2, interaction_scale, "/") * altitude_1000_c
  log_estimate <- as.numeric(design %*% fixed$beta)
  standard_error <- sqrt(pmax(rowSums((design %*% fixed$variance) * design), 0))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

curve_rows <- vector("list", nrow(altitude_anchors) * 2)
row_index <- 1L
for (i in seq_len(nrow(altitude_anchors))) {
  anchor <- altitude_anchors[i, ]
  support <- local_support |>
    filter(altitude_anchor == anchor$altitude_anchor)
  temperatures <- sort(unique(c(
    seq(
      floor(support$p05[[1]] * 4) / 4,
      ceiling(support$p95[[1]] * 4) / 4,
      by = .25
    ),
    reference_temperature
  )))
  basis <- basis_matrix_cumulative(temperatures, reference_temperature)
  for (spec in list(
    list(outcome = "Laboratory positivity", fit = fit_lab),
    list(outcome = "ILI consultation rate", fit = fit_ili)
  )) {
    curve_rows[[row_index]] <- bind_cols(
      tibble(
        outcome = spec$outcome,
        altitude_anchor = anchor$altitude_anchor,
        altitude = anchor$altitude,
        altitude_label = anchor$altitude_label,
        reference_temperature = reference_temperature,
        support_p05 = support$p05[[1]],
        support_p95 = support$p95[[1]],
        local_county_weeks = support$local_county_weeks[[1]],
        local_counties = support$local_counties[[1]],
        temperature = as.numeric(rownames(basis))
      ),
      estimate_curve(spec$fit, basis, anchor$altitude_1000_c[[1]])
    )
    row_index <- row_index + 1L
  }
}
curves <- bind_rows(curve_rows) |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("Laboratory positivity", "ILI consultation rate")
    ),
    altitude_label = factor(
      altitude_label, levels = levels(altitude_anchors$altitude_label)
    )
  )

key_contrasts <- curves |>
  filter(temperature %in% c(-5, 5)) |>
  mutate(
    comparison = if_else(
      temperature == -5, "-5 C versus 0 C", "5 C versus 0 C"
    )
  ) |>
  select(
    outcome, altitude_anchor, altitude, comparison,
    estimate, low, high, standard_error
  )

# Source data and model summary -------------------------------------------
model_summary <- bind_rows(
  tibble(
    outcome = c("Laboratory positivity", "ILI consultation rate"),
    observations = c(nobs(fit_lab), nobs(fit_ili)),
    counties = n_distinct(matched_dat$county),
    events = c(sum(matched_dat$positive), sum(matched_dat$ILI)),
    denominators = c(sum(matched_dat$tests), sum(matched_dat$visits)),
    convergence = c(fit_lab$fit$convergence, fit_ili$fit$convergence),
    pdHess = c(fit_lab$sdr$pdHess, fit_ili$sdr$pdHess),
    AIC = c(AIC(fit_lab), AIC(fit_ili)),
    BIC = c(BIC(fit_lab), BIC(fit_ili))
  ) |>
    left_join(bind_rows(wald_lab, wald_ili), by = "outcome")
)

write.csv(
  curves,
  file.path(source_dir, "supplementary_figure_s1_matched_temperature_curves.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  local_support,
  file.path(source_dir, "supplementary_figure_s1_local_support.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  model_summary,
  file.path(source_dir, "supplementary_figure_s1_model_summary.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  key_contrasts,
  file.path(source_dir, "supplementary_figure_s1_key_contrasts.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# Plot --------------------------------------------------------------------
palette <- c(
  "Laboratory positivity" = "#2E6FAE",
  "ILI consultation rate" = "#C6842C"
)
fill_palette <- c(
  "Laboratory positivity" = "#C7D9E9",
  "ILI consultation rate" = "#ECD6B8"
)

plot_dat <- curves |>
  mutate(
    low_plot = low,
    high_plot = high
  )

support_labels <- plot_dat |>
  distinct(outcome, altitude_label, support_p05, support_p95) |>
  mutate(
    label = sprintf(
      "Local P5-P95: %.1f to %.1f %s",
      support_p05, support_p95, degree_c
    )
  )

figure_plot <- ggplot(
  plot_dat,
  aes(x = temperature, y = estimate, group = outcome)
) +
  geom_hline(yintercept = 1, linewidth = .35, linetype = "dashed",
             colour = "#7A8793") +
  geom_vline(xintercept = 0, linewidth = .35, linetype = "dotted",
             colour = "#444444") +
  geom_ribbon(
    aes(ymin = low_plot, ymax = high_plot, fill = outcome),
    alpha = .72, linewidth = 0, show.legend = FALSE
  ) +
  geom_line(aes(colour = outcome), linewidth = .8, show.legend = FALSE) +
  geom_text(
    data = support_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.04, vjust = 1.55,
    size = 2.15, colour = "#6B7B8C"
  ) +
  facet_grid(
    rows = vars(outcome), cols = vars(altitude_label),
    scales = "free", switch = "y"
  ) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = fill_palette) +
  scale_y_log10(
    breaks = c(.20, .25, .5, .75, 1, 1.5, 2, 4, 8, 16),
    labels = label_number(accuracy = .01)
  ) +
  labs(
    title = "Temperature-response curves in the matched ILI surveillance footprint",
    subtitle = paste0(
      "Identical 18-county and county-week sample; 2-m air temperature; ",
      "cumulative lag 0-2 weeks; reference: 0 ", degree_c,
      "; row-specific y-axis ranges"
    ),
    x = paste0("Weekly mean 2-m air temperature (", degree_c, ")"),
    y = "Cumulative odds ratio"
  ) +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 10.5, face = "bold", colour = "#172A3A", margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 6.8, colour = "#657789", margin = margin(b = 8)
    ),
    axis.title = element_text(size = 7.2, colour = "#202020"),
    axis.text = element_text(size = 6.4, colour = "#202020"),
    axis.line = element_line(linewidth = .4, colour = "#303030"),
    axis.ticks = element_line(linewidth = .35, colour = "#303030"),
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 7.2, face = "bold", colour = "#1F2D38", margin = margin(b = 3)
    ),
    strip.text.y.left = element_text(
      size = 7.4, face = "bold", angle = 90, colour = "#1F2D38",
      margin = margin(r = 5)
    ),
    panel.spacing.x = unit(10, "pt"),
    panel.spacing.y = unit(11, "pt"),
    plot.margin = margin(8, 10, 7, 8)
  )

base_name <- file.path(
  figure_dir,
  "Supplementary_Figure_S1_Result2_ILI_aligned_temperature_response_20260902"
)
width_mm <- 183
height_mm <- 122
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(base_name, ".svg"), width = width_in, height = height_in)
print(figure_plot)
dev.off()

grDevices::cairo_pdf(
  paste0(base_name, ".pdf"), width = width_in, height = height_in,
  family = "Arial"
)
print(figure_plot)
dev.off()

grDevices::tiff(
  paste0(base_name, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, bg = "white", compression = "lzw",
  type = "cairo"
)
print(figure_plot)
dev.off()

grDevices::png(
  paste0(base_name, ".png"), width = width_in, height = height_in,
  units = "in", res = 300, bg = "white", type = "cairo"
)
print(figure_plot)
dev.off()

# Legend, contract and QA --------------------------------------------------
wald_lab_text <- sprintf("chi-square(%d)=%.2f, P=%.3f",
                         wald_lab$df, wald_lab$chisq, wald_lab$p_value)
wald_ili_text <- sprintf("chi-square(%d)=%.2f, P=%.3f",
                         wald_ili$df, wald_ili$chisq, wald_ili$p_value)

legend_lines <- c(
  "Supplementary Figure S1. Temperature-response curves for laboratory influenza positivity and influenza-like illness consultation rates within the matched ILI surveillance footprint.",
  paste0("The analysis was restricted to the identical 18 counties and 5,062 county-weeks with valid laboratory and ILI denominators. Columns show continuous-altitude model slices at the P10 (2,080 m), P50 (2,563 m), and P90 (3,423 m) altitudes of the ILI-covered counties; these are display anchors rather than altitude strata. Rows show laboratory influenza positivity and ILI consultation rate, respectively. Curves are cumulative odds ratios for weekly mean 2-m air temperature over lags 0-2 weeks, relative to 0 ", degree_c, ". Shaded bands denote 95% Wald confidence intervals. Curves are displayed within the local P5-P95 temperature support among county-weeks within +/-250 m of each altitude anchor. Row-specific y-axis ranges are used because the two surveillance outcomes are not numerically interchangeable. Both outcomes were fitted using hierarchical beta-binomial models with identical temperature-lag basis functions, continuous-altitude interactions, seasonal harmonics, long-term trend, pandemic-period adjustment, and random intercepts for prefecture and county. The spring 2023 influenza peak was retained and no event indicator was included. The joint 12-df Wald tests for temperature-by-altitude interaction were:"),
  paste0("laboratory positivity, ", wald_lab_text, "; ILI consultation rate, ", wald_ili_text, "."),
  "Laboratory positivity and ILI consultation rate are complementary surveillance outcomes and their odds ratios should not be interpreted as numerically interchangeable."
)
writeLines(
  legend_lines,
  file.path(figure_dir, "Supplementary_Figure_S1_legend.md"),
  useBytes = TRUE
)

contract_lines <- c(
  "# Supplementary Figure S1 figure contract",
  "",
  "- Core conclusion: within the same 18-county and county-week footprint, laboratory positivity and ILI consultation rates have outcome-specific nonlinear temperature-response curves; altitude comparisons are restricted to the observed ILI support.",
  "- Figure archetype: quantitative 2 x 3 grid.",
  "- Backend: R only.",
  "- Final size: 183 x 122 mm.",
  "- Evidence hierarchy: matched laboratory row establishes the virological comparator; matched ILI row provides complementary syndromic evidence.",
  "- Statistical interval: 95% Wald CI.",
  "- Reviewer risk controlled: exact county-week matching; ILI-specific P10/P50/P90 anchors; no unsupported 4,000-m slice; local P5-P95 temperature support displayed; interaction tests and sample denominators reported."
)
writeLines(
  contract_lines,
  file.path(figure_dir, "Supplementary_Figure_S1_contract.md"),
  useBytes = TRUE
)

reference_check <- curves |>
  filter(abs(temperature - reference_temperature) < 1e-12) |>
  summarise(max_reference_error = max(abs(estimate - 1))) |>
  pull(max_reference_error)

qa <- tibble(
  check = c(
    "Exact matched sample is shared by both outcomes",
    "Expected 18 counties retained",
    "All fitted models converged",
    "Positive-definite Hessians",
    "Reference curves equal OR 1 at 0 C",
    "All curve estimates and intervals are finite",
    "Confidence intervals are ordered",
    "ILI display anchors lie inside observed altitude range",
    "Unsupported 4,000-m ILI anchor excluded",
    "0 C is inside local P10-P90 temperature support"
  ),
  pass = c(
    nobs(fit_lab) == nobs(fit_ili) && nobs(fit_lab) == nrow(matched_dat),
    n_distinct(matched_dat$county) == 18,
    fit_lab$fit$convergence == 0 && fit_ili$fit$convergence == 0,
    isTRUE(fit_lab$sdr$pdHess) && isTRUE(fit_ili$sdr$pdHess),
    is.finite(reference_check) && reference_check < 1e-8,
    all(is.finite(curves$estimate) & is.finite(curves$low) & is.finite(curves$high)),
    all(curves$low <= curves$estimate & curves$estimate <= curves$high),
    all(anchor_values >= min(county_altitudes$altitude) &
        anchor_values <= max(county_altitudes$altitude)),
    max(anchor_values) < 4000,
    all(reference_temperature >= local_support$p10 &
        reference_temperature <= local_support$p90)
  ),
  detail = c(
    sprintf("n=%s county-weeks for each outcome", format(nrow(matched_dat), big.mark = ",")),
    sprintf("counties=%d", n_distinct(matched_dat$county)),
    sprintf("convergence: laboratory=%d; ILI=%d", fit_lab$fit$convergence, fit_ili$fit$convergence),
    sprintf("pdHess: laboratory=%s; ILI=%s", fit_lab$sdr$pdHess, fit_ili$sdr$pdHess),
    sprintf("maximum absolute reference error=%.3g", reference_check),
    sprintf("curve rows=%d", nrow(curves)),
    "low <= estimate <= high",
    sprintf("observed altitude range: %.0f-%.0f m", min(county_altitudes$altitude), max(county_altitudes$altitude)),
    sprintf("highest display anchor=%.0f m", max(anchor_values)),
    sprintf("local common P10-P90: %.1f to %.1f C", max(local_support$p10), min(local_support$p90))
  )
)
write.csv(
  qa, file.path(figure_dir, "Supplementary_Figure_S1_QA.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

if (any(!qa$pass)) {
  stop("At least one Supplementary Figure S1 QA check failed.")
}

message("Supplementary Figure S1 completed: ", base_name)
