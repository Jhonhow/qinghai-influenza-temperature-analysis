suppressPackageStartupMessages({
  library(dlnm)
  library(glmmTMB)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Base CSV I/O avoids a Windows/readr normalization issue when the physical
# project path contains non-ASCII characters and is accessed through a drive
# alias. The files remain UTF-8 encoded and otherwise unchanged.
read_csv <- function(file, ..., show_col_types = FALSE) {
  utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(x, file, ...) {
  utils::write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
}

project_root <- Sys.getenv("R4_PROJECT_ROOT", unset = ".")

bundle_path <- file.path(project_root,
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "confirmatory_anomaly_shape_20260831", "confirmatory_anomaly_bundle.rds"
)
checkpoint_dir <- file.path(project_root,
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "confirmatory_anomaly_shape_20260831"
)
source_dir <- file.path(project_root,
  "result0715", "source_data", "v4_confirmatory_anomaly_shape_20260831"
)
gradient_path <- file.path(project_root,
  "result0715", "source_data", "v4_result4_climate_context_phase_20260826",
  "temperature_altitude_gradient_summary.csv"
)
out_dir <- Sys.getenv(
  "R4_EXPORT_ROOT",
  unset = file.path(project_root, "outputs", "result4_supplementary_FINAL_20260910")
)
source_out <- file.path(out_dir, "source")
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
cb_spline_grid <- bundle$cb_spline_grid
spline_names <- bundle$spline_names
spline_int_names <- bundle$spline_int_names
fit_spline_alt <- bundle$fit_spline_alt

extract_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, vcov = as.matrix(vcov(fit)$cond))
}

cumulative_basis <- function(values, centre = 0) {
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

curve_design <- function(values, altitude_1000_c, fit) {
  basis <- cumulative_basis(values, 0)
  fixed <- extract_fixed(fit)
  design <- matrix(
    0, nrow(basis), length(fixed$beta),
    dimnames = list(rownames(basis), names(fixed$beta))
  )
  design[, spline_names] <- basis
  design[, spline_int_names] <- basis * altitude_1000_c
  design
}

estimate_curve <- function(values, altitude_1000_c, fit) {
  design <- curve_design(values, altitude_1000_c, fit)
  fixed <- extract_fixed(fit)
  eta <- as.numeric(design %*% fixed$beta)
  tibble(anomaly = as.numeric(rownames(design)), estimate = exp(eta))
}

full_turning <- read_csv(
  file.path(source_dir, "confirmatory_full_support_turning_points.csv"),
  show_col_types = FALSE
)
full_shift_summary <- read_csv(
  file.path(source_dir, "confirmatory_full_support_peak_shift.csv"),
  show_col_types = FALSE
)
loyo_interaction <- read_csv(
  file.path(source_dir, "confirmatory_leave_one_year_out_interaction.csv"),
  show_col_types = FALSE
)
model_tests <- read_csv(
  file.path(source_dir, "confirmatory_model_tests.csv"),
  show_col_types = FALSE
)
gee_test <- read_csv(
  file.path(source_dir, "confirmatory_twoway_fixed_effects_gee.csv"),
  show_col_types = FALSE
)
humidity_test <- read_csv(
  file.path(source_dir, "confirmatory_absolute_humidity_adjustment.csv"),
  show_col_types = FALSE
)
gradient <- read_csv(gradient_path, show_col_types = FALSE)

# Reconstruct the common P10-P90 support used for the confirmatory maxima.
curve_grid <- sort(unique(bundle$curves_confirm$anomaly))
anchor_specs <- bundle$curves_confirm |>
  distinct(altitude_anchor, altitude, altitude_1000_c) |>
  arrange(altitude)

# Full-support maxima from each leave-one-year-out refit. These are point
# estimates from the same model and support used in the main confirmatory fit.
loyo_full_peaks <- bind_rows(lapply(2017:2024, function(omitted) {
  fit <- readRDS(file.path(
    checkpoint_dir, sprintf("loyo_%d_spline_alt.rds", omitted)
  ))
  bind_rows(lapply(seq_len(nrow(anchor_specs)), function(i) {
    spec <- anchor_specs[i, ]
    curve <- estimate_curve(curve_grid, spec$altitude_1000_c, fit)
    top <- curve[which.max(curve$estimate), ]
    tibble(
      omitted_year = omitted,
      altitude_anchor = spec$altitude_anchor,
      altitude = spec$altitude,
      peak_anomaly = top$anomaly,
      peak_or = top$estimate,
      at_lower_boundary = top$anomaly == min(curve_grid),
      at_upper_boundary = top$anomaly == max(curve_grid)
    )
  }))
}))

# Reproduce the coefficient-draw distribution of the altitude-related shift in
# the full-support response maximum.
set.seed(20260831)
n_draws <- 4000L
fixed_full <- extract_fixed(fit_spline_alt)
beta_draws <- MASS::mvrnorm(
  n = n_draws, mu = fixed_full$beta, Sigma = fixed_full$vcov
)
colnames(beta_draws) <- names(fixed_full$beta)

peak_matrix <- matrix(
  NA_real_, nrow(anchor_specs), n_draws,
  dimnames = list(anchor_specs$altitude_anchor, NULL)
)
for (i in seq_len(nrow(anchor_specs))) {
  spec <- anchor_specs[i, ]
  design <- curve_design(curve_grid, spec$altitude_1000_c, fit_spline_alt)
  eta <- design %*% t(beta_draws)
  peak_index <- max.col(t(eta), ties.method = "first")
  peak_matrix[i, ] <- curve_grid[peak_index]
}
alt_km <- anchor_specs$altitude / 1000
alt_centered <- alt_km - mean(alt_km)
shift_draws <- colSums(peak_matrix * alt_centered) / sum(alt_centered^2)
shift_data <- tibble(draw = seq_along(shift_draws), shift_c_per_1000m = shift_draws)
shift_check <- shift_data |>
  summarise(
    median = median(shift_c_per_1000m),
    low = quantile(shift_c_per_1000m, 0.025),
    high = quantile(shift_c_per_1000m, 0.975),
    probability_warmer = mean(shift_c_per_1000m > 0)
  )

# Coefficient-draw intervals for the peak shift in each leave-one-year-out
# refit. This is the direct stability estimand; showing these intervals avoids
# over-emphasizing individual peak coordinates from the three altitude anchors.
simulate_shift_for_fit <- function(fit, omitted_year, n = 2000L) {
  fixed <- extract_fixed(fit)
  set.seed(20260902 + omitted_year)
  draws <- MASS::mvrnorm(n = n, mu = fixed$beta, Sigma = fixed$vcov)
  colnames(draws) <- names(fixed$beta)
  peak_values <- matrix(NA_real_, nrow(anchor_specs), n)
  boundary_draw <- matrix(FALSE, nrow(anchor_specs), n)
  point_peaks <- numeric(nrow(anchor_specs))
  for (i in seq_len(nrow(anchor_specs))) {
    spec <- anchor_specs[i, ]
    design <- curve_design(curve_grid, spec$altitude_1000_c, fit)
    eta_draw <- design %*% t(draws)
    peak_index <- max.col(t(eta_draw), ties.method = "first")
    peak_values[i, ] <- curve_grid[peak_index]
    boundary_draw[i, ] <- peak_index == 1L | peak_index == length(curve_grid)
    point_curve <- estimate_curve(curve_grid, spec$altitude_1000_c, fit)
    point_peaks[i] <- point_curve$anomaly[which.max(point_curve$estimate)]
  }
  shift_draw <- colSums(peak_values * alt_centered) / sum(alt_centered^2)
  point_shift <- sum(point_peaks * alt_centered) / sum(alt_centered^2)
  tibble(
    omitted_year = omitted_year,
    point_estimate = point_shift,
    median = median(shift_draw),
    low = quantile(shift_draw, 0.025),
    high = quantile(shift_draw, 0.975),
    probability_warmer = mean(shift_draw > 0),
    probability_any_boundary = mean(colSums(boundary_draw) > 0),
    draws = n
  )
}

loyo_shift_summary <- bind_rows(lapply(2017:2024, function(omitted) {
  fit <- readRDS(file.path(
    checkpoint_dir, sprintf("loyo_%d_spline_alt.rds", omitted)
  ))
  simulate_shift_for_fit(fit, omitted)
}))

global_interaction_full <- model_tests |>
  filter(test == "Global spline-by-altitude interaction") |>
  transmute(
    sample = "Full sample",
    omitted_year = NA_integer_,
    likelihood_ratio,
    df = df_difference,
    p_value
  )
global_interaction_loyo <- loyo_interaction |>
  transmute(
    sample = paste("Omit", omitted_year),
    omitted_year,
    likelihood_ratio,
    df = df_difference,
    p_value
  )
interaction_summary <- bind_rows(
  global_interaction_full,
  global_interaction_loyo
)

write_csv(loyo_full_peaks, file.path(source_out, "figure_s3_loyo_full_support_peaks.csv"))
write_csv(shift_data, file.path(source_out, "figure_s3_peak_shift_draws.csv"))
write_csv(shift_check, file.path(source_out, "figure_s3_peak_shift_summary.csv"))
write_csv(loyo_shift_summary, file.path(source_out, "figure_s3_loyo_peak_shift_summary.csv"))
write_csv(interaction_summary, file.path(source_out, "figure_s3_global_interaction_summary.csv"))

# Data for Supplementary Table S3.
table_turning <- full_turning |>
  transmute(
    altitude_anchor, altitude_m = altitude,
    peak_anomaly_c = median_peak_anomaly,
    peak_anomaly_low_c = peak_anomaly_low,
    peak_anomaly_high_c = peak_anomaly_high,
    peak_or = median_peak_or,
    peak_or_low, peak_or_high,
    probability_internal_peak,
    probability_peak_below_zero,
    probability_peak_above_reference
  )

global_interaction <- model_tests |>
  filter(test == "Global spline-by-altitude interaction")
linear_interaction <- model_tests |>
  filter(test == "Linear anomaly-by-altitude interaction")

table_robustness <- bind_rows(
  tibble(
    analysis = "Hierarchical nonlinear anomaly-by-altitude interaction",
    estimate = NA_real_, low = NA_real_, high = NA_real_,
    statistic = global_interaction$likelihood_ratio,
    df = global_interaction$df_difference,
    p_value = global_interaction$p_value,
    details = "Likelihood-ratio test"
  ),
  tibble(
    analysis = "Linear anomaly-by-altitude interaction",
    estimate = NA_real_, low = NA_real_, high = NA_real_,
    statistic = linear_interaction$likelihood_ratio,
    df = linear_interaction$df_difference,
    p_value = linear_interaction$p_value,
    details = "Likelihood-ratio test"
  ),
  tibble(
    analysis = "Two-way fixed-effects GEE nonlinear interaction",
    estimate = NA_real_, low = NA_real_, high = NA_real_,
    statistic = gee_test$statistic,
    df = gee_test$df,
    p_value = gee_test$p_value,
    details = "County and province-week fixed effects"
  ),
  tibble(
    analysis = "Absolute-humidity-adjusted hierarchical interaction",
    estimate = NA_real_, low = NA_real_, high = NA_real_,
    statistic = humidity_test$likelihood_ratio,
    df = humidity_test$df_difference,
    p_value = humidity_test$p_value,
    details = "Likelihood-ratio test"
  ),
  tibble(
    analysis = "Leave-one-year-out interaction",
    estimate = NA_real_, low = min(loyo_interaction$p_value),
    high = max(loyo_interaction$p_value),
    statistic = sum(loyo_interaction$p_value < 0.01),
    df = nrow(loyo_interaction), p_value = NA_real_,
    details = "Folds with P<0.01 / total; estimate interval gives P-value range"
  ),
  tibble(
    analysis = "Shift in full-support response maximum per 1,000 m",
    estimate = shift_check$median,
    low = shift_check$low, high = shift_check$high,
    statistic = shift_check$probability_warmer,
    df = n_draws, p_value = NA_real_,
    details = "Coefficient-draw median (95% interval); statistic is Pr(shift>0)"
  )
)

write_csv(gradient, file.path(source_out, "table_s3_temperature_altitude_gradient.csv"))
write_csv(table_turning, file.path(source_out, "table_s3_response_maxima.csv"))
write_csv(table_robustness, file.path(source_out, "table_s3_model_robustness.csv"))

# Figure styling.
blue <- "#2F6FAE"
blue_dark <- "#184E81"
blue_light <- "#DCEAF5"
gray_line <- "#A8B2BC"
ink <- "#162433"

theme_supp <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 9.7, color = ink),
    plot.subtitle = element_text(size = 7.2, color = "#657586", lineheight = 1.02),
    axis.title = element_text(size = 8.8, color = ink),
    axis.text = element_text(size = 7.8, color = ink),
    axis.line = element_line(linewidth = 0.45, color = ink),
    axis.ticks = element_line(linewidth = 0.4, color = ink),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "#E5EBF0"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 7, 4, 5)
  )

display_order <- c("Full sample", paste("Omit", 2017:2024))
scenario_levels <- rev(display_order)
interaction_plot <- interaction_summary |>
  mutate(
    scenario = factor(sample, levels = scenario_levels),
    is_full = sample == "Full sample"
  )
interaction_full_plot <- interaction_plot |> filter(is_full)
interaction_loyo_plot <- interaction_plot |> filter(!is_full)
interaction_threshold <- qchisq(0.95, df = 6)

p_a <- ggplot() +
  geom_vline(
    xintercept = interaction_threshold,
    linetype = "dashed", linewidth = 0.5, color = "#7C8996"
  ) +
  geom_segment(
    data = interaction_loyo_plot,
    aes(x = interaction_threshold, xend = likelihood_ratio,
        y = scenario, yend = scenario),
    linewidth = 0.7, color = "#C2CBD3", lineend = "round"
  ) +
  geom_point(
    data = interaction_loyo_plot,
    aes(likelihood_ratio, scenario),
    shape = 21, size = 2.1, stroke = 0.55,
    color = "#8997A5", fill = "white"
  ) +
  geom_segment(
    data = interaction_full_plot,
    aes(x = interaction_threshold, xend = likelihood_ratio,
        y = scenario, yend = scenario),
    linewidth = 1.25, color = blue_dark, lineend = "round"
  ) +
  geom_point(
    data = interaction_full_plot,
    aes(likelihood_ratio, scenario),
    shape = 21, size = 3.0, stroke = 0.9,
    color = blue_dark, fill = "white"
  ) +
  scale_x_continuous(
    limits = c(0, max(interaction_plot$likelihood_ratio) * 1.08),
    breaks = c(0, interaction_threshold, 25, 35, 45),
    labels = c("0", "12.6", "25", "35", "45")
  ) +
  labs(
    title = "a  Global interaction across yearly refits",
    subtitle = "LRT; dashed line, P = 0.05 (6 df); all refits P < 0.01",
    x = expression(paste("Likelihood-ratio ", chi^2, " (6 df)")),
    y = NULL
  ) +
  theme_supp +
  theme(
    axis.text.y = element_text(size = 7.4),
    panel.grid.major.y = element_line(linewidth = 0.35, color = "#E7ECF0"),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5, 7, 4, 5)
  )

shift_full_plot <- tibble(
  scenario = factor("Full sample", levels = scenario_levels),
  median = shift_check$median,
  low = shift_check$low,
  high = shift_check$high,
  probability_warmer = shift_check$probability_warmer
)
shift_loyo_plot <- loyo_shift_summary |>
  mutate(
    scenario = factor(paste("Omit", omitted_year), levels = scenario_levels),
    boundary_sensitive = probability_any_boundary > 0.25
  )

p_b <- ggplot() +
  geom_vline(
    xintercept = 0, linetype = "dashed", linewidth = 0.5,
    color = "#7C8996"
  ) +
  geom_segment(
    data = shift_loyo_plot,
    aes(x = low, xend = high, y = scenario, yend = scenario),
    linewidth = 0.72, color = gray_line, lineend = "round"
  ) +
  geom_point(
    data = shift_loyo_plot,
    aes(median, scenario, shape = boundary_sensitive),
    size = 2.05, stroke = 0.55,
    color = "#8997A5", fill = "white"
  ) +
  geom_segment(
    data = shift_full_plot,
    aes(x = low, xend = high, y = scenario, yend = scenario),
    linewidth = 1.3, color = blue_dark, lineend = "round"
  ) +
  geom_point(
    data = shift_full_plot,
    aes(median, scenario),
    shape = 21, size = 2.9, stroke = 0.9,
    color = blue_dark, fill = "white"
  ) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.06))) +
  labs(
    title = "b  Peak displacement across yearly refits",
    subtitle = "Positive = warmer maximum at higher altitude",
    x = expression(paste("Peak shift per 1,000 m (", degree, "C)")),
    y = NULL,
    caption = "Bars, 95% intervals; triangles, boundary-sensitive refits"
  ) +
  theme_supp +
  theme(
    axis.text.y = element_text(size = 7.4),
    panel.grid.major.y = element_line(linewidth = 0.35, color = "#E7ECF0"),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(
      family = "Arial", size = 6.7, color = "#657586", hjust = 0,
      margin = margin(t = 3)
    ),
    plot.margin = margin(5, 5, 2, 7)
  )

degree_c_text <- paste0(intToUtf8(176), "C")

figure <- p_a + p_b +
  plot_layout(widths = c(0.95, 1.05)) +
  plot_annotation(
    title = "Leave-one-year-out robustness of altitude-related response features",
    subtitle = paste0(
      "Full sample and eight leave-one-year-out refits; cumulative lag 0-2 weeks; ",
      "common P10-P90 anomaly support: ",
      sprintf("%.2f to %.2f ", min(curve_grid), max(curve_grid)),
      degree_c_text
    ),
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 11.5, color = ink,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        family = "Arial", size = 7.4, color = "#657586",
        margin = margin(b = 7)
      ),
      plot.margin = margin(7, 7, 5, 7)
    )
  )

stem <- file.path(out_dir, "Supplementary_Figure_S3_Result4_Robustness_FINAL_20260910")
png_path <- paste0(stem, ".png")
pdf_path <- paste0(stem, ".pdf")
svg_path <- paste0(stem, ".svg")
tiff_path <- paste0(stem, ".tiff")
ggsave(
  png_path, figure, width = 183, height = 104, units = "mm", dpi = 600,
  device = grDevices::png, type = "cairo-png", bg = "white"
)
ggsave(pdf_path, figure, width = 183, height = 104, units = "mm", device = cairo_pdf, bg = "white")
ggsave(svg_path, figure, width = 183, height = 104, units = "mm", device = svg, bg = "white")
ggsave(
  tiff_path, figure, width = 183, height = 104, units = "mm", dpi = 600,
  device = grDevices::tiff, compression = "lzw", bg = "white"
)

# Figure legend and machine-checkable QA.
legend_text <- paste(
  "Supplementary Figure S3. Leave-one-year-out robustness of altitude-related response features.",
  "a, Likelihood-ratio statistics for the global nonlinear temperature-anomaly-by-altitude interaction in the full sample and after omitting each calendar year in turn. The dashed line is the chi-square threshold for P = 0.05 with 6 degrees of freedom; every refit had P < 0.01.",
  "b, Change in the fitted response maximum per 1,000-m increase in altitude; positive values indicate a warmer maximum at higher altitude. Points and horizontal bars show medians and 95% fixed-coefficient intervals. Triangles identify boundary-sensitive refits.",
  "All maxima were searched only within the anomaly range jointly supported between the P10 and P90 observed distributions (−2.58 to 2.64 °C).",
  "A maximum on a support boundary indicates that its location is weakly identified and should not be interpreted as a confirmed biological threshold.",
  "The simulations quantify fixed-coefficient uncertainty conditional on each fitted hierarchical model and do not include model-selection uncertainty."
)
writeLines(legend_text, file.path(out_dir, "Supplementary_Figure_S3_legend.txt"))

qa <- tibble(
  check = c(
    "All eight LOYO fits represented at three altitude anchors",
    "All eight LOYO peak-shift summaries are finite",
    "All LOYO full models converged with positive-definite Hessian",
    "All LOYO maxima lie within common support",
    "Reconstructed shift median matches saved result",
    "Reconstructed shift interval matches saved result",
    "All LOYO global interaction tests have P < 0.01",
    "PNG export exists",
    "PDF export exists",
    "SVG export exists",
    "TIFF export exists"
  ),
  passed = c(
    nrow(loyo_full_peaks) == 24L && n_distinct(loyo_full_peaks$omitted_year) == 8L,
    nrow(loyo_shift_summary) == 8L &&
      all(is.finite(unlist(loyo_shift_summary[c("median", "low", "high")]))),
    all(loyo_interaction$convergence_full == 0) && all(loyo_interaction$pdHess_full),
    all(loyo_full_peaks$peak_anomaly >= min(curve_grid)) &&
      all(loyo_full_peaks$peak_anomaly <= max(curve_grid)),
    abs(shift_check$median - full_shift_summary$median) < 1e-8,
    abs(shift_check$low - full_shift_summary$low) < 1e-8 &&
      abs(shift_check$high - full_shift_summary$high) < 1e-8,
    all(loyo_interaction$p_value < 0.01),
    file.exists(png_path), file.exists(pdf_path),
    file.exists(svg_path), file.exists(tiff_path)
  ),
  details = c(
    sprintf("%d rows; %d omitted years", nrow(loyo_full_peaks), n_distinct(loyo_full_peaks$omitted_year)),
    sprintf("%d finite LOYO shift summaries", nrow(loyo_shift_summary)),
    sprintf("%d/%d converged and pdHess", sum(loyo_interaction$convergence_full == 0 & loyo_interaction$pdHess_full), nrow(loyo_interaction)),
    paste0(
      sprintf("support %.2f to %.2f ", min(curve_grid), max(curve_grid)),
      degree_c_text
    ),
    sprintf("reconstructed %.9f; saved %.9f", shift_check$median, full_shift_summary$median),
    sprintf("reconstructed %.9f to %.9f; saved %.9f to %.9f", shift_check$low, shift_check$high, full_shift_summary$low, full_shift_summary$high),
    sprintf("P range %.3g to %.3g", min(loyo_interaction$p_value), max(loyo_interaction$p_value)),
    png_path, pdf_path, svg_path, tiff_path
  )
)
write_csv(qa, file.path(out_dir, "Supplementary_Figure_S3_QA.csv"))
if (!all(qa$passed)) stop("One or more Supplementary Figure S3 QA checks failed.")

message("Created: ", png_path)
message("Created: ", pdf_path)
message("Created: ", svg_path)
message("Created: ", tiff_path)
