suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
  library(dlnm)
  library(glmmTMB)
})

# Figure contract
# Core conclusion: the local temperature-anomaly response is nonlinear, but
# the fitted maximum shifts from a negative anomaly at low altitude to a
# positive anomaly at high altitude, with substantial location uncertainty.
# Archetype: asymmetric quantitative grid; panel b is the hero evidence.
# Output: 183-mm double-column figure with editable SVG/PDF and 600-dpi TIFF.

analysis_id <- "confirmatory_anomaly_shape_20260831"
checkpoint_dir <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", analysis_id
)
source_dir <- file.path(
  "result0715", "source_data", paste0("v4_", analysis_id)
)
context_source_dir <- file.path(
  "result0715", "source_data",
  "v4_result4_climate_context_phase_20260826"
)
output_dir <- file.path(
  "result0715", "main_text_figures", "v4_t2m_20260825",
  "result4_confirmatory"
)
qa_dir <- file.path(
  "result0715", "figure_qa", "v4_t2m_20260825",
  "result4_confirmatory"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(file.path(checkpoint_dir, "confirmatory_anomaly_bundle.rds"))
support <- read.csv(
  file.path(context_source_dir, "temperature_anomaly_support.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
gradient <- read.csv(
  file.path(context_source_dir, "temperature_altitude_gradient_summary.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
turning <- read.csv(
  file.path(source_dir, "confirmatory_full_support_turning_points.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
model_tests <- read.csv(
  file.path(source_dir, "confirmatory_model_tests.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
gee_test <- read.csv(
  file.path(source_dir, "confirmatory_twoway_fixed_effects_gee.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

fit <- bundle$fit_spline_alt
cb <- bundle$cb_spline_grid
spline_names <- bundle$spline_names
interaction_names <- bundle$spline_int_names

extract_fixed <- function(model) {
  list(beta = fixef(model)$cond, vcov = as.matrix(vcov(model)$cond))
}

cumulative_basis <- function(values, centre = 0) {
  values <- sort(unique(c(values, centre)))
  original_names <- colnames(cb)
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
      cb, coef = coefficient, vcov = zero_vcov,
      model.link = "logit", at = values, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

estimate_curve <- function(values, altitude_1000_c) {
  basis <- cumulative_basis(values, 0)
  fixed <- extract_fixed(fit)
  design <- matrix(
    0, nrow(basis), length(fixed$beta),
    dimnames = list(rownames(basis), names(fixed$beta))
  )
  design[, spline_names] <- basis
  design[, interaction_names] <- basis * altitude_1000_c
  log_est <- as.numeric(design %*% fixed$beta)
  se <- sqrt(pmax(rowSums((design %*% fixed$vcov) * design), 0))
  tibble(
    anomaly = as.numeric(rownames(design)),
    estimate = exp(log_est),
    low = exp(log_est - 1.96 * se),
    high = exp(log_est + 1.96 * se)
  )
}

p05_low <- max(support$p05)
p95_high <- min(support$p95)
p10_low <- max(support$p10)
p90_high <- min(support$p90)
plot_low <- ceiling(p05_low * 10) / 10
plot_high <- floor(p95_high * 10) / 10
values <- seq(plot_low, plot_high, by = 0.02)

anchor_specs <- bundle$curves_confirm |>
  distinct(altitude_anchor, altitude, altitude_1000_c) |>
  arrange(factor(altitude_anchor, levels = c("P10", "P50", "P90"))) |>
  mutate(
    altitude_anchor = factor(altitude_anchor, levels = c("P10", "P50", "P90")),
    altitude_label = factor(
      sprintf("%s m", format(round(altitude), big.mark = ",", scientific = FALSE)),
      levels = sprintf(
        "%s m", format(round(altitude), big.mark = ",", scientific = FALSE)
      )
    )
  )

curves <- bind_rows(lapply(seq_len(nrow(anchor_specs)), function(i) {
  spec <- anchor_specs[i, ]
  estimate_curve(values, spec$altitude_1000_c) |>
    mutate(
      altitude_anchor = spec$altitude_anchor,
      altitude_label = spec$altitude_label,
      altitude = spec$altitude,
      support_region = case_when(
        anomaly < p10_low ~ "Lower-support extension",
        anomaly > p90_high ~ "Upper-support extension",
        TRUE ~ "Common P10-P90"
      ),
      support_segment = case_when(
        anomaly < p10_low ~ "cold_edge",
        anomaly > p90_high ~ "warm_edge",
        TRUE ~ "central"
      )
    )
}))
common_support_upper <- max(
  curves$anomaly[curves$support_region == "Common P10-P90"]
)

anchor_colours <- c(P10 = "#78A8CC", P50 = "#3F82B8", P90 = "#1F5F96")
ink <- "#17202A"
muted <- "#6D7C8B"
grid <- "#D9E0E6"
font_family <- "Microsoft YaHei"

turning_plot <- turning |>
  mutate(
    altitude_anchor = factor(altitude_anchor, levels = c("P10", "P50", "P90")),
    altitude_label = factor(
      sprintf("%s m", format(round(altitude), big.mark = ",", scientific = FALSE)),
      levels = rev(sprintf(
        "%s m", format(round(anchor_specs$altitude), big.mark = ",", scientific = FALSE)
      ))
    ),
    point_colour = unname(anchor_colours[as.character(altitude_anchor)]),
    probability_percent = 100 * probability_peak_below_zero,
    probability_label = sprintf("%.1f%%", probability_percent),
    location_label = sprintf(
      "%.2f (%.2f to %.2f) \u00B0C",
      median_peak_anomaly, peak_anomaly_low, peak_anomaly_high
    )
  )

peak_points <- turning_plot |>
  select(altitude_anchor, altitude_label, altitude, median_peak_anomaly) |>
  left_join(
    curves |>
      group_by(altitude_anchor, altitude_label, altitude) |>
      slice_min(abs(anomaly - first(
        turning_plot$median_peak_anomaly[
          turning_plot$altitude_anchor == first(altitude_anchor)
        ]
      )), n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(
        altitude_anchor, altitude_label, altitude,
        peak_x = anomaly, peak_y = estimate
      ),
    by = c("altitude_anchor", "altitude_label", "altitude")
  )

theme_manuscript <- theme_classic(base_size = 7.2, base_family = font_family) +
  theme(
    plot.title = element_text(
      face = "bold", colour = ink, size = 8.4, margin = margin(b = 1.5)
    ),
    plot.subtitle = element_text(
      colour = muted, size = 6.45, lineheight = 1.05, margin = margin(b = 4)
    ),
    axis.title = element_text(colour = ink, size = 6.9),
    axis.text = element_text(colour = ink, size = 6.25),
    axis.line = element_line(colour = ink, linewidth = 0.36),
    axis.ticks = element_line(colour = ink, linewidth = 0.28),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = ink, size = 6.8),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(4, 5, 4, 5)
  )

annual_gradient <- gradient |> filter(period_type == "Annual")
monthly_gradient <- gradient |>
  filter(period_type == "Month") |>
  mutate(period = factor(period, levels = month.abb))

panel_a <- ggplot(monthly_gradient, aes(period, estimate, group = 1)) +
  geom_hline(
    yintercept = annual_gradient$estimate,
    colour = muted, linewidth = 0.42, linetype = "22"
  ) +
  geom_errorbar(
    aes(ymin = low, ymax = high), width = 0,
    colour = anchor_colours[["P50"]], linewidth = 0.42
  ) +
  geom_line(colour = anchor_colours[["P50"]], linewidth = 0.56) +
  geom_point(
    shape = 21, fill = "white", colour = anchor_colours[["P50"]],
    size = 1.65, stroke = 0.4
  ) +
  scale_y_continuous(breaks = seq(-6.5, -3.5, 0.5)) +
  scale_x_discrete(breaks = month.abb[c(1, 3, 5, 7, 9, 11)]) +
  labs(
    title = "Empirical air-temperature gradient",
    subtitle = sprintf(
      "Monthly slopes (95%% CI); annual %.2f \u00B0C/1,000 m",
      annual_gradient$estimate
    ),
    x = NULL, y = "Gradient (\u00B0C per 1,000 m)"
  ) +
  coord_cartesian(ylim = c(-6.75, -3.45), clip = "off") +
  theme_manuscript +
  theme(axis.text.x = element_text(size = 5.85))

panel_b <- ggplot(curves, aes(anomaly, estimate, colour = altitude_anchor)) +
  geom_hline(yintercept = 1, colour = muted, linewidth = 0.36, linetype = "22") +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.34, linetype = "22") +
  geom_ribbon(
    data = filter(curves, support_region == "Common P10-P90"),
    aes(x = anomaly, ymin = low, ymax = high, fill = altitude_anchor),
    inherit.aes = FALSE, alpha = 0.18, colour = NA
  ) +
  geom_ribbon(
    data = filter(curves, support_region != "Common P10-P90"),
    aes(x = anomaly, ymin = low, ymax = high, fill = altitude_anchor,
        group = interaction(altitude_anchor, support_segment)),
    inherit.aes = FALSE, alpha = 0.07, colour = NA
  ) +
  geom_line(
    data = filter(curves, support_region == "Common P10-P90"),
    linewidth = 0.72
  ) +
  geom_line(
    data = filter(curves, support_region != "Common P10-P90"),
    aes(group = interaction(altitude_anchor, support_segment)),
    linewidth = 0.62, linetype = "22", alpha = 0.8
  ) +
  geom_point(
    data = peak_points,
    aes(x = peak_x, y = peak_y, colour = altitude_anchor),
    inherit.aes = FALSE, shape = 21, fill = "white",
    size = 1.8, stroke = 0.55
  ) +
  facet_wrap(~ altitude_label, nrow = 1) +
  scale_colour_manual(values = anchor_colours) +
  scale_fill_manual(values = anchor_colours) +
  scale_x_continuous(breaks = c(-3, -2, 0, 2, 3)) +
  scale_y_log10(
    breaks = c(0.6, 0.8, 1, 1.2, 1.5),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Temperature-anomaly responses",
    subtitle = "Solid: common support; dashed: extended support; points: maxima",
    x = "Temperature anomaly (\u00B0C)", y = "Cumulative OR"
  ) +
  coord_cartesian(xlim = c(plot_low, plot_high), ylim = c(0.58, 1.55)) +
  theme_manuscript +
  theme(
    panel.spacing.x = unit(6, "pt"),
    panel.grid.major.y = element_line(colour = grid, linewidth = 0.22)
  )

global_p <- model_tests$p_value[
  model_tests$test == "Global spline-by-altitude interaction"
]
gee_p <- gee_test$p_value[[1]]

panel_c <- ggplot(turning_plot, aes(median_peak_anomaly, altitude_label)) +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.42, linetype = "22") +
  geom_vline(
    xintercept = common_support_upper, colour = "#B8C3CC",
    linewidth = 0.40, linetype = "33"
  ) +
  geom_errorbar(
    aes(
      xmin = peak_anomaly_low, xmax = peak_anomaly_high,
      colour = point_colour
    ),
    orientation = "y", width = 0, linewidth = 0.72
  ) +
  geom_point(
    aes(colour = point_colour), shape = 21, fill = "white",
    size = 2.35, stroke = 0.72
  ) +
  geom_text(
    aes(x = 2.80, label = location_label),
    hjust = 0, colour = ink, family = font_family, size = 2.05
  ) +
  scale_colour_identity() +
  scale_x_continuous(
    breaks = c(-1, 0, 1, 2, common_support_upper),
    labels = c(
      "-1", "0", "1", "2", sprintf("%.2f", common_support_upper)
    )
  ) +
  labs(
    title = "Estimated response maximum",
    subtitle = "Median (95% simulation interval); dotted: support limit",
    x = "Temperature anomaly at maximum (\u00B0C)",
    y = NULL
  ) +
  coord_cartesian(xlim = c(-1.5, 4.15), clip = "off") +
  theme_manuscript +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.23),
    plot.margin = margin(4, 45, 4, 5)
  )

panel_d <- ggplot(
  turning_plot,
  aes(probability_percent, altitude_label, fill = point_colour)
) +
  geom_col(width = 0.34, alpha = 0.88) +
  geom_text(
    aes(
      x = pmin(probability_percent + 4, 105),
      label = probability_label
    ),
    hjust = 0, colour = ink, family = font_family, size = 2.1
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = c(0, 50, 100),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Probability of a negative-anomaly peak",
    subtitle = "4,000 simulations",
    x = NULL, y = NULL
  ) +
  coord_cartesian(xlim = c(0, 113), clip = "off") +
  theme_manuscript +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.23),
    axis.text.y = element_text(size = 5.9),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 14, 4, 5)
  )

top_row <- panel_a + panel_b + plot_layout(widths = c(0.75, 1.45))
bottom_row <- panel_c + panel_d + plot_layout(widths = c(1.18, 0.82))
combined <- (top_row / bottom_row) +
  plot_layout(heights = c(1.16, 0.84)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        family = font_family, face = "bold", colour = ink, size = 8.4
      ),
      plot.margin = margin(6, 8, 5, 7)
    )
  )

base_path <- file.path(
  output_dir, "Figure5_AGU_no_title_20261008"
)
width_in <- 183 / 25.4
height_in <- 126 / 25.4

svglite::svglite(paste0(base_path, ".svg"), width = width_in, height = height_in)
print(combined)
grDevices::dev.off()

grDevices::cairo_pdf(
  paste0(base_path, ".pdf"), width = width_in, height = height_in,
  family = font_family, bg = "white"
)
print(combined)
grDevices::dev.off()

grDevices::png(
  paste0(base_path, ".png"), width = width_in, height = height_in,
  units = "in", res = 450, type = "cairo-png", bg = "white"
)
print(combined)
grDevices::dev.off()

grDevices::tiff(
  paste0(base_path, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw", type = "cairo", bg = "white"
)
print(combined)
grDevices::dev.off()

legend <- paste0(
  "Figure 5. Local thermal context and temperature-anomaly responses across altitude. ",
  "Laboratory positivity was modelled against weekly mean 2-m air-temperature anomaly ",
  "using cumulative lag 0-2 weeks during 2017-2024. Local climatology was defined using ",
  "a leave-one-year-out approach, and the spring 2023 epidemic was retained. ",
  "a, Monthly empirical air-temperature gradients with 95% confidence intervals; the ",
  "dashed line denotes the annual gradient. b, Temperature-anomaly response curves at ",
  "2,150, 2,980 and 4,000 m. Solid curves show common support, dashed extensions show ",
  "extended support, and points mark fitted maxima; shading denotes 95% confidence intervals. ",
  "c, Estimated temperature anomaly at the fitted response maximum, shown as simulation ",
  "medians and 95% intervals. d, Proportion of simulations in which the fitted maximum ",
  "occurred below 0 \u00B0C."
)
writeLines(
  enc2utf8(legend),
  file.path(output_dir, "Figure_5_AGU_legend_20261008.txt"),
  useBytes = TRUE
)

write.csv(
  curves,
  file.path(source_dir, "figure_result4_panel_b_confirmatory_broad_curves.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  turning_plot,
  file.path(source_dir, "figure_result4_panel_c_parametric_simulation_maxima.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  turning_plot |>
    select(
      altitude_anchor, altitude, probability_internal_peak,
      probability_peak_below_zero, probability_percent, probability_label
    ),
  file.path(source_dir, "figure_result4_panel_d_peak_side_simulation_proportion.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

qa <- tibble(
  check = c(
    "Panel b is based on the unpenalized confirmatory model",
    "Solid curve range is restricted to common P10-P90 support",
    "Dashed extensions remain within common P05-P95 support",
    "All three altitude anchors are present",
    "High-altitude fitted maximum is above zero anomaly",
    "Panel d negative-side probabilities are valid",
    "High-altitude negative-side maximum probability is below 50 percent",
    "All confidence limits are finite and positive",
    "All four export formats exist",
    "Raster exports exceed minimum usable file size"
  ),
  pass = c(
    is.null(fit$modelInfo$priors) || nrow(fit$modelInfo$priors) == 0,
    min(filter(curves, support_region == "Common P10-P90")$anomaly) >= p10_low - 0.03 &&
      max(filter(curves, support_region == "Common P10-P90")$anomaly) <= p90_high + 0.03,
    min(curves$anomaly) >= p05_low - 0.03 && max(curves$anomaly) <= p95_high + 0.03,
    n_distinct(curves$altitude_anchor) == 3,
    turning_plot$median_peak_anomaly[
      turning_plot$altitude_anchor == "P90"
    ] > 0,
    all(turning_plot$probability_peak_below_zero >= 0 &
          turning_plot$probability_peak_below_zero <= 1),
    turning_plot$probability_peak_below_zero[
      turning_plot$altitude_anchor == "P90"
    ] < 0.5,
    all(is.finite(c(curves$low, curves$high))) &&
      all(c(curves$low, curves$high) > 0),
    all(file.exists(paste0(base_path, c(".svg", ".pdf", ".png", ".tiff")))),
    all(file.info(paste0(base_path, c(".png", ".tiff")))$size > 100000)
  )
)
write.csv(
  qa,
  file.path(
    qa_dir,
    "Result4_T2m_local_thermal_context_REVISED_20260910_qa.csv"
  ),
  row.names = FALSE, fileEncoding = "UTF-8"
)
if (!all(qa$pass)) warning("At least one confirmatory Result 4 QA item failed.")

message("Created confirmatory final Result 4 figure.")
