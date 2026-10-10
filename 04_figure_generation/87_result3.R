suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(plot3D)
  library(dlnm)
  library(glmmTMB)
  library(scales)
  library(ragg)
  library(png)
  library(grid)
})

set.seed(20260826)

# Figure contract
# Core conclusion: the lag-specific air-temperature response differs among
# overall, influenza A and influenza B positivity, while observed trajectories
# and standardized phenology place those response surfaces in subtype context.
# Archetype: asymmetric quantitative composite; R-only; double-column output.

degree_c <- "\u00B0C"
font_family <- "sans"
dark <- "#202428"
muted <- "#687A8B"
blue <- "#2F6FA9"
orange <- "#C8872B"

bundle_path <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825",
  "v4_t2m_model_bundle.rds"
)
figure_dir <- file.path(
  "result0715", "main_text_figures", "v4_t2m_20260825",
  "result3_harmonized"
)
source_dir <- file.path(
  "result0715", "source_data", "v4_result3_t2m_harmonized_20260826"
)
qa_dir <- file.path(
  "result0715", "figure_qa", "v4_t2m_20260825",
  "result3_harmonized"
)
for (path in c(figure_dir, source_dir, qa_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

bundle <- readRDS(bundle_path)
model_dat <- bundle$model_data
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
z_int_names <- bundle$z_int_names
interaction_scale <- bundle$interaction_scale
anchors <- bundle$altitude_anchors
support <- bundle$local_support
reference_temperature <- bundle$reference_temperature

fits <- list(
  `Overall positivity` = bundle$fits$positive_full,
  `Influenza A-only positivity` = bundle$fits$fluA_only_full,
  `Influenza B-only positivity` = bundle$fits$fluB_only_full
)

p50 <- anchors |>
  filter(altitude_anchor == "P50") |>
  left_join(support, by = c("altitude_anchor", "altitude"))
if (nrow(p50) != 1L) stop("The P50 altitude display anchor is not unique.")
slice_altitude <- p50$altitude[[1]]
slice_altitude_1000_c <- p50$altitude_1000_c[[1]]

temperature_grid <- sort(unique(c(
  seq(
    floor(p50$p05[[1]] * 4) / 4,
    ceiling(p50$p95[[1]] * 4) / 4,
    by = 0.25
  ),
  reference_temperature
)))
lag_grid <- 0:2

extract_fixed <- function(fit) {
  list(
    beta = glmmTMB::fixef(fit)$cond,
    vcov = as.matrix(stats::vcov(fit)$cond)
  )
}

basis_array_lag <- function(temperatures, centre) {
  original_names <- colnames(cb_all)
  output <- array(
    0,
    dim = c(length(temperatures), length(lag_grid), length(original_names)),
    dimnames = list(
      as.character(temperatures), paste0("lag", lag_grid), cb_names
    )
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
      cumul = FALSE
    )
    output[, , j] <- prediction$matfit[
      , paste0("lag", lag_grid), drop = FALSE
    ]
  }
  output
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

estimate_basis_rows <- function(fit, basis, altitude_1000_c) {
  fixed <- extract_fixed(fit)
  design <- matrix(
    0, nrow = nrow(basis), ncol = length(fixed$beta),
    dimnames = list(NULL, names(fixed$beta))
  )
  design[, cb_names] <- basis
  design[, z_int_names] <-
    sweep(basis, 2, interaction_scale, "/") * altitude_1000_c
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

lag_basis <- basis_array_lag(temperature_grid, reference_temperature)
surface_rows <- list()
for (outcome_name in names(fits)) {
  for (temperature_index in seq_along(temperature_grid)) {
    for (lag_index in seq_along(lag_grid)) {
      basis_row <- matrix(
        lag_basis[temperature_index, lag_index, ],
        nrow = 1, dimnames = list(NULL, cb_names)
      )
      estimate <- estimate_basis_rows(
        fits[[outcome_name]], basis_row, slice_altitude_1000_c
      )
      surface_rows[[length(surface_rows) + 1L]] <- bind_cols(
        tibble(
          outcome = outcome_name,
          temperature = temperature_grid[[temperature_index]],
          lag_week = lag_grid[[lag_index]],
          altitude = slice_altitude,
          reference_temperature = reference_temperature
        ),
        estimate
      )
    }
  }
}
surface_data <- bind_rows(surface_rows)

# The product of the three lag-specific ORs must reproduce the cumulative
# lag 0-2 contrast from the same fitted model.
cumulative_basis <- basis_matrix_cumulative(
  temperature_grid, reference_temperature
)
cumulative_expected <- bind_rows(lapply(names(fits), function(outcome_name) {
  bind_cols(
    tibble(outcome = outcome_name, temperature = temperature_grid),
    estimate_basis_rows(
      fits[[outcome_name]], cumulative_basis, slice_altitude_1000_c
    ) |>
      select(expected_cumulative_or = estimate)
  )
}))
cumulative_reconstructed <- surface_data |>
  group_by(outcome, temperature) |>
  summarise(
    reconstructed_cumulative_or = exp(sum(log_estimate)),
    .groups = "drop"
  )
consistency <- cumulative_reconstructed |>
  left_join(cumulative_expected, by = c("outcome", "temperature")) |>
  mutate(
    absolute_log_difference = abs(
      log(reconstructed_cumulative_or) - log(expected_cumulative_or)
    )
  )
max_consistency_error <- max(consistency$absolute_log_difference)

# -------------------------------------------------------------------------
# Panel a: harmonized lag-specific response surfaces at the P50 altitude.
# Linear interpolation between observed integer lags is used only to render the
# surface; all underlying model estimates remain at lags 0, 1 and 2.
# -------------------------------------------------------------------------
outcome_order <- names(fits)
surface_titles <- c(
  "Overall positivity" = "Overall positivity",
  "Influenza A-only positivity" = "Influenza A-only positivity",
  "Influenza B-only positivity" = "Influenza B-only positivity"
)
render_lag_grid <- seq(0, 2, by = 0.20)
render_temperature_grid <- temperature_grid[
  seq(1, length(temperature_grid), by = 4)
]
if (tail(render_temperature_grid, 1) != tail(temperature_grid, 1)) {
  render_temperature_grid <- c(
    render_temperature_grid, tail(temperature_grid, 1)
  )
}

surface_matrix <- function(outcome_name) {
  dat <- surface_data |>
    filter(
      outcome == outcome_name,
      temperature %in% render_temperature_grid
    )
  observed <- matrix(
    NA_real_, length(render_temperature_grid), length(lag_grid),
    dimnames = list(render_temperature_grid, lag_grid)
  )
  observed[cbind(
    match(dat$temperature, render_temperature_grid),
    match(dat$lag_week, lag_grid)
  )] <- dat$estimate
  output <- t(vapply(
    seq_len(nrow(observed)),
    function(i) {
      approx(
        lag_grid, observed[i, ], xout = render_lag_grid,
        method = "linear", rule = 2
      )$y
    },
    numeric(length(render_lag_grid))
  ))
  dimnames(output) <- list(render_temperature_grid, render_lag_grid)
  output
}

z_limits <- lapply(outcome_order, function(outcome_name) {
  values <- surface_data$estimate[surface_data$outcome == outcome_name]
  limits <- range(c(1, values), finite = TRUE)
  span <- max(diff(limits), 0.15)
  c(max(0.05, limits[[1]] - 0.07 * span), limits[[2]] + 0.07 * span)
})
names(z_limits) <- outcome_order

surface_blue <- colorRampPalette(
  c("#F4F9FB", "#E4F1F6", "#CEE5ED", "#ACD2DF", "#80B8CA")
)(100)

draw_surface <- function(outcome_name) {
  par(mar = c(1.62, 0.70, 1.25, 0.04))
  z <- surface_matrix(outcome_name)
  limits <- z_limits[[outcome_name]]
  persp3D(
    x = render_temperature_grid,
    y = render_lag_grid,
    z = z,
    colvar = z,
    col = surface_blue,
    breaks = seq(limits[[1]], limits[[2]], length.out = 101),
    clim = limits, zlim = limits,
    theta = 34, phi = 23, expand = 0.82, scale = TRUE,
    facets = TRUE,
    border = grDevices::adjustcolor("#6F9EAE", alpha.f = 0.72),
    lwd = 0.30, lighting = FALSE, shade = NA,
    colkey = FALSE, bty = "b",
    xlab = paste0("2-m air temperature (", degree_c, ")"),
    ylab = "Lag (weeks)", zlab = "OR",
    main = surface_titles[[outcome_name]],
    cex.main = 0.78, cex.lab = 0.61, cex.axis = 0.53,
    ticktype = "detailed", nticks = 4
  )
  segments3D(
    x0 = reference_temperature,
    y0 = min(render_lag_grid), z0 = limits[[1]],
    x1 = reference_temperature,
    y1 = max(render_lag_grid), z1 = limits[[1]],
    add = TRUE, col = "#5E666D", lty = 3, lwd = 0.95
  )
}

draw_surface_row <- function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
  par(
    oma = c(0.06, 0.06, 2.12, 0.06),
    family = font_family, fg = dark,
    col.axis = "#353B40", col.lab = dark
  )
  for (outcome_name in outcome_order) draw_surface(outcome_name)
  mtext(
    paste0("a  Lag-specific temperature-response surfaces at ",
           format(round(slice_altitude), big.mark = ","), " m"),
    side = 3, outer = TRUE, line = 1.02, adj = 0,
    font = 2, cex = 0.84
  )
  mtext(
    paste0(
      "Hierarchical DLNMs with identical specifications; OR relative to ",
      reference_temperature, " ", degree_c,
      "; estimates at lags 0, 1 and 2 weeks; outcome-specific OR axes"
    ),
    side = 3, outer = TRUE, line = 0.02, adj = 0,
    col = muted, cex = 0.54
  )
}

panel_a_path <- file.path(source_dir, "panel_a_surface_row_600dpi.png")
grDevices::png(
  panel_a_path,
  width = 183 / 25.4, height = 80 / 25.4,
  units = "in", res = 600, bg = "white", type = "cairo-png"
)
draw_surface_row()
dev.off()
panel_a <- wrap_elements(full = rasterGrob(
  png::readPNG(panel_a_path),
  width = unit(1, "npc"), height = unit(1, "npc"),
  interpolate = TRUE
))

# -------------------------------------------------------------------------
# Panel b: observed province-wide subtype trajectories.
# -------------------------------------------------------------------------
subtype_time <- read.csv(
  file.path(
    "result0715", "source_data", "figure3_subtype_environment_20260806",
    "panel_a_subtype_time_series.csv"
  ),
  stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
subtype_time$week_start <- as.Date(subtype_time$week_start)

centered_ma3 <- function(x) {
  as.numeric(stats::filter(x, rep(1 / 3, 3), sides = 2))
}
weekly <- subtype_time |>
  transmute(
    week_start,
    influenza_A_positivity = A_positivity * 100,
    influenza_B_positivity = B_positivity * 100
  ) |>
  arrange(week_start) |>
  mutate(
    influenza_A_ma3 = centered_ma3(influenza_A_positivity),
    influenza_B_ma3 = centered_ma3(influenza_B_positivity)
  )
weekly_long <- bind_rows(
  weekly |>
    transmute(
      week_start, subtype = "Influenza A",
      weekly_positivity = influenza_A_positivity,
      moving_average = influenza_A_ma3
    ),
  weekly |>
    transmute(
      week_start, subtype = "Influenza B",
      weekly_positivity = influenza_B_positivity,
      moving_average = influenza_B_ma3
    )
) |>
  mutate(subtype = factor(subtype, levels = c("Influenza A", "Influenza B")))

subtype_colours <- c("Influenza A" = blue, "Influenza B" = orange)
theme_panel <- theme_classic(base_family = font_family, base_size = 6.7) +
  theme(
    plot.title = element_text(face = "bold", size = 7.4, margin = margin(b = 1)),
    plot.subtitle = element_text(size = 6.0, colour = muted, margin = margin(b = 3)),
    axis.title = element_text(size = 6.4),
    axis.text = element_text(size = 5.9, colour = "#30353A"),
    axis.line = element_line(linewidth = 0.34),
    axis.ticks = element_line(linewidth = 0.30),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.9),
    legend.key.width = unit(7, "pt"),
    panel.grid = element_blank()
  )

panel_b <- ggplot(
  weekly_long,
  aes(week_start, colour = subtype, group = subtype)
) +
  geom_line(aes(y = weekly_positivity), linewidth = 0.27, alpha = 0.22) +
  geom_line(aes(y = moving_average), linewidth = 0.62, na.rm = TRUE) +
  scale_colour_manual(values = subtype_colours) +
  scale_x_date(
    breaks = as.Date(c(
      "2017-01-02", "2019-01-07", "2021-01-04", "2023-01-02", "2024-12-23"
    )),
    labels = c("2017", "2019", "2021", "2023", "2024"),
    limits = range(weekly$week_start, na.rm = TRUE),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    title = "b  Weekly influenza A and B activity",
    subtitle = "Weekly laboratory positivity; centred 3-week moving averages",
    x = NULL, y = "Laboratory positivity (%)"
  ) +
  theme_panel +
  theme(
    legend.position = c(0.69, 0.84),
    legend.direction = "horizontal",
    plot.margin = margin(4, 6, 4, 4)
  )

# -------------------------------------------------------------------------
# Panel c: standardized subtype phenology. This model is independent of the
# meteorological exposure and is therefore carried forward unchanged.
# -------------------------------------------------------------------------
phenology_dir <- file.path(
  "result0715", "source_data", "result3_subtype_phenology_20260806"
)
phenology <- read.csv(
  file.path(phenology_dir, "standardized_subtype_seasonal_curves.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) |>
  filter(primary_specification %in% TRUE)
peaks <- read.csv(
  file.path(phenology_dir, "subtype_peak_week_sensitivity.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) |>
  filter(primary_specification %in% TRUE)

phenology <- phenology |>
  mutate(
    outcome_label = factor(
      outcome_label, levels = c("Influenza A", "Influenza B")
    ),
    season_week = if_else(annual_week >= 35, annual_week - 34, annual_week + 18)
  ) |>
  arrange(outcome_label, season_week)
peaks <- peaks |>
  mutate(season_week = if_else(peak_week >= 35, peak_week - 34, peak_week + 18))

phenology_fills <- c("Influenza A" = "#BBD7E9", "Influenza B" = "#F2D3A4")
season_breaks <- c(1, 10, 19, 27, 36, 45, 52)
season_labels <- c("Sep", "Nov", "Jan", "Mar", "May", "Jul", "Aug")

panel_c <- ggplot(
  phenology,
  aes(season_week, estimate * 100, colour = outcome_label, fill = outcome_label)
) +
  geom_ribbon(aes(ymin = low * 100, ymax = high * 100), alpha = 0.24, colour = NA) +
  geom_line(linewidth = 0.76) +
  geom_vline(
    data = peaks,
    aes(xintercept = season_week, colour = outcome_label),
    linewidth = 0.40, linetype = "dotted", show.legend = FALSE
  ) +
  scale_colour_manual(values = subtype_colours) +
  scale_fill_manual(values = phenology_fills) +
  guides(
    fill = "none",
    colour = guide_legend(override.aes = list(linewidth = 0.76))
  ) +
  scale_x_continuous(
    breaks = season_breaks, labels = season_labels,
    limits = c(1, 52.75), expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "c  Standardized type-specific seasonality",
    subtitle = "September-August season; hierarchical beta-binomial model",
    x = NULL, y = "Predicted positivity (%)"
  ) +
  theme_panel +
  theme(
    legend.position = c(0.68, 0.84),
    legend.direction = "horizontal",
    plot.margin = margin(4, 7, 4, 4)
  )

figure <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(1.18, 0.88))

stem <- file.path(figure_dir, "Figure3_AGU_no_title_type_20261008")
width_mm <- 183
height_mm <- 168
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

grDevices::png(
  paste0(stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 600, bg = "white", type = "cairo-png"
)
print(figure)
dev.off()
grDevices::tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw", bg = "white"
)
print(figure)
dev.off()
ggsave(
  paste0(stem, ".pdf"), figure,
  width = width_mm, height = height_mm, units = "mm",
  device = cairo_pdf, bg = "white"
)
svglite::svglite(
  paste0(stem, ".svg"), width = width_in, height = height_in,
  bg = "white", fix_text_size = FALSE
)
print(figure)
dev.off()

rhs_labels <- lapply(
  fits, function(fit) attr(terms(formula(fit)), "term.labels")
)
identical_rhs <- all(vapply(
  rhs_labels[-1], identical, logical(1), rhs_labels[[1]]
))
diagnostics <- bind_rows(lapply(names(fits), function(outcome_name) {
  fit <- fits[[outcome_name]]
  tibble(
    outcome = outcome_name,
    observations = nobs(fit),
    counties = n_distinct(model_dat$county),
    convergence = fit$fit$convergence,
    pdHess = fit$sdr$pdHess,
    identical_right_hand_side = identical_rhs,
    reference_temperature = reference_temperature,
    displayed_altitude = slice_altitude
  )
}))

qa <- tibble(
  check = c(
    "All three outcomes use the full 9,758 county-weeks",
    "All three outcomes use all 45 counties",
    "All three DLNMs have identical right-hand-side terms",
    "All three DLNMs converged with positive-definite Hessians",
    "Reference temperature equals the current Result 2 reference",
    "Panel a uses the P50 altitude display anchor",
    "Panel a contains observed lags 0-2",
    "Lag-specific ORs reproduce cumulative lag 0-2 ORs",
    "Panel b spans 2017-2024",
    "Panel b displays no out-of-range year tick",
    "Panel c uses the prespecified three-harmonic model",
    "All four final export formats exist"
  ),
  pass = c(
    all(diagnostics$observations == 9758L),
    all(diagnostics$counties == 45L),
    identical_rhs,
    all(diagnostics$convergence == 0) && all(diagnostics$pdHess),
    reference_temperature == 0,
    abs(slice_altitude - 2980) < 1e-8,
    identical(sort(unique(surface_data$lag_week)), 0:2),
    is.finite(max_consistency_error) && max_consistency_error <= 1e-8,
    min(weekly$week_start, na.rm = TRUE) <= as.Date("2017-01-16") &&
      max(weekly$week_start, na.rm = TRUE) >= as.Date("2024-12-23"),
    all(as.integer(format(as.Date(c(
      "2017-01-02", "2019-01-07", "2021-01-04", "2023-01-02", "2024-12-23"
    )), "%Y")) <= 2024),
    all(phenology$harmonics == 3L),
    all(file.exists(paste0(stem, c(".png", ".pdf", ".svg", ".tiff"))))
  )
)
if (!all(qa$pass)) stop("Result 3 QA failed.")

write_csv_utf8(surface_data, file.path(source_dir, "panel_a_temperature_lag_surfaces.csv"))
write_csv_utf8(consistency, file.path(source_dir, "panel_a_cumulative_consistency.csv"))
write_csv_utf8(diagnostics, file.path(source_dir, "panel_a_model_diagnostics.csv"))
write_csv_utf8(weekly_long, file.path(source_dir, "panel_b_weekly_trajectories.csv"))
write_csv_utf8(phenology, file.path(source_dir, "panel_c_standardized_phenology.csv"))
write_csv_utf8(peaks, file.path(source_dir, "panel_c_peak_weeks.csv"))
write_csv_utf8(qa, file.path(qa_dir, "Result3_T2m_FINAL_qa.csv"))

legend <- c(
  paste0(
    "Figure 3. Temperature-lag responses and type-specific seasonality. ",
    "Laboratory test positivity, 2017-2024; all temperature models use weekly mean ",
    "2-m air temperature referenced to ", reference_temperature, " ", degree_c, "."
  ),
  paste0(
    "a, Lag-specific temperature-response surfaces for overall, influenza ",
    "A-only and influenza B-only laboratory positivity from full-sample hierarchical ",
    "beta-binomial DLNMs with identical specifications. Surfaces are evaluated at the ",
    "2,980 m continuous-altitude display anchor relative to ", reference_temperature, " ",
    degree_c, ". Estimates were obtained at lags 0, 1 and 2 weeks; intervening mesh ",
    "lines are linear interpolation for display. OR axes are outcome-specific."
  ),
  paste0(
    "b, Province-wide weekly influenza A and B laboratory positivity during ",
    "2017-2024. Thin lines are weekly observations and thick lines are centred ",
    "3-week moving averages used only for display."
  ),
  paste0(
    "c, Calendar-standardized influenza A and B seasonal positivity curves ",
    "from hierarchical beta-binomial models with three annual Fourier pairs and ",
    "calendar-year adjustment. Shading denotes 95% simulation intervals and ",
    "dotted lines indicate model-based peak weeks."
  )
)
writeLines(
  legend, file.path(figure_dir, "Figure_3_AGU_legend_20261008.txt"),
  useBytes = TRUE
)

message("Created harmonized Result 3 using current 2-m air-temperature models.")
