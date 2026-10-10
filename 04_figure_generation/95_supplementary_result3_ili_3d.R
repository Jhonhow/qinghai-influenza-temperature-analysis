suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(glmmTMB)
  library(dlnm)
  library(plot3D)
})

set.seed(20260902)

# Supplementary Figure S2
# Temperature-lag response surfaces for the full ILI surveillance sample at
# supported P10/P50/P90 continuous-altitude display anchors.

degree_c <- paste0(intToUtf8(176), "C")
analysis_id <- "result3_ili_temperature_lag_3d_final_20260909"
model_root <- file.path("model_outputs", "reanalysis_v4_t2m_20260825")
panel_path <- file.path(model_root, "county_week_panel_v4_t2m_2016_2025.rds")
bundle_path <- file.path(
  model_root, "primary_result2_result4_models_20260825",
  "v4_t2m_model_bundle.rds"
)
fit_path <- file.path(
  model_root, "result5_robustness_20260901", "ili_hierarchical.rds"
)
figure_dir <- file.path(
  "result0715", "supplementary_figures", "v4_t2m_20260825", analysis_id
)
source_dir <- file.path(
  "result0715", "source_data", paste0("v4_", analysis_id)
)
qa_dir <- file.path(
  "result0715", "figure_qa", "v4_t2m_20260825", analysis_id
)
for (path in c(figure_dir, source_dir, qa_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

panel <- readRDS(panel_path)
bundle <- readRDS(bundle_path)
fit_ili <- readRDS(fit_path)
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
z_int_names <- bundle$z_int_names
interaction_scale <- bundle$interaction_scale
altitude_mean <- bundle$altitude_mean
reference_temperature <- 0
lag_grid <- 0:2

# Reconstruct the exact fitted ILI rows in the saved cross-basis order.
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
    iso_year = isoyear(week_start),
    iso_week = isoweek(week_start),
    week_id = sprintf("%04d-W%02d", iso_year, iso_week)
  ) |>
  left_join(
    panel |> filter(week_start >= week_min, week_start <= week_max),
    by = c("county", "week_start", "iso_year", "iso_week", "week_id"),
    suffix = c("", "_panel")
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
    week_index = as.integer((week_start - week_min) / 7) + 1L
  ) |>
  select(
    -ends_with("_panel"), -city_attribute, -altitude_attribute,
    -altitude_1000_c_attribute
  ) |>
  arrange(county, week_start)
if (nrow(weekly_grid) != nrow(cb_all)) stop("Saved cross-basis does not align.")

cb_df <- as.data.frame(cb_all)
names(cb_df) <- cb_names
ili_dat <- bind_cols(
  weekly_grid |> mutate(exposure_row_id = row_number()), cb_df
) |>
  filter(
    ili_available %in% TRUE,
    is.finite(ILI), is.finite(visits), visits >= 1,
    ILI >= 0, ILI <= visits,
    complete.cases(across(all_of(cb_names))),
    !is.na(altitude_1000_c), !is.na(city)
  ) |>
  arrange(county, week_index)
if (nrow(ili_dat) != nobs(fit_ili)) stop("ILI model rows were not reproduced.")

county_altitudes <- ili_dat |> distinct(county, altitude) |> arrange(altitude)
if (nrow(county_altitudes) != 18L) stop("Expected 18 ILI counties.")
anchor_values <- as.numeric(
  quantile(county_altitudes$altitude, c(.10, .50, .90), names = FALSE)
)
anchors <- tibble(
  altitude_anchor = factor(c("P10", "P50", "P90"), levels = c("P10", "P50", "P90")),
  altitude = anchor_values,
  altitude_display = floor(anchor_values + .5),
  altitude_1000_c = (anchor_values - altitude_mean) / 1000
)

local_support <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  local <- ili_dat |> filter(abs(altitude - anchors$altitude[[i]]) <= 250)
  tibble(
    altitude_anchor = anchors$altitude_anchor[[i]],
    altitude = anchors$altitude[[i]],
    county_weeks = nrow(local),
    counties = n_distinct(local$county),
    p05 = as.numeric(quantile(local$temperature, .05)),
    p95 = as.numeric(quantile(local$temperature, .95))
  )
})) |>
  mutate(
    display_lower = ceiling(p05 * 2) / 2,
    display_upper = floor(p95 * 2) / 2
  )
common_lower <- ceiling(max(local_support$p05) * 2) / 2
common_upper <- floor(min(local_support$p95) * 2) / 2
if (!(common_lower < 0 && common_upper > 0)) stop("Invalid common support.")
if (any(local_support$display_lower >= reference_temperature) ||
    any(local_support$display_upper <= reference_temperature)) {
  stop("A panel-specific support range does not contain the reference temperature.")
}

temperature_grid_for <- function(anchor_name) {
  support <- local_support |>
    filter(as.character(altitude_anchor) == anchor_name)
  sort(unique(c(
    seq(support$display_lower[[1]], support$display_upper[[1]], by = .25),
    reference_temperature
  )))
}

basis_array_lag <- function(temperatures, centre) {
  original_names <- colnames(cb_all)
  output <- array(
    0,
    dim = c(length(temperatures), length(lag_grid), length(original_names)),
    dimnames = list(as.character(temperatures), paste0("lag", lag_grid), cb_names)
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
    output[, , j] <- prediction$matfit[, paste0("lag", lag_grid), drop = FALSE]
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
      model.link = "logit", at = temperatures, cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  output
}

extract_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, variance = as.matrix(vcov(fit)$cond))
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
    rowSums((design %*% fixed$variance) * design), 0
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

surface_rows <- list()
for (i in seq_len(nrow(anchors))) {
  anchor_name <- as.character(anchors$altitude_anchor[[i]])
  temperature_grid <- temperature_grid_for(anchor_name)
  lag_basis <- basis_array_lag(temperature_grid, reference_temperature)
  support <- local_support |>
    filter(as.character(altitude_anchor) == anchor_name)
  for (temperature_index in seq_along(temperature_grid)) {
    for (lag_index in seq_along(lag_grid)) {
      basis_row <- matrix(
        lag_basis[temperature_index, lag_index, ],
        nrow = 1, dimnames = list(NULL, cb_names)
      )
      estimate <- estimate_basis_rows(
        fit_ili, basis_row, anchors$altitude_1000_c[[i]]
      )
      surface_rows[[length(surface_rows) + 1L]] <- bind_cols(
        tibble(
          altitude_anchor = anchors$altitude_anchor[[i]],
          altitude = anchors$altitude[[i]],
          altitude_display = anchors$altitude_display[[i]],
          temperature = temperature_grid[[temperature_index]],
          lag_week = lag_grid[[lag_index]],
          reference_temperature = reference_temperature,
          common_support_lower = common_lower,
          common_support_upper = common_upper,
          local_support_lower = support$display_lower[[1]],
          local_support_upper = support$display_upper[[1]]
        ),
        estimate
      )
    }
  }
}
surface_data <- bind_rows(surface_rows)

# QA: lag-specific products reproduce cumulative estimates from the same model.
cumulative_expected <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  anchor_name <- as.character(anchors$altitude_anchor[[i]])
  temperature_grid <- temperature_grid_for(anchor_name)
  cumulative_basis <- basis_matrix_cumulative(
    temperature_grid, reference_temperature
  )
  bind_cols(
    tibble(
      altitude_anchor = anchors$altitude_anchor[[i]],
      temperature = temperature_grid
    ),
    estimate_basis_rows(
      fit_ili, cumulative_basis, anchors$altitude_1000_c[[i]]
    ) |> select(expected_cumulative_or = estimate)
  )
}))
cumulative_reconstructed <- surface_data |>
  group_by(altitude_anchor, temperature) |>
  summarise(
    reconstructed_cumulative_or = exp(sum(log_estimate)),
    .groups = "drop"
  )
consistency <- cumulative_reconstructed |>
  left_join(cumulative_expected, by = c("altitude_anchor", "temperature")) |>
  mutate(
    absolute_log_difference = abs(
      log(reconstructed_cumulative_or) - log(expected_cumulative_or)
    )
  )
max_consistency_error <- max(consistency$absolute_log_difference)

fixed <- extract_fixed(fit_ili)
beta_int <- fixed$beta[z_int_names]
variance_int <- fixed$variance[z_int_names, z_int_names, drop = FALSE]
wald_chisq <- as.numeric(t(beta_int) %*% solve(variance_int, beta_int))
wald_df <- length(beta_int)
wald_p <- pchisq(wald_chisq, df = wald_df, lower.tail = FALSE)

# Render each surface within its altitude-specific local P5-P95 temperature
# support. The z axis remains shared; x ranges are deliberately panel-specific
# to avoid truncating supported cold exposure at the P90 anchor.
render_temperature_grid_for <- function(anchor_name) {
  temperature_grid <- temperature_grid_for(anchor_name)
  output <- temperature_grid[seq(1, length(temperature_grid), by = 2)]
  if (tail(output, 1) != tail(temperature_grid, 1)) {
    output <- c(output, tail(temperature_grid, 1))
  }
  output
}
render_lag_grid <- seq(0, 2, by = .1)
surface_matrix <- function(anchor_name) {
  render_temperature_grid <- render_temperature_grid_for(anchor_name)
  dat <- surface_data |>
    filter(
      altitude_anchor == anchor_name,
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
    function(i) approx(
      lag_grid, observed[i, ], xout = render_lag_grid,
      method = "linear", rule = 2
    )$y,
    numeric(length(render_lag_grid))
  ))
  dimnames(output) <- list(render_temperature_grid, render_lag_grid)
  output
}
all_rendered <- unlist(lapply(levels(anchors$altitude_anchor), surface_matrix))
z_limits <- range(c(1, all_rendered), finite = TRUE)
z_span <- diff(z_limits)
z_limits <- c(
  max(.05, z_limits[[1]] - .08 * z_span),
  z_limits[[2]] + .14 * z_span
)
surface_colours <- colorRampPalette(
  c("#FFF9F0", "#F7E5C6", "#EBC98D", "#D79B45", "#B96C18")
)(100)

draw_surface <- function(i) {
  anchor_name <- as.character(anchors$altitude_anchor[[i]])
  render_temperature_grid <- render_temperature_grid_for(anchor_name)
  z <- surface_matrix(anchor_name)
  # Symmetric inner margins keep the perspective box inside each layout cell.
  # The additional right margin is especially important for the P90 panel,
  # whose rear upper edge otherwise reaches the device boundary when exported.
  par(mar = c(2.18, 1.18, 1.42, 1.18), xpd = FALSE)
  persp3D(
    x = render_temperature_grid,
    y = render_lag_grid,
    z = z,
    colvar = z,
    col = surface_colours,
    breaks = seq(z_limits[[1]], z_limits[[2]], length.out = 101),
    clim = z_limits, zlim = z_limits,
    theta = 31, phi = 24, expand = .75, scale = TRUE,
    facets = TRUE,
    border = grDevices::adjustcolor("#A36A2D", alpha.f = .68),
    lwd = .30, lighting = FALSE, shade = NA,
    colkey = FALSE, bty = "b",
    xlab = paste0("2-m air temperature (", degree_c, ")"),
    ylab = "Lag (weeks)", zlab = "OR",
    main = sprintf(
      "%s: %s m", anchor_name,
      format(anchors$altitude_display[[i]], big.mark = ",")
    ),
    cex.main = .80, cex.lab = .60, cex.axis = .52,
    ticktype = "detailed", nticks = 4
  )
  segments3D(
    x0 = reference_temperature,
    y0 = min(render_lag_grid), z0 = z_limits[[1]],
    x1 = reference_temperature,
    y1 = max(render_lag_grid), z1 = z_limits[[1]],
    add = TRUE, col = "#5E666D", lty = 3, lwd = .95
  )
}

stem <- file.path(figure_dir, "Supplementary_Figure_S2_ILI_temperature_lag_surfaces")
width_mm <- 183
height_mm <- 88
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

draw_figure <- function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
  par(
    oma = c(.10, .12, 2.42, .22), family = "Arial",
    fg = "#202428", col.axis = "#353B40", col.lab = "#202428"
  )
  for (i in seq_len(nrow(anchors))) draw_surface(i)
  mtext(
    "ILI temperature-lag response surfaces at supported altitudes",
    side = 3, outer = TRUE, line = 1.25, adj = 0,
    font = 2, cex = .92, col = "#172A3A"
  )
  mtext(
    paste0(
      "Full 2017-2024 ILI sample; panel-specific local P5-P95 support; ",
      "reference 0 ", degree_c,
      "; lags 0-2 weeks; shared OR axis; global interaction P = ",
      sprintf("%.3f", wald_p)
    ),
    side = 3, outer = TRUE, line = .20, adj = 0,
    cex = .57, col = "#687A8B"
  )
}

grDevices::png(
  paste0(stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 600, bg = "white", type = "cairo-png"
)
draw_figure()
dev.off()
grDevices::tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw", bg = "white"
)
draw_figure()
dev.off()
grDevices::cairo_pdf(
  paste0(stem, ".pdf"), width = width_in, height = height_in,
  family = "Arial"
)
draw_figure()
dev.off()
grDevices::svg(
  paste0(stem, ".svg"), width = width_in, height = height_in,
  family = "Arial", bg = "white"
)
draw_figure()
dev.off()

write.csv(
  surface_data, file.path(source_dir, "supplementary_figure_s2_ili_lag_surfaces.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  local_support, file.path(source_dir, "supplementary_figure_s2_local_support.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  consistency, file.path(source_dir, "supplementary_figure_s2_cumulative_consistency.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

legend <- paste0(
  "Supplementary Figure S2. ILI temperature-lag response surfaces at supported altitudes. ",
  "Lag-specific odds ratios for the ILI consultation rate are shown at the P10 (2,080 m), P50 (2,563 m), and P90 (3,423 m) altitude anchors. Weekly mean 2-m air temperature is referenced to 0 ", degree_c, ". Each panel is restricted to its local P5-P95 temperature support (P10: ", sprintf("%.1f to %.1f", local_support$display_lower[[1]], local_support$display_upper[[1]]), " ", degree_c, "; P50: ", sprintf("%.1f to %.1f", local_support$display_lower[[2]], local_support$display_upper[[2]]), " ", degree_c, "; P90: ", sprintf("%.1f to %.1f", local_support$display_lower[[3]], local_support$display_upper[[3]]), " ", degree_c, "), while all panels use the same OR axis. Estimates at observed lags 0, 1, and 2 weeks are joined by linear interpolation for display only. The hierarchical beta-binomial distributed-lag nonlinear model included seasonal harmonics, long-term trend, pandemic period, continuous altitude, and temperature-by-altitude interaction, with prefecture and county random intercepts. The analysis included ", format(nrow(ili_dat), big.mark = ","), " county-weeks from 18 counties; the global interaction test gave P = ", sprintf("%.3f", wald_p), ". ILI is a syndromic outcome and is not equivalent to laboratory-confirmed influenza positivity."
)
writeLines(
  legend,
  file.path(figure_dir, "Supplementary_Figure_S2_legend.md"),
  useBytes = TRUE
)

qa <- tibble(
  check = c(
    "Full fitted ILI sample reconstructed",
    "Expected 18-county coverage",
    "ILI model converged",
    "Positive-definite Hessian",
    "All anchors lie inside observed altitude range",
    "Each panel-specific support contains the reference temperature",
    "Lag-specific ORs reproduce cumulative ORs",
    "Shared OR axis contains all rendered estimates",
    "All four export formats exist"
  ),
  pass = c(
    nrow(ili_dat) == nobs(fit_ili),
    n_distinct(ili_dat$county) == 18L,
    fit_ili$fit$convergence == 0,
    isTRUE(fit_ili$sdr$pdHess),
    all(anchor_values >= min(county_altitudes$altitude) &
          anchor_values <= max(county_altitudes$altitude)),
    all(local_support$display_lower < reference_temperature &
          local_support$display_upper > reference_temperature),
    max_consistency_error <= 1e-8,
    min(all_rendered) >= z_limits[[1]] && max(all_rendered) <= z_limits[[2]],
    all(file.exists(paste0(stem, c(".png", ".pdf", ".svg", ".tiff"))))
  ),
  detail = c(
    sprintf("n=%s county-weeks", format(nrow(ili_dat), big.mark = ",")),
    sprintf("counties=%d", n_distinct(ili_dat$county)),
    sprintf("convergence=%d", fit_ili$fit$convergence),
    sprintf("pdHess=%s", fit_ili$sdr$pdHess),
    sprintf("observed altitude %.0f to %.0f m", min(county_altitudes$altitude), max(county_altitudes$altitude)),
    paste(sprintf(
      "%s %.1f to %.1f C",
      as.character(local_support$altitude_anchor),
      local_support$display_lower,
      local_support$display_upper
    ), collapse = "; "),
    sprintf("maximum log error %.3g", max_consistency_error),
    sprintf("shared OR axis %.3f to %.3f", z_limits[[1]], z_limits[[2]]),
    "PNG, PDF, SVG, and TIFF"
  )
)
write.csv(
  qa, file.path(qa_dir, "Supplementary_Figure_S2_QA.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
if (any(!qa$pass)) stop("Supplementary Figure S2 failed QA.")

message("Supplementary Figure S2 completed: ", stem)
