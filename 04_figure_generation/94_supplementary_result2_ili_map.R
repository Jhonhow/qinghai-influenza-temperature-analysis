suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(glmmTMB)
  library(dlnm)
  library(sf)
  library(svglite)
  library(scales)
})

set.seed(20260902)

# Supplementary Figure S1 final
# a: geographic coverage of the 18-county ILI network.
# b-d: ILI temperature-response curves at supported P10/P50/P90 altitudes.
# The full ILI sample is used because laboratory curves are not displayed in
# this figure; no ILI observations are discarded merely to match laboratory
# availability.

analysis_id <- "result2_ili_s1_final_20260909"
model_root <- file.path("model_outputs", "reanalysis_v4_t2m_20260825")
panel_path <- file.path(model_root, "county_week_panel_v4_t2m_2016_2025.rds")
bundle_path <- file.path(
  model_root, "primary_result2_result4_models_20260825",
  "v4_t2m_model_bundle.rds"
)
ili_fit_path <- file.path(
  model_root, "result5_robustness_20260901", "ili_hierarchical.rds"
)
figure_dir <- file.path(
  "result0715", "supplementary_figures", "v4_t2m_20260825", analysis_id
)
source_dir <- file.path(
  "result0715", "source_data", paste0("v4_", analysis_id)
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

degree_c <- paste0(intToUtf8(176), "C")
reference_temperature <- 0

panel <- readRDS(panel_path)
bundle <- readRDS(bundle_path)
fit_ili <- readRDS(ili_fit_path)
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
z_int_names <- bundle$z_int_names
interaction_scale <- bundle$interaction_scale
altitude_mean <- bundle$altitude_mean

# Reconstruct the complete ILI model rows in the same order as the saved
# cross-basis. This reproduces the data construction in the fitted ILI model.
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

if (nrow(weekly_grid) != nrow(cb_all)) {
  stop("Weekly grid and saved temperature cross-basis do not align.")
}
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

if (nrow(ili_dat) != nobs(fit_ili)) {
  stop(sprintf(
    "Reconstructed ILI rows (%d) differ from fitted observations (%d).",
    nrow(ili_dat), nobs(fit_ili)
  ))
}

county_altitudes <- ili_dat |>
  distinct(county, altitude) |>
  arrange(altitude)
if (nrow(county_altitudes) != 18L) stop("Expected 18 ILI-covered counties.")

anchor_values <- as.numeric(
  quantile(county_altitudes$altitude, c(.10, .50, .90), names = FALSE)
)
anchor_display <- floor(anchor_values + .5)
altitude_anchors <- tibble(
  altitude_anchor = factor(c("P10", "P50", "P90"),
                           levels = c("P10", "P50", "P90")),
  altitude = anchor_values,
  altitude_display = anchor_display,
  altitude_1000_c = (anchor_values - altitude_mean) / 1000,
  altitude_label = factor(
    sprintf(
      "%s: %s m", c("P10", "P50", "P90"),
      format(anchor_display, big.mark = ",", scientific = FALSE)
    ),
    levels = sprintf(
      "%s: %s m", c("P10", "P50", "P90"),
      format(anchor_display, big.mark = ",", scientific = FALSE)
    )
  )
)

local_support <- bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
  anchor <- altitude_anchors$altitude[[i]]
  z <- ili_dat |> filter(abs(altitude - anchor) <= 250)
  if (nrow(z) < 100) stop("Insufficient local temperature support.")
  tibble(
    altitude_anchor = altitude_anchors$altitude_anchor[[i]],
    altitude = anchor,
    local_county_weeks = nrow(z),
    local_counties = n_distinct(z$county),
    p05 = as.numeric(quantile(z$temperature, .05)),
    p10 = as.numeric(quantile(z$temperature, .10)),
    median = median(z$temperature),
    p90 = as.numeric(quantile(z$temperature, .90)),
    p95 = as.numeric(quantile(z$temperature, .95))
  )
}))

common_lower <- ceiling(max(local_support$p05) * 2) / 2
common_upper <- floor(min(local_support$p95) * 2) / 2
if (!(common_lower < -5 && common_upper > 5 &&
      reference_temperature >= common_lower &&
      reference_temperature <= common_upper)) {
  stop("The derived common support is unsuitable for the planned contrasts.")
}

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

extract_fixed <- function(fit) {
  list(beta = fixef(fit)$cond, variance = as.matrix(vcov(fit)$cond))
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

temperature_grid <- sort(unique(c(
  seq(common_lower, common_upper, by = .25),
  reference_temperature, -5
)))
basis <- basis_matrix_cumulative(temperature_grid, reference_temperature)

curves <- bind_rows(lapply(seq_len(nrow(altitude_anchors)), function(i) {
  anchor <- altitude_anchors[i, ]
  bind_cols(
    tibble(
      altitude_anchor = anchor$altitude_anchor,
      altitude = anchor$altitude,
      altitude_display = anchor$altitude_display,
      altitude_label = anchor$altitude_label,
      temperature = as.numeric(rownames(basis)),
      reference_temperature = reference_temperature,
      common_support_lower = common_lower,
      common_support_upper = common_upper
    ),
    estimate_curve(fit_ili, basis, anchor$altitude_1000_c[[1]])
  )
})) |>
  mutate(
    altitude_anchor = factor(altitude_anchor, levels = c("P10", "P50", "P90")),
    altitude_label = factor(
      altitude_label, levels = levels(altitude_anchors$altitude_label)
    )
  )

key_contrasts <- curves |>
  filter(abs(temperature + 5) < 1e-12) |>
  mutate(
    label = sprintf(
      "-5 vs 0 %s\nOR %.2f (%.2f-%.2f)",
      degree_c, estimate, low, high
    )
  )

fixed <- extract_fixed(fit_ili)
beta_int <- fixed$beta[z_int_names]
variance_int <- fixed$variance[z_int_names, z_int_names, drop = FALSE]
wald_chisq <- as.numeric(t(beta_int) %*% solve(variance_int, beta_int))
wald_df <- length(beta_int)
wald_p <- pchisq(wald_chisq, df = wald_df, lower.tail = FALSE)

# Administrative-boundary map --------------------------------------------
# GDAL cannot reliably open the non-ASCII project path on this Windows/R
# locale, so copy the shapefile components into an ASCII temporary cache.
map_source_dir <- file.path(
  "Codex_transfer_20260828", "external_inputs", "geography",
  "qinghai_shp"
)
map_cache <- file.path(Sys.getenv("TEMP"), "qh_result2_supp_map")
dir.create(map_cache, recursive = TRUE, showWarnings = FALSE)
map_extensions <- c("shp", "shx", "dbf", "prj", "cpg")
for (extension in map_extensions) {
  source_file <- file.path(map_source_dir, paste0("china_county.", extension))
  cached_file <- file.path(map_cache, paste0("china_county.", extension))
  if (!file.exists(cached_file)) {
    if (!file.copy(source_file, cached_file, overwrite = TRUE)) {
      stop("Could not create the temporary administrative-boundary cache.")
    }
  }
}

county_map_all <- suppressWarnings(st_read(
  file.path(map_cache, "china_county.shp"), quiet = TRUE
))
# Positional field 6 is the province code in the supplied boundary dataset.
qinghai_map <- county_map_all[
  as.character(county_map_all[[6]]) == "630000",
] |>
  st_make_valid()

ili_map_names <- c(
  "Minhehuizuandtuzu", "Ledu", "Pingan", "Guide",
  "Chengzhong", "Chengxi", "Chengbei", "Chengdong",
  "Tongren", "Huangzhong", "Datonghuizutuzu", "Geermu",
  "Gonghe", "Delingha", "Haiyan", "Gangcha", "Yushu", "Maqin"
)
qinghai_map$coverage <- ifelse(
  qinghai_map$name %in% ili_map_names,
  "ILI surveillance county", "Other county"
)
if (sum(qinghai_map$coverage == "ILI surveillance county") != 18L ||
    nrow(qinghai_map) != 45L) {
  stop("The Qinghai boundary map did not resolve 18 of 45 ILI counties.")
}
province_outline <- st_union(st_geometry(qinghai_map))

# Plot --------------------------------------------------------------------
# Panel b uses the same visual grammar as main-text Result 2 panel a, but is
# restricted to the ILI analytical sample.
heatmap_data <- ili_dat |>
  mutate(
    temperature_bin = floor(temperature) + .5,
    altitude_bin = floor(altitude / 100) * 100 + 50
  ) |>
  count(temperature_bin, altitude_bin, name = "county_weeks")

theme_axes <- theme_classic(base_size = 6.7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = .35, colour = "#303030"),
    axis.ticks = element_line(linewidth = .3, colour = "#303030"),
    axis.text = element_text(size = 5.8, colour = "#222222"),
    axis.title = element_text(size = 6.4, colour = "#222222"),
    plot.title = element_text(size = 7.2, face = "bold", colour = "#1D2C38"),
    plot.subtitle = element_text(size = 5.8, colour = "#6B7B8C"),
    plot.margin = margin(4, 5, 4, 4)
  )

p_map <- ggplot() +
  geom_sf(
    data = qinghai_map,
    aes(fill = coverage), colour = "white", linewidth = .20
  ) +
  geom_sf(
    data = province_outline, fill = NA, colour = "#53616D", linewidth = .48
  ) +
  scale_fill_manual(
    values = c(
      "ILI surveillance county" = "#9B5D34",
      "Other county" = "#EADBC8"
    ),
    breaks = c("ILI surveillance county", "Other county")
  ) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(
    title = "ILI surveillance footprint",
    subtitle = "18 of 45 study counties",
    fill = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 7.2, face = "bold", colour = "#1D2C38"),
    plot.subtitle = element_text(size = 5.9, colour = "#6B7B8C", margin = margin(b = 2)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.text = element_text(size = 5.2),
    legend.key.height = unit(5, "pt"),
    legend.key.width = unit(8, "pt"),
    legend.spacing.x = unit(2, "pt"),
    plot.margin = margin(4, 7, 3, 3)
  )

p_heat <- ggplot(
  heatmap_data,
  aes(x = temperature_bin, y = altitude_bin, fill = county_weeks)
) +
  geom_tile(width = .96, height = 96, colour = NA) +
  scale_fill_gradientn(
    colours = c("#F5F8F9", "#D8E7EC", "#9FC5D1", "#4B8798"),
    trans = "sqrt",
    breaks = pretty_breaks(n = 4),
    guide = guide_colourbar(
      title.position = "top", title.hjust = 0,
      barwidth = unit(34, "mm"), barheight = unit(2.2, "mm")
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(.01, .01))) +
  scale_y_continuous(
    breaks = c(2000, 2500, 3000, 3500),
    expand = expansion(mult = c(.01, .01))
  ) +
  labs(
    title = "Observed 2-m air-temperature support in ILI surveillance counties",
    subtitle = "Binned county-week observations; darker tiles indicate more observations",
    x = paste0("Weekly mean 2-m air temperature (", degree_c, ")"),
    y = "County altitude (m)",
    fill = "County-weeks"
  ) +
  theme_axes +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.title = element_text(size = 5.3, face = "bold"),
    legend.text = element_text(size = 4.9),
    legend.margin = margin(-3, 0, -1, 0),
    plot.margin = margin(4, 4, 3, 5)
  )

p_curves <- ggplot(curves, aes(x = temperature, y = estimate)) +
  geom_hline(
    yintercept = 1, linewidth = .35, linetype = "dashed", colour = "#7B8995"
  ) +
  geom_vline(
    xintercept = 0, linewidth = .35, linetype = "dotted", colour = "#3F464C"
  ) +
  geom_ribbon(
    aes(ymin = low, ymax = high), fill = "#EBD4B2", alpha = .72, linewidth = 0
  ) +
  geom_line(colour = "#C67D1D", linewidth = .8) +
  geom_pointrange(
    data = key_contrasts,
    aes(x = temperature, y = estimate, ymin = low, ymax = high),
    inherit.aes = FALSE, colour = "#9F5E12", fill = "white",
    shape = 21, size = .25, linewidth = .42
  ) +
  geom_label(
    data = key_contrasts,
    aes(x = common_lower + .2, y = 1.47, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1,
    size = 2.0, lineheight = .95, linewidth = 0,
    label.padding = unit(1.3, "pt"), fill = alpha("white", .86),
    colour = "#5B4530"
  ) +
  facet_wrap(~ altitude_label, nrow = 1) +
  scale_x_continuous(
    limits = c(common_lower, common_upper),
    breaks = c(-5, 0, 5, 10)
  ) +
  scale_y_log10(
    limits = c(.54, 1.50),
    breaks = c(.6, .75, 1, 1.25, 1.5),
    labels = label_number(accuracy = .01)
  ) +
  labs(
    title = "ILI temperature-response curves at supported altitudes",
    subtitle = paste0(
      "Cumulative lag 0-2 ORs; reference: 0 ", degree_c,
      "; common P5-P95 temperature support: ",
      sprintf("%.1f to %.1f ", common_lower, common_upper), degree_c
    ),
    x = paste0("Weekly mean 2-m air temperature (", degree_c, ")"),
    y = "Cumulative odds ratio"
  ) +
  theme_axes +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6.6, face = "bold", colour = "#1D2C38"),
    panel.spacing.x = unit(9, "pt"),
    plot.margin = margin(4, 6, 4, 4)
  )

top_row <- p_map | p_heat
top_row <- top_row + plot_layout(widths = c(.95, 2.05))
figure_plot <- top_row / p_curves
figure_plot <- figure_plot +
  plot_layout(heights = c(.78, 1.18)) +
  plot_annotation(
    title = "ILI surveillance coverage, temperature support and response",
    subtitle = paste0(
      "Full 2017-2024 ILI sample (18 counties; ",
      format(nrow(ili_dat), big.mark = ","),
      " county-weeks); weekly mean 2-m air temperature; global interaction P = ",
      sprintf("%.3f", wald_p)
    ),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(
        family = "Arial", size = 10.3, face = "bold", colour = "#172A3A"
      ),
      plot.subtitle = element_text(
        family = "Arial", size = 6.4, colour = "#66798A", margin = margin(b = 4)
      ),
      plot.tag = element_text(
        family = "Arial", size = 8, face = "bold", colour = "#111111"
      ),
      plot.margin = margin(7, 8, 6, 7)
    )
  )

base_name <- file.path(
  figure_dir,
  "Supplementary_Figure_S1_ILI_coverage_support_response_FINAL"
)
width_mm <- 183
height_mm <- 132
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
  units = "in", res = 600, bg = "white", compression = "lzw", type = "cairo"
)
print(figure_plot)
dev.off()

grDevices::png(
  paste0(base_name, ".png"), width = width_in, height = height_in,
  units = "in", res = 300, bg = "white", type = "cairo"
)
print(figure_plot)
dev.off()

# Source data, legend, and QA ---------------------------------------------
write.csv(
  curves,
  file.path(source_dir, "supplementary_figure_s1_ili_curves.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  key_contrasts,
  file.path(source_dir, "supplementary_figure_s1_ili_minus5_vs_0.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  local_support,
  file.path(source_dir, "supplementary_figure_s1_ili_local_support.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  heatmap_data,
  file.path(source_dir, "supplementary_figure_s1_ili_temperature_altitude_support.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  st_drop_geometry(qinghai_map) |>
    select(county_code = 2, map_name = name, coverage),
  file.path(source_dir, "supplementary_figure_s1_ili_map_coverage.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

legend_lines <- c(
  "Supplementary Figure S1. ILI surveillance coverage, temperature support and response in Qinghai.",
  paste0(
    "a, Geographic coverage of the 18 counties contributing ILI surveillance data. ",
    "b, Distribution of county-week observations by weekly mean 2-m air temperature and county altitude; darker tiles indicate more observations. ",
    "c, Cumulative lag 0-2 associations between temperature and ILI consultation rate at the P10 (2,080 m), P50 (2,563 m) and P90 (3,423 m) continuous-altitude display anchors, referenced to 0 ", degree_c, " and restricted to common P5-P95 temperature support. Shading indicates 95% confidence intervals; points show the -5 versus 0 ", degree_c, " comparison. The global temperature-by-altitude interaction was not significant (P=", sprintf("%.3f", wald_p), "). ILI is a syndromic rather than laboratory-confirmed outcome."
  )
)
writeLines(
  legend_lines,
  file.path(figure_dir, "Supplementary_Figure_S1_legend_FINAL.md"),
  useBytes = TRUE
)

qa <- tibble(
  check = c(
    "Full fitted ILI sample reconstructed",
    "Expected ILI county coverage",
    "ILI model converged",
    "Positive-definite Hessian",
    "Reference equals OR 1 at 0 C",
    "Common support contains -5, 0, and 5 C",
    "All displayed CIs lie inside the chosen y-axis",
    "All three altitude anchors are observed",
    "Map resolves 18 of 45 counties",
    "Heatmap counts reproduce fitted ILI sample",
    "No unsupported 4,000-m ILI slice"
  ),
  pass = c(
    nrow(ili_dat) == nobs(fit_ili),
    n_distinct(ili_dat$county) == 18,
    fit_ili$fit$convergence == 0,
    isTRUE(fit_ili$sdr$pdHess),
    max(abs(curves$estimate[curves$temperature == 0] - 1)) < 1e-8,
    common_lower < -5 && common_upper > 5,
    min(curves$low) >= .54 && max(curves$high) <= 1.50,
    all(anchor_values >= min(county_altitudes$altitude) &
        anchor_values <= max(county_altitudes$altitude)),
    nrow(qinghai_map) == 45 &&
      sum(qinghai_map$coverage == "ILI surveillance county") == 18,
    sum(heatmap_data$county_weeks) == nrow(ili_dat),
    max(anchor_values) < 4000
  ),
  detail = c(
    sprintf("n=%s county-weeks", format(nrow(ili_dat), big.mark = ",")),
    sprintf("counties=%d", n_distinct(ili_dat$county)),
    sprintf("convergence=%d", fit_ili$fit$convergence),
    sprintf("pdHess=%s", fit_ili$sdr$pdHess),
    "reference error < 1e-8",
    sprintf("common support %.1f to %.1f C", common_lower, common_upper),
    sprintf("CI range %.3f to %.3f", min(curves$low), max(curves$high)),
    sprintf("altitude range %.0f to %.0f m", min(county_altitudes$altitude), max(county_altitudes$altitude)),
    "18 highlighted; 27 other counties",
    sprintf("heatmap n=%s county-weeks", format(sum(heatmap_data$county_weeks), big.mark = ",")),
    sprintf("highest display anchor %.0f m", max(anchor_values))
  )
)
write.csv(
  qa, file.path(figure_dir, "Supplementary_Figure_S1_QA_FINAL.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
if (any(!qa$pass)) stop("Supplementary Figure S1 candidate failed QA.")

message("Supplementary Figure S1 final completed: ", base_name)
