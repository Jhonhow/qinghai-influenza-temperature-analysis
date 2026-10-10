suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(dlnm)
  library(glmmTMB)
})

# Figure contract
# Core conclusion: the absolute-temperature response is nonlinear, but the
# full-sample altitude modification is limited and imprecisely estimated.
# Archetype: quantitative grid; R-only; double-column manuscript figure.

degree_c <- "\u00B0C"
font_family <- "sans"
dark <- "#202428"
muted <- "#687A8B"
blue <- "#2F6FA9"
blue_mid <- "#7FA9CA"
blue_light <- "#DDE9F4"
tile_low <- "#E8F1F4"
tile_high <- "#4E839E"
orange <- "#C8872B"
orange_light <- "#F2DFC3"
grey <- "#A7B0B8"
null_grey <- "#8F9CA7"

model_dir <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825"
)
figure_dir <- file.path("outputs", "result2_revision_20260908")
source_dir <- file.path(figure_dir, "source_data")
qa_dir <- file.path(figure_dir, "qa")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(file.path(model_dir, "v4_t2m_model_bundle.rds"))
model_dat <- bundle$model_data
fit <- bundle$fits$positive_full
cb_all <- bundle$cb_all
cb_names <- bundle$cb_names
z_int_names <- bundle$z_int_names
interaction_scale <- bundle$interaction_scale
anchors <- bundle$altitude_anchors
support <- bundle$local_support
altitude_mean <- bundle$altitude_mean
reference_temperature <- 0

fixed_beta <- fixef(fit)$cond
fixed_vcov <- as.matrix(vcov(fit)$cond)
altitude_span_km <- diff(anchors$altitude_1000_c[c(1, 3)])

common_p05 <- max(support$p05)
common_p95 <- min(support$p95)
common_p10 <- max(support$p10)
common_p90 <- min(support$p90)

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

basis_cumulative <- function(temperatures, centre = reference_temperature) {
  temperatures <- sort(unique(c(temperatures, centre)))
  original_names <- colnames(cb_all)
  output <- matrix(
    NA_real_, length(temperatures), length(original_names),
    dimnames = list(as.character(temperatures), original_names)
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
      model.link = "logit", at = temperatures,
      cen = centre, cumul = TRUE
    )
    output[, j] <- prediction$cumfit[, "lag2"]
  }
  colnames(output) <- cb_names
  output
}

estimate_basis <- function(basis, altitude_1000_c) {
  design <- matrix(
    0, nrow = nrow(basis), ncol = length(fixed_beta),
    dimnames = list(NULL, names(fixed_beta))
  )
  design[, cb_names] <- basis
  design[, z_int_names] <-
    sweep(basis, 2, interaction_scale, "/") * altitude_1000_c
  log_estimate <- as.numeric(design %*% fixed_beta)
  standard_error <- sqrt(pmax(
    rowSums((design %*% fixed_vcov) * design), 0
  ))
  tibble(
    estimate = exp(log_estimate),
    low = exp(log_estimate - 1.96 * standard_error),
    high = exp(log_estimate + 1.96 * standard_error),
    log_estimate = log_estimate,
    standard_error = standard_error
  )
}

estimate_contrast_at_altitude <- function(colder, warmer, altitude) {
  basis <- basis_cumulative(c(colder, warmer))
  basis_difference <- basis[as.character(colder), ] -
    basis[as.character(warmer), ]
  repeated_basis <- matrix(
    basis_difference, nrow = length(altitude), ncol = length(basis_difference),
    byrow = TRUE, dimnames = list(NULL, names(basis_difference))
  )
  estimate_basis(
    repeated_basis,
    (altitude - altitude_mean) / 1000
  )
}

estimate_moving_modification <- function(width = 5) {
  colder <- seq(-15, 10, by = .25)
  warmer <- colder + width
  temperatures <- sort(unique(c(colder, warmer, reference_temperature)))
  basis <- basis_cumulative(temperatures)
  get_basis <- function(x) {
    basis[match(as.character(x), rownames(basis)), , drop = FALSE]
  }
  basis_difference <- get_basis(colder) - get_basis(warmer)
  design <- matrix(
    0, nrow = length(colder), ncol = length(fixed_beta),
    dimnames = list(NULL, names(fixed_beta))
  )
  design[, z_int_names] <-
    sweep(basis_difference, 2, interaction_scale, "/") * altitude_span_km
  log_raw <- as.numeric(design %*% fixed_beta)
  se_raw <- sqrt(pmax(rowSums((design %*% fixed_vcov) * design), 0))
  log_per_1000m <- log_raw / altitude_span_km
  se_per_1000m <- se_raw / altitude_span_km

  tibble(
    contrast_width = width,
    colder_temperature = colder,
    warmer_temperature = warmer,
    raw_p90_p10_ratio = exp(log_raw),
    estimate_per_1000m = exp(log_per_1000m),
    low_per_1000m = exp(log_per_1000m - 1.96 * se_per_1000m),
    high_per_1000m = exp(log_per_1000m + 1.96 * se_per_1000m),
    common_p05_p95 =
      colder >= common_p05 & warmer <= common_p95,
    common_p10_p90 =
      colder >= common_p10 & warmer <= common_p90
  ) |>
    mutate(
      support_class = case_when(
        common_p10_p90 ~ "P10-P90",
        common_p05_p95 ~ "P05-P95",
        TRUE ~ "Extrapolation"
      )
    )
}

# Panel a source data
density_dat <- model_dat |>
  transmute(
    temperature_bin = floor(temperature / 1.25) * 1.25 + .625,
    altitude_bin = floor(altitude / 100) * 100 + 50
  ) |>
  count(temperature_bin, altitude_bin, name = "county_weeks")

# Panel b source data
curves <- bind_rows(lapply(seq_len(nrow(anchors)), function(i) {
  spec <- anchors[i, ] |>
    left_join(support, by = c("altitude_anchor", "altitude"))
  temperatures <- sort(unique(c(
    seq(
      floor(spec$p05[[1]] * 4) / 4,
      ceiling(spec$p95[[1]] * 4) / 4,
      by = .25
    ),
    reference_temperature
  )))
  basis <- basis_cumulative(temperatures)
  bind_cols(
    tibble(
      altitude_anchor = spec$altitude_anchor[[1]],
      altitude = spec$altitude[[1]],
      support_p05 = spec$p05[[1]],
      support_p95 = spec$p95[[1]],
      temperature = as.numeric(rownames(basis)),
      reference_temperature = reference_temperature
    ),
    estimate_basis(basis, spec$altitude_1000_c[[1]])
  )
}))

# Panel c source data: exact legacy definition, updated exposure.
moving <- estimate_moving_modification(width = 5)

# Panel d source data: two fixed supported contrasts across altitude.
altitude_grid <- seq(
  floor(min(model_dat$altitude) / 25) * 25,
  ceiling(max(model_dat$altitude) / 25) * 25,
  by = 25
)
fixed_specs <- tribble(
  ~contrast, ~colder_temperature, ~warmer_temperature,
  "-5 vs 0 C", -5, 0,
  "0 vs 5 C", 0, 5
)
fixed <- bind_rows(lapply(seq_len(nrow(fixed_specs)), function(i) {
  bind_cols(
    fixed_specs[i, ][rep(1, length(altitude_grid)), ],
    tibble(altitude = altitude_grid),
    estimate_contrast_at_altitude(
      fixed_specs$colder_temperature[[i]],
      fixed_specs$warmer_temperature[[i]],
      altitude_grid
    )
  )
}))

write_csv_utf8(density_dat, file.path(source_dir, "source_result2_panel_a_density.csv"))
write_csv_utf8(curves, file.path(source_dir, "source_result2_panel_b_curves.csv"))
write_csv_utf8(moving, file.path(source_dir, "source_result2_panel_c_5C_per1000m.csv"))
write_csv_utf8(fixed, file.path(source_dir, "source_result2_panel_d_fixed_contrasts.csv"))

theme_paper <- function(base_size = 6.7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(
        face = "bold", colour = dark, size = rel(1.06),
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        colour = muted, size = rel(.87), margin = margin(b = 5)
      ),
      axis.title = element_text(colour = dark),
      axis.text = element_text(colour = dark),
      axis.line = element_line(linewidth = .35, colour = dark),
      axis.ticks = element_line(linewidth = .35, colour = dark),
      legend.text = element_text(size = 5.3),
      legend.key.height = unit(2.8, "mm"),
      legend.margin = margin(0, 0, 1, 0),
      plot.tag = element_text(face = "bold", size = 7.7),
      panel.grid = element_blank()
    )
}

p_a_base <- ggplot(
  density_dat,
  aes(temperature_bin, altitude_bin, fill = county_weeks)
) +
  geom_tile(width = 1.24, height = 96, colour = "white", linewidth = .08) +
  scale_fill_gradient(
    low = tile_low, high = tile_high, name = "County-weeks",
    guide = guide_colorbar(
      direction = "horizontal", title.position = "left",
      title.hjust = .5, label.position = "bottom",
      barwidth = unit(16, "mm"), barheight = unit(1.5, "mm"),
      ticks = TRUE
    )
  ) +
  scale_x_continuous(breaks = seq(-20, 20, 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(2000, 4000, 500), expand = c(.01, .01)) +
  labs(
    x = paste0("Weekly mean air temperature (", degree_c, ")"),
    y = "County altitude (m)"
  ) +
  theme_paper() +
  theme(
    plot.margin = margin(0, 4, 1, 2)
  )

p_a_legend <- cowplot::get_legend(
  p_a_base +
    theme(
      legend.position = "top",
      legend.justification = "right",
      legend.box.just = "right",
      legend.direction = "horizontal",
      legend.title = element_text(size = 5.6, face = "bold"),
      legend.text = element_text(size = 5.1),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

p_a_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Observed 2-m air-temperature support across altitude",
    x = 0, y = .72, hjust = 0, vjust = .5,
    fontfamily = font_family, fontface = "bold", size = 7.1,
    colour = dark
  ) +
  cowplot::draw_label(
    "Binned county-week observations; darker tiles indicate more observations",
    x = 0, y = .18, hjust = 0, vjust = .5,
    fontfamily = font_family, size = 5.8, colour = muted
  )

p_a_header <- cowplot::plot_grid(
  p_a_title, p_a_legend,
  nrow = 1, rel_widths = c(.73, .27), align = "h", axis = "tb"
)

p_a_plot <- p_a_base +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 4, 1, 2)
  )

p_a <- cowplot::plot_grid(
  p_a_header, p_a_plot,
  ncol = 1, rel_heights = c(.20, 1), align = "v", axis = "lr"
)

curves <- curves |>
  mutate(
    altitude_anchor = factor(altitude_anchor, levels = c("P10", "P50", "P90")),
    slice = factor(
      altitude_anchor,
      levels = c("P10", "P50", "P90"),
      labels = paste0(
        format(round(support$altitude), big.mark = ","),
        " m model slice\nreference: 0 ", degree_c
      )
    )
  )

p_b <- ggplot(curves, aes(temperature, estimate)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = blue_light, alpha = .82) +
  geom_hline(yintercept = 1, linetype = 2, colour = null_grey, linewidth = .35) +
  geom_vline(
    xintercept = reference_temperature,
    linetype = 3, colour = dark, linewidth = .35
  ) +
  geom_line(colour = blue, linewidth = .74) +
  facet_wrap(~slice, nrow = 1, scales = "free_x") +
  scale_y_log10(
    breaks = c(.3, .5, 1, 2, 3),
    labels = c("0.3", "0.5", "1.0", "2.0", "3.0")
  ) +
  coord_cartesian(ylim = c(.25, 3.2)) +
  labs(
    title = "Nonlinear absolute-temperature responses across altitude",
    subtitle = paste0(
      "Cumulative lag 0\u20132 OR; curves restricted to local P5\u2013P95 ",
      "temperature support"
    ),
    x = paste0("Weekly mean air temperature (", degree_c, ")"),
    y = "Cumulative OR"
  ) +
  theme_paper() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7.5),
    panel.spacing.x = unit(3.0, "mm")
  )

outer_moving <- moving |> filter(common_p05_p95)
core_moving <- moving |> filter(common_p10_p90)

p_c <- ggplot(moving, aes(colder_temperature, estimate_per_1000m)) +
  geom_ribbon(
    data = outer_moving,
    aes(ymin = low_per_1000m, ymax = high_per_1000m),
    fill = blue_light, alpha = .62
  ) +
  geom_hline(yintercept = 1, linetype = 2, colour = null_grey, linewidth = .35) +
  geom_line(
    aes(colour = "Extrapolation", linetype = "Extrapolation"),
    linewidth = .56
  ) +
  geom_line(
    data = outer_moving,
    aes(colour = "P05-P95", linetype = "P05-P95"),
    linewidth = .66
  ) +
  geom_line(
    data = core_moving,
    aes(colour = "P10-P90", linetype = "P10-P90"),
    linewidth = .80
  ) +
  scale_colour_manual(
    values = c("P10-P90" = blue, "P05-P95" = blue_mid, "Extrapolation" = grey),
    breaks = c("P10-P90", "P05-P95", "Extrapolation"),
    labels = c("Common P10\u2013P90", "Common P5\u2013P95", "Extrapolation"), name = NULL
  ) +
  scale_linetype_manual(
    values = c("P10-P90" = "solid", "P05-P95" = "solid", "Extrapolation" = "dashed"),
    breaks = c("P10-P90", "P05-P95", "Extrapolation"),
    labels = c("Common P10\u2013P90", "Common P5\u2013P95", "Extrapolation"), name = NULL
  ) +
  scale_x_continuous(limits = c(-15, 10), breaks = seq(-15, 10, 5), expand = c(0, 0)) +
  scale_y_log10(
    breaks = c(.65, .8, 1, 1.2, 1.4),
    labels = c("0.65", "0.80", "1.00", "1.20", "1.40")
  ) +
  coord_cartesian(ylim = c(.60, 1.42)) +
  labs(
    title = paste0("Altitude modification across moving 5 ", degree_c, " contrasts"),
    subtitle = "ROR per 1,000-m increase; global interaction P = 0.602",
    x = paste0("Lower temperature in each 5 ", degree_c, " comparison (", degree_c, ")"),
    y = "ROR per 1,000-m increase"
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 1, nrow = 1)
  ) +
  theme_paper() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.key.width = unit(6.0, "mm")
  )

fixed <- fixed |>
  mutate(
    contrast_display = factor(
      contrast,
      levels = c("-5 vs 0 C", "0 vs 5 C"),
      labels = c(paste0("\u22125 vs 0 ", degree_c), paste0("0 vs 5 ", degree_c))
    )
  )

p_d <- ggplot(
  fixed,
  aes(altitude, estimate, colour = contrast_display, fill = contrast_display)
) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = .16, colour = NA) +
  geom_hline(yintercept = 1, linetype = 2, colour = null_grey, linewidth = .35) +
  geom_line(linewidth = .76) +
  scale_colour_manual(
    values = c(blue, orange), name = NULL
  ) +
  scale_fill_manual(
    values = c(blue_light, orange_light), name = NULL
  ) +
  scale_x_continuous(breaks = seq(2000, 4000, 500), expand = expansion(mult = c(.01, .01))) +
  scale_y_log10(
    breaks = c(.5, .7, 1, 1.3, 1.7),
    labels = c("0.5", "0.7", "1.0", "1.3", "1.7")
  ) +
  coord_cartesian(ylim = c(.48, 1.75)) +
  labs(
    title = "Prespecified temperature contrasts across altitude",
    subtitle = "Cumulative lag 0\u20132 ORs; both comparisons lie within common P10\u2013P90 support",
    x = "County altitude (m)", y = "Cumulative OR"
  ) +
  guides(
    fill = "none",
    colour = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(linewidth = .78)
    )
  ) +
  theme_paper() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.key.width = unit(6.0, "mm")
  )

# cowplot is used for the final assembly because the locally installed
# patchwork/grid combination clips lower-row axes when a faceted plot is placed
# above two independent plots. The individual ggplot objects are unchanged.
bottom_row <- cowplot::plot_grid(
  p_c, p_d, nrow = 1,
  labels = c("c", "d"),
  label_fontfamily = font_family,
  label_fontface = "bold", label_size = 7.7,
  label_x = c(0.003, 0.003), label_y = c(0.995, 0.995),
  hjust = 0, vjust = 1
)

body <- cowplot::plot_grid(
  p_a, p_b, bottom_row, ncol = 1,
  labels = c("a", "b", ""),
  label_fontfamily = font_family,
  label_fontface = "bold", label_size = 7.7,
  label_x = 0.003, label_y = 0.995,
  hjust = 0, vjust = 1,
  rel_heights = c(0.84, 1.33, 1.22)
)

figure2 <- body

stem <- file.path(figure_dir, "Figure2_AGU_no_title_20261008")
width_mm <- 183
height_mm <- 190
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

grDevices::png(
  paste0(stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 600, bg = "white", type = "cairo-png"
)
print(figure2)
grDevices::dev.off()
grDevices::tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw", bg = "white"
)
print(figure2)
grDevices::dev.off()
ggsave(
  paste0(stem, ".pdf"), figure2,
  width = width_mm, height = height_mm, units = "mm",
  device = cairo_pdf, bg = "white"
)
svglite::svglite(
  paste0(stem, ".svg"), width = width_in, height = height_in,
  bg = "white", fix_text_size = FALSE
)
print(figure2)
grDevices::dev.off()

reference_rows <- curves |>
  filter(abs(temperature - reference_temperature) < 1e-9)
qa <- tibble(
  check = c(
    "Full-sample model used",
    "No spring 2023 event term",
    "All three altitude slices present",
    "Panel b reference OR equals 1",
    "Panel c per-1000 m standardisation is reversible",
    "Panel c contains all three support classes",
    "Panel d temperatures share common P10-P90 support",
    "All plotted estimates finite"
  ),
  pass = c(
    nobs(fit) == nrow(model_dat),
    !grepl("spring2023_event", paste(deparse(formula(fit)), collapse = " ")),
    n_distinct(curves$altitude_anchor) == 3,
    nrow(reference_rows) == 3 &
      max(abs(log(reference_rows$estimate))) < 1e-8,
    max(abs(
      moving$estimate_per_1000m^altitude_span_km -
        moving$raw_p90_p10_ratio
    )) < 1e-10,
    n_distinct(moving$support_class) == 3,
    all(c(-5, 0, 5) >= common_p10 & c(-5, 0, 5) <= common_p90),
    all(is.finite(unlist(curves[c("estimate", "low", "high")]))) &
      all(is.finite(unlist(moving[c(
        "estimate_per_1000m", "low_per_1000m", "high_per_1000m"
      )]))) &
      all(is.finite(unlist(fixed[c("estimate", "low", "high")])) )
  )
)
write_csv_utf8(qa, file.path(qa_dir, "Result2_T2m_FINAL_qa.csv"))
if (!all(qa$pass)) stop("Unified Result 2 QA failed.")

legend_lines <- c(
  paste0(
    "Figure 2. Ambient-temperature exposure, response and altitude modification. ",
    "All panels use weekly mean 2-m air temperature and the full 2017-2024 sample."
  ),
  paste0(
    "a, Distribution of county-week observations across county altitude and weekly mean ",
    "2-m air temperature; darker tiles indicate more observations."
  ),
  paste0(
    "b, Cumulative lag 0-2 temperature-response curves from the hierarchical ",
    "beta-binomial distributed lag nonlinear model at the P10, P50 and P90 continuous-altitude ",
    "display anchors. The spring 2023 epidemic was retained without an event-specific ",
    "indicator. Curves are referenced to 0 ", degree_c,
    " and restricted to local P5-P95 temperature support."
  ),
  paste0(
    "c, Altitude modification across moving 5 ", degree_c,
    " comparisons, expressed as the ratio of odds ratios (ROR) per 1,000-m increase in altitude. ",
    "Dark and light solid segments denote common P10-P90 and P5-P95 support, respectively; ",
    "dashed grey segments denote extrapolation. Shading shows 95% confidence intervals within ",
    "common P5-P95 support. The global temperature-by-altitude interaction was assessed using ",
    "a 12-df joint Wald test (P=0.602)."
  ),
  paste0(
    "d, Cumulative lag 0-2 ORs for the prespecified -5 ", degree_c,
    " versus 0 ", degree_c, " and 0 ", degree_c, " versus 5 ", degree_c,
    " comparisons across continuous altitude. The warmer temperature is the reference in each ",
    "comparison; both lie within common P10-P90 temperature support. Shading shows 95% ",
    "confidence intervals."
  ),
  "P10, P50 and P90 are display anchors on the continuous altitude scale, not altitude strata."
)
writeLines(
  legend_lines,
  file.path(figure_dir, "Figure_2_AGU_legend_20261008.txt"),
  useBytes = TRUE
)

message("Created unified Result 2 with legacy-method Panels c and d.")
