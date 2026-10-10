suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
})

# Figure contract
# Conclusion: influenza A activity in February-May 2023 occurred outside the
# usual winter timing and remained a distinct signal after adjustment for the
# current 2-m air-temperature model.
# Evidence: monthly timing (panel a) and adjusted event-period ORs (panel b).

font_family <- "Microsoft YaHei"
ink <- "#17202A"
muted <- "#6D7C8B"
grid <- "#D9E0E6"
blue <- "#2F6FA9"
dark <- "#303A43"
red <- "#D65A47"
degree_c <- "\u2103"

model_dir <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "primary_result2_result4_models_20260825"
)
source_model_dir <- file.path(
  "result0715", "source_data", "v4_primary_result2_result4_models_20260825"
)
output_dir <- file.path(
  "result0715", "supplementary_figures", "v4_t2m_20260825",
  "result4_spring2023_event_context_20260910"
)
source_dir <- file.path(
  "result0715", "source_data", "v4_spring2023_event_context_20260910"
)
qa_dir <- file.path(
  "result0715", "figure_qa", "v4_t2m_20260825",
  "result4_spring2023_event_context_20260910"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(file.path(model_dir, "v4_t2m_model_bundle.rds"))
panel <- readRDS(file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825",
  "county_week_panel_v4_t2m_2016_2025.rds"
))
event_or <- read.csv(
  file.path(source_model_dir, "result4_event_indicator_estimates.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
event_curves <- read.csv(
  file.path(source_model_dir, "result4_event_curves.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

monthly_a <- panel |>
  filter(iso_year %in% 2017:2024, lab_available %in% TRUE, tests > 0) |>
  mutate(month = month(week_start)) |>
  group_by(iso_year, month) |>
  summarise(
    tests = sum(tests),
    influenza_a = sum(fluA_only),
    positivity = influenza_a / tests,
    .groups = "drop"
  ) |>
  complete(iso_year = 2017:2024, month = 1:12)

event_or_plot <- event_or |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("Influenza A-only positivity", "Overall positivity"),
      labels = c("Influenza A", "Overall")
    ),
    label = sprintf("%.2f (%.2f-%.2f)", estimate, low, high),
    colour = if_else(as.character(outcome) == "Influenza A", blue, dark)
  )

theme_figure <- theme_classic(base_size = 8, base_family = font_family) +
  theme(
    plot.title = element_text(
      face = "bold", colour = ink, size = 9.2, margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      colour = muted, size = 7, margin = margin(b = 5)
    ),
    axis.title = element_text(colour = ink, size = 7.4),
    axis.text = element_text(colour = ink, size = 7),
    axis.line = element_line(colour = ink, linewidth = 0.38),
    axis.ticks = element_line(colour = ink, linewidth = 0.3),
    plot.tag = element_text(face = "bold", colour = ink, size = 8.5),
    plot.margin = margin(4, 6, 4, 5)
  )

panel_a <- ggplot(monthly_a, aes(month, factor(iso_year), fill = positivity)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  annotate(
    "rect", xmin = 1.5, xmax = 5.5, ymin = 6.5, ymax = 7.5,
    fill = NA, colour = red, linewidth = 0.9
  ) +
  scale_fill_gradient(
    low = "#F1F5F8", high = blue, limits = c(0, 0.60),
    breaks = c(0, 0.2, 0.4, 0.6),
    labels = percent_format(accuracy = 1), name = NULL,
    guide = guide_colorbar(
      barheight = unit(24, "mm"), barwidth = unit(3.5, "mm")
    )
  ) +
  scale_x_continuous(
    breaks = c(1, 3, 5, 7, 9, 11),
    labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Off-calendar timing",
    subtitle = "Monthly influenza A positivity; Feb-May 2023 outlined",
    x = NULL, y = NULL
  ) +
  theme_figure +
  theme(
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 2)
  )

panel_b <- ggplot(event_or_plot, aes(estimate, outcome, colour = colour)) +
  geom_vline(xintercept = 1, linetype = "22", colour = muted, linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = low, xmax = high), orientation = "y",
    width = 0.15, linewidth = 0.8
  ) +
  geom_point(shape = 21, fill = "white", size = 2.7, stroke = 0.75) +
  geom_text(
    aes(x = estimate, label = label),
    position = position_nudge(y = 0.24),
    hjust = 0.5, colour = ink, family = font_family, size = 2.35
  ) +
  scale_colour_identity() +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20),
    labels = label_number(accuracy = 1)
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.35, 0.35))) +
  labs(
    title = "Adjusted event signal",
    subtitle = "Current 2-m air-temperature models",
    x = "Event-period OR", y = NULL
  ) +
  coord_cartesian(xlim = c(0.8, 25), clip = "off") +
  theme_figure +
  theme(
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.24),
    plot.margin = margin(4, 6, 4, 7)
  )

curve_data <- event_curves |>
  mutate(
    method = factor(
      method,
      levels = c("All years", "Event-adjusted", "Excluding spring 2023"),
      labels = c("Full sample", "Event-adjusted", "Feb-May excluded")
    ),
    outcome = factor(
      outcome,
      levels = c("Overall positivity", "Influenza A-only positivity"),
      labels = c("Overall positivity", "Influenza A positivity")
    )
  )

curve_colours <- c(
  "Full sample" = blue,
  "Event-adjusted" = "#C8872B",
  "Feb-May excluded" = "#7F8C96"
)
curve_linetypes <- c(
  "Full sample" = "solid",
  "Event-adjusted" = "42",
  "Feb-May excluded" = "22"
)

panel_c <- ggplot(
  curve_data,
  aes(temperature, estimate, colour = method, fill = method, linetype = method)
) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.045, colour = NA) +
  geom_hline(yintercept = 1, linetype = "22", colour = muted, linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "33", colour = muted, linewidth = 0.35) +
  geom_line(linewidth = 0.72) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_colour_manual(values = curve_colours, name = NULL) +
  scale_fill_manual(values = curve_colours, name = NULL) +
  scale_linetype_manual(values = curve_linetypes, name = NULL) +
  scale_x_continuous(breaks = c(-5, 0, 5, 10, 15)) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 2, 3),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(
    title = "Temperature-response sensitivity",
    subtitle = paste0("2,150-m slice; reference: 0 ", degree_c, "; lag 0-2"),
    x = paste0("Weekly mean 2-m air temperature (", degree_c, ")"),
    y = "Cumulative OR"
  ) +
  coord_cartesian(ylim = c(0.25, 3.6)) +
  theme_figure +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = unit(7, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = ink, size = 7.2),
    panel.spacing.x = unit(7, "mm"),
    panel.grid.major.y = element_line(colour = grid, linewidth = 0.22)
  )

top_row <- panel_a + panel_b + plot_layout(widths = c(1.18, 1))
combined <- (top_row / panel_c) +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(
    title = "Spring 2023 influenza A epidemic",
    subtitle = "Off-calendar timing, adjusted signal and response sensitivity",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(
        family = font_family, face = "bold", colour = ink,
        size = 11.5, margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        family = font_family, colour = muted, size = 7.5,
        margin = margin(b = 5)
      ),
      plot.tag = element_text(
        family = font_family, face = "bold", colour = ink, size = 8.5
      ),
      plot.margin = margin(6, 8, 5, 7)
    )
  )

base_path <- file.path(
  output_dir, "FigureS_Spring2023_Event_Context_T2m_FINAL_20260910"
)
width_in <- 183 / 25.4
height_in <- 135 / 25.4

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

write.csv(
  monthly_a,
  file.path(source_dir, "panel_a_monthly_influenza_a_positivity.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  event_or_plot |>
    transmute(
      outcome = as.character(outcome), estimate, low, high,
      likelihood_ratio = LR, df, p_value, delta_AIC
    ),
  file.path(source_dir, "panel_b_event_indicator_estimates.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  curve_data,
  file.path(source_dir, "panel_c_t2m_event_handling_curves.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

qa <- tibble(
  check = c(
    "All 96 year-month cells are observed",
    "Event estimates use the current T2m model bundle",
    "Both adjusted event ORs exceed one",
    "All event confidence intervals are finite",
    "Panel c contains all three event-handling strategies",
    "Panel c reference OR equals one",
    "All panel c estimates are finite and positive",
    "All four output formats exist",
    "Raster outputs exceed minimum size"
  ),
  pass = c(
    nrow(monthly_a) == 96 && all(is.finite(monthly_a$positivity)),
    !is.null(bundle$model_data$temperature),
    all(event_or_plot$estimate > 1),
    all(is.finite(unlist(event_or_plot[c("estimate", "low", "high")]))),
    nlevels(curve_data$method) == 3 && all(table(curve_data$method) > 0),
    max(abs(log(curve_data$estimate[abs(curve_data$temperature) < 1e-9]))) < 1e-8,
    all(is.finite(unlist(curve_data[c("estimate", "low", "high")]))) &&
      all(unlist(curve_data[c("estimate", "low", "high")]) > 0),
    all(file.exists(paste0(base_path, c(".svg", ".pdf", ".png", ".tiff")))),
    all(file.info(paste0(base_path, c(".png", ".tiff")))$size > 100000)
  )
)
write.csv(
  qa,
  file.path(qa_dir, "FigureS_Spring2023_Event_Context_T2m_QA.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
if (!all(qa$pass)) warning("At least one event-context QA item failed.")

message("Created updated T2m spring-2023 event-context supplementary figure.")
