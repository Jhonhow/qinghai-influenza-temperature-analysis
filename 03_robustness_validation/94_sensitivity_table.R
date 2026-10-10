suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(svglite)
  library(ragg)
})

workspace <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(workspace, "outputs", "result5_temporal_supplement_20260902")
source_dir <- file.path(output_dir, "source_data")

phase_results <- read_csv(
  file.path(source_dir, "table_s4_china_phase_altitude_results.csv"),
  show_col_types = FALSE
)
phase_heterogeneity <- read_csv(
  file.path(source_dir, "table_s4_china_phase_heterogeneity.csv"),
  show_col_types = FALSE
)
event_results <- read_csv(
  file.path(source_dir, "table_s4_spring2023_event_sensitivity.csv"),
  show_col_types = FALSE
)

colours <- c(
  ink = "#17212B",
  text = "#303841",
  muted = "#667786",
  rule = "#D8E0E6",
  grey = "#56616C",
  grey_mid = "#97A3AD",
  grey_pale = "#F1F4F6",
  blue = "#2F73AE",
  blue_pale = "#EAF2F8",
  orange = "#C78A30",
  red = "#C85647"
)
font_family <- "Microsoft YaHei"
degree_c <- "\u2103"

theme_r5_refined <- function() {
  theme_classic(base_size = 8.0, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = colours[["text"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = colours[["text"]]),
      axis.title = element_text(size = 7.6, colour = colours[["ink"]]),
      axis.text = element_text(size = 7.0, colour = colours[["text"]]),
      plot.title = element_text(size = 8.7, face = "bold", colour = colours[["ink"]]),
      plot.subtitle = element_text(size = 6.9, colour = colours[["muted"]]),
      plot.margin = margin(3, 3, 3, 3),
      panel.grid = element_blank()
    )
}

timeline <- tibble::tribble(
  ~track, ~label, ~start, ~end, ~segment,
  "Primary periods", "Pre-pandemic", as.Date("2017-01-01"), as.Date("2019-12-31"), "Pre-pandemic",
  "Primary periods", "Pandemic disruption", as.Date("2020-01-01"), as.Date("2022-12-31"), "Pandemic disruption",
  "Primary periods", "Post-transition", as.Date("2023-01-01"), as.Date("2024-12-31"), "Post-transition",
  "Spring 2023 epidemic", "Feb-May 2023", as.Date("2023-02-01"), as.Date("2023-05-31"), "Spring 2023 epidemic",
  "WHO PHEIC", "WHO PHEIC", as.Date("2020-01-30"), as.Date("2023-05-05"), "WHO PHEIC"
) |>
  mutate(
    track = factor(track, levels = c("WHO PHEIC", "Spring 2023 epidemic", "Primary periods")),
    midpoint = start + (end - start) / 2
  )

timeline_colours <- c(
  "Pre-pandemic" = colours[["grey_mid"]],
  "Pandemic disruption" = colours[["orange"]],
  "Post-transition" = colours[["blue"]],
  "Spring 2023 epidemic" = colours[["red"]],
  "WHO PHEIC" = colours[["grey"]]
)

p_a <- ggplot() +
  geom_segment(
    data = timeline |> filter(track == "Primary periods"),
    aes(x = start, xend = end, y = track, yend = track, colour = segment),
    linewidth = 6.4, lineend = "butt"
  ) +
  geom_segment(
    data = timeline |> filter(track != "Primary periods"),
    aes(x = start, xend = end, y = track, yend = track, colour = segment),
    linewidth = 4.6, lineend = "butt"
  ) +
  geom_text(
    data = timeline |> filter(track == "Primary periods"),
    aes(x = midpoint, y = track, label = label),
    colour = "white", size = 2.0, fontface = "bold"
  ) +
  geom_text(
    data = timeline |> filter(segment == "Spring 2023 epidemic"),
    aes(x = end + 45, y = track, label = "Feb-May 2023"),
    hjust = 0, colour = colours[["red"]], size = 2.0, fontface = "bold"
  ) +
  scale_colour_manual(values = timeline_colours, guide = "none") +
  scale_x_date(
    date_breaks = "1 year", date_labels = "%Y",
    limits = as.Date(c("2016-11-01", "2025-02-01")),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "a   Prespecified windows",
    subtitle = "Primary periods, spring 2023 epidemic and WHO PHEIC",
    x = NULL, y = NULL
  ) +
  theme_r5_refined() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 6.8),
    plot.margin = margin(3, 5, 5, 3)
  )

phase_het_p <- phase_heterogeneity |>
  filter(estimand == "Altitude modification of -5 vs 0 C") |>
  pull(p_value)

phase_plot <- phase_results |>
  mutate(
    y = recode(
      phase,
      "Pre-pandemic (2017-2019)" = 3,
      "Pandemic disruption (2020-2022)" = 2,
      "Post-transition (2023-2024)" = 1
    ),
    highlight = if_else(grepl("Post-transition", phase), "Post-transition", "Other periods"),
    estimate_label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, low, high)
  )

p_b_core <- ggplot(phase_plot, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.42, colour = colours[["grey_mid"]]) +
  geom_errorbar(
    aes(xmin = low, xmax = high, colour = highlight),
    orientation = "y", width = 0, linewidth = 0.8
  ) +
  geom_point(
    aes(colour = highlight), shape = 21, fill = "white",
    stroke = 0.85, size = 2.65
  ) +
  scale_colour_manual(
    values = c("Post-transition" = colours[["blue"]], "Other periods" = colours[["grey"]]),
    guide = "none"
  ) +
  scale_x_log10(
    limits = c(0.55, 2.30), breaks = c(0.6, 0.8, 1, 1.5, 2),
    labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = c(3, 2, 1),
    labels = c("Pre-pandemic\n2017-2019", "Pandemic disruption\n2020-2022", "Post-transition\n2023-2024"),
    limits = c(0.55, 3.45), expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "b   Period-specific modification",
    subtitle = sprintf("Period heterogeneity: P = %.3f", phase_het_p),
    x = "Ratio of odds ratios per 1,000 m", y = NULL
  ) +
  theme_r5_refined() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

p_b_text <- ggplot(phase_plot, aes(x = 0, y = y, label = estimate_label)) +
  geom_text(hjust = 0, size = 2.2, colour = colours[["text"]]) +
  annotate(
    "text", x = 0, y = 3.40, label = "Estimate (95% CI)",
    hjust = 0, vjust = 1, size = 2.05, fontface = "bold", colour = colours[["muted"]]
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0.55, 3.45), expand = expansion(mult = c(0, 0))) +
  theme_void(base_family = font_family) +
  theme(plot.margin = margin(3, 2, 20, 1))

b_block <- (p_b_core | p_b_text) +
  plot_layout(widths = c(1.65, 0.88))

event_plot <- event_results |>
  mutate(
    scope_short = recode(
      scope,
      "Post-transition 2023-2024" = "Post-transition 2023-2024",
      "Full 2017-2024 sample" = "Full 2017-2024 sample"
    ),
    strategy = recode(
      event_handling,
      "Peak retained" = "Retained, unadjusted",
      "Event indicator adjusted" = "Retained + event indicator",
      "February-May 2023 excluded" = "Feb-May excluded"
    ),
    y = case_when(
      scope_short == "Post-transition 2023-2024" & strategy == "Retained, unadjusted" ~ 6,
      scope_short == "Post-transition 2023-2024" & strategy == "Retained + event indicator" ~ 5,
      scope_short == "Post-transition 2023-2024" & strategy == "Feb-May excluded" ~ 4,
      scope_short == "Full 2017-2024 sample" & strategy == "Retained, unadjusted" ~ 2,
      scope_short == "Full 2017-2024 sample" & strategy == "Retained + event indicator" ~ 1,
      TRUE ~ 0
    ),
    group = if_else(grepl("Post-transition", scope_short), "Post-transition", "Full sample"),
    estimate_label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, low, high)
  )

p_c_core <- ggplot(event_plot, aes(x = estimate, y = y)) +
  annotate("rect", xmin = 0.55, xmax = 2.30, ymin = 3.55, ymax = 7.45, fill = colours[["blue_pale"]]) +
  annotate("rect", xmin = 0.55, xmax = 2.30, ymin = -0.45, ymax = 3.35, fill = colours[["grey_pale"]]) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.42, colour = colours[["grey_mid"]]) +
  geom_errorbar(
    aes(xmin = low, xmax = high, colour = group),
    orientation = "y", width = 0, linewidth = 0.78
  ) +
  geom_point(
    aes(colour = group), shape = 21, fill = "white",
    stroke = 0.85, size = 2.55
  ) +
  annotate(
    "text", x = 0.57, y = 7.15, label = "Post-transition",
    hjust = 0, size = 2.1, fontface = "bold", colour = colours[["blue"]]
  ) +
  annotate(
    "text", x = 0.57, y = 3.05, label = "Full sample",
    hjust = 0, size = 2.1, fontface = "bold", colour = colours[["grey"]]
  ) +
  scale_colour_manual(
    values = c("Post-transition" = colours[["blue"]], "Full sample" = colours[["grey"]]),
    guide = "none"
  ) +
  scale_x_log10(
    limits = c(0.55, 2.30), breaks = c(0.6, 0.8, 1, 1.5, 2),
    labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = c(6, 5, 4, 2, 1, 0),
    labels = c(
      "Retained, unadjusted", "Retained + indicator", "Feb-May excluded",
      "Retained, unadjusted", "Retained + indicator", "Feb-May excluded"
    ),
    limits = c(-0.45, 7.45), expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "c   Sensitivity to spring 2023",
    subtitle = "Three event-handling strategies",
    x = "Ratio of odds ratios per 1,000 m", y = NULL
  ) +
  theme_r5_refined() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

p_c_text <- ggplot(event_plot, aes(x = 0, y = y, label = estimate_label)) +
  geom_text(hjust = 0, size = 2.15, colour = colours[["text"]]) +
  annotate(
    "text", x = 0, y = 7.25, label = "Estimate (95% CI)",
    hjust = 0, size = 2.05, fontface = "bold", colour = colours[["muted"]]
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(-0.45, 7.45), expand = expansion(mult = c(0, 0))) +
  theme_void(base_family = font_family) +
  theme(plot.margin = margin(3, 2, 20, 1))

c_block <- (p_c_core | p_c_text) +
  plot_layout(widths = c(1.7, 0.9))

figure_design <- "
AA
BC
"

figure_s4 <- p_a + b_block + c_block +
  plot_layout(
    design = figure_design,
    heights = c(0.62, 1.38),
    widths = c(0.91, 1.09)
  ) +
  plot_annotation(
    title = "Temporal robustness of altitude modification",
    subtitle = paste0(
      "Overall laboratory positivity; weekly mean 2-m air temperature; cumulative lag 0-2; ",
      "contrast: -5 ", degree_c, " versus 0 ", degree_c
    ),
    theme = theme(
      plot.title = element_text(family = font_family, size = 10.8, face = "bold", colour = colours[["ink"]]),
      plot.subtitle = element_text(family = font_family, size = 7.2, colour = colours[["muted"]]),
      plot.margin = margin(4, 5, 3, 5)
    )
  )

base <- file.path(output_dir, "Supplementary_Figure_S5_Temporal_Robustness_FINAL_20260910")
width_in <- 183 / 25.4
height_in <- 126 / 25.4

temp_svg <- file.path(tempdir(), "Supplementary_Figure_S5_Temporal_Robustness_FINAL_20260910.svg")
svglite::svglite(temp_svg, width = width_in, height = height_in)
print(figure_s4)
dev.off()
stopifnot(file.copy(temp_svg, paste0(base, ".svg"), overwrite = TRUE))

temp_pdf <- file.path(tempdir(), "Supplementary_Figure_S5_Temporal_Robustness_FINAL_20260910.pdf")
grDevices::cairo_pdf(
  temp_pdf, width = width_in, height = height_in, family = "Arial"
)
print(figure_s4)
dev.off()
stopifnot(file.copy(temp_pdf, paste0(base, ".pdf"), overwrite = TRUE))

temp_tiff <- file.path(tempdir(), "Supplementary_Figure_S5_Temporal_Robustness_FINAL_20260910.tiff")
ragg::agg_tiff(
  temp_tiff, width = width_in, height = height_in, units = "in", res = 600,
  compression = "lzw"
)
print(figure_s4)
dev.off()
stopifnot(file.copy(temp_tiff, paste0(base, ".tiff"), overwrite = TRUE))

temp_png <- file.path(tempdir(), "Supplementary_Figure_S5_Temporal_Robustness_FINAL_20260910.png")
ragg::agg_png(temp_png, width = width_in, height = height_in, units = "in", res = 300)
print(figure_s4)
dev.off()
stopifnot(file.copy(temp_png, paste0(base, ".png"), overwrite = TRUE))

write_csv(
  bind_rows(
    phase_plot |>
      transmute(panel = "b", group = phase, analysis = "Period-specific", estimate, low, high),
    event_plot |>
      transmute(panel = "c", group = scope_short, analysis = event_handling, estimate, low, high)
  ),
  file.path(source_dir, "Supplementary_Figure_S5_Temporal_Robustness_source_data.csv")
)

message("Refined Supplementary Figure S5 exported successfully.")
