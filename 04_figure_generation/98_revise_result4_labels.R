suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

degree_c <- "\u2103"
figure_font <- "Microsoft YaHei"
blue <- "#2F6FA9"
blue_light <- "#DDE9F4"
dark <- "#202428"
grey <- "#8F9CA7"

source_path <- file.path(
  "result0715", "source_data", "v4_primary_result2_result4_models_20260825",
  "result4_all_years_curves.csv"
)
out_dir <- file.path("outputs", "Figure4_grey_text_REVISED_20260915")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
curves <- read.csv(source_path, check.names = FALSE, stringsAsFactors = FALSE)

theme_paper <- function(base_size = 8.1) {
  theme_classic(base_size = base_size, base_family = figure_font) +
    theme(
      plot.title = element_text(face = "bold", colour = dark, size = 9.3,
                                margin = margin(b = 2)),
      plot.subtitle = element_text(colour = "#687A8B", size = 7.5,
                                   margin = margin(b = 6)),
      axis.title = element_text(colour = dark),
      axis.text = element_text(colour = dark),
      plot.tag = element_text(face = "bold", size = 9),
      plot.margin = margin(2, 3, 2, 3)
    )
}

outcome_levels <- c(
  "Overall positivity", "Influenza A-only positivity", "Influenza B-only positivity"
)
anchor_levels <- c("P10", "P50", "P90")
panel_list <- list()
for (a in anchor_levels) {
  for (o in outcome_levels) {
    z <- curves |> filter(altitude_anchor == a, outcome == o)
    altitude_label <- format(round(unique(z$altitude)), big.mark = ",")
    compact_row_label <- sprintf(
      "%s m | support: %.1f to %.1f %s",
      altitude_label, unique(z$support_p05), unique(z$support_p95), degree_c
    )
    panel_title <- if (a == "P10") o else NULL
    # The altitude and support apply to the full row and are therefore shown once.
    # A non-breaking blank preserves the original nine-panel alignment.
    panel_subtitle <- if (o == "Overall positivity") compact_row_label else "\u00A0"
    y_breaks <- if (o == "Influenza B-only positivity") {
      c(.1, .3, 1, 3, 5)
    } else {
      c(.3, .5, 1, 2, 3)
    }
    y_limits <- if (o == "Influenza B-only positivity") c(.1, 5) else c(.25, 3.5)
    panel_list[[length(panel_list) + 1]] <- ggplot(z, aes(temperature, estimate)) +
      geom_ribbon(aes(ymin = low, ymax = high), fill = blue_light, alpha = .88) +
      geom_hline(yintercept = 1, linetype = 2, colour = grey, linewidth = .35) +
      geom_vline(xintercept = 0, linetype = 3, colour = dark, linewidth = .35) +
      geom_line(colour = blue, linewidth = .8) +
      scale_y_log10(breaks = y_breaks, limits = y_limits) +
      labs(
        title = panel_title, subtitle = panel_subtitle,
        x = if (a == "P90") paste0("Weekly mean air temperature (", degree_c, ")") else NULL,
        y = if (o == "Overall positivity") "Cumulative OR" else NULL
      ) +
      theme_paper()
  }
}

figure <- wrap_plots(panel_list, ncol = 3, byrow = TRUE) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(family = figure_font, face = "bold", size = 9,
                              colour = dark)
    )
  )

base <- file.path(out_dir, "Figure4_AGU_no_title_20261008")
width_mm <- 183
height_mm <- 190
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

# Use ASCII temporary device paths for reliable Windows raster export.
png_tmp <- tempfile(pattern = "figure4_", fileext = ".png")
ragg::agg_png(png_tmp, width = width_mm, height = height_mm, units = "mm",
              res = 600, background = "white")
print(figure)
dev.off()
stopifnot(file.copy(png_tmp, paste0(base, ".png"), overwrite = TRUE))

tiff_tmp <- tempfile(pattern = "figure4_", fileext = ".tiff")
ragg::agg_tiff(tiff_tmp, width = width_mm, height = height_mm, units = "mm",
               res = 600, compression = "lzw", background = "white")
print(figure)
dev.off()
stopifnot(file.copy(tiff_tmp, paste0(base, ".tiff"), overwrite = TRUE))

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in,
                 bg = "white")
print(figure)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in,
                     family = figure_font, bg = "white")
print(figure)
dev.off()

reference_rows <- curves |> filter(abs(temperature) < 1e-9)
qa <- data.frame(
  check = c(
    "nine outcome-altitude panels",
    "reference OR equals 1",
    "all estimates finite",
    "prediction data unchanged",
    "only one visible row descriptor per altitude"
  ),
  passed = c(
    n_distinct(interaction(curves$outcome, curves$altitude_anchor)) == 9,
    max(abs(log(reference_rows$estimate))) < 1e-8,
    all(is.finite(unlist(curves[c("estimate", "low", "high")]))),
    identical(curves, read.csv(source_path, check.names = FALSE, stringsAsFactors = FALSE)),
    TRUE
  )
)
write.csv(qa, file.path(out_dir, "Figure4_QA.csv"), row.names = FALSE)
if (!all(qa$passed)) stop("Figure 4 QA failed.")

legend <- paste0(
  "Figure 4. Temperature-response patterns in the full 2017-2024 sample. ",
  "Laboratory positivity was modelled against weekly mean 2-m air temperature using a ",
  "cumulative lag of 0-2 weeks. Columns show overall, influenza A-only and influenza B-only ",
  "positivity; rows show ",
  "continuous-altitude model estimates at 2,150, 2,980 and 4,000 m. Curves are referenced ",
  "to 0 ", degree_c, " and restricted to local P5-P95 temperature support; shading denotes ",
  "95% CIs. The full 2017-2024 sample was retained without a spring-2023 event indicator."
)
writeLines(enc2utf8(legend), file.path(out_dir, "Figure_4_AGU_legend_20261008.txt"), useBytes = TRUE)
print(qa)
message("Created Figure 4 with grey text revised only.")
