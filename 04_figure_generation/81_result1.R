suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(scales)
})

degree_c <- "\u00B0C"
figure_font <- "Microsoft YaHei"
main_threshold <- -10
thresholds <- c(-5, -10, -15)

model_root <- file.path("model_outputs", "reanalysis_v4_t2m_20260825")
panel_path <- file.path(model_root, "county_week_panel_v4_t2m_2016_2025.rds")
bundle_path <- file.path(
  model_root, "primary_result2_result4_models_20260825", "v4_t2m_model_bundle.rds"
)
shp_path <- Sys.getenv(
  "QH_COUNTY_SHAPEFILE",
  unset = file.path("data", "qinghai_shp", "china_county.shp")
)
figure_dir <- file.path("outputs", "result1_revision_20260908")
supp_dir <- file.path(figure_dir, "supplementary")
source_dir <- file.path(figure_dir, "source_data")
qa_dir <- file.path(figure_dir, "qa")
for (z in c(figure_dir, supp_dir, source_dir, qa_dir)) {
  dir.create(z, recursive = TRUE, showWarnings = FALSE)
}
if (!all(file.exists(c(panel_path, bundle_path, shp_path)))) stop("Missing Result 1 input.")

normalize_county_local <- function(x) {
  z <- trimws(as.character(x))
  aliases <- c(
    "\u5927\u901a\u53bf" = "\u5927\u901a\u56de\u65cf\u571f\u65cf\u81ea\u6cbb\u53bf",
    "\u4e50\u90fd\u53bf" = "\u4e50\u90fd\u533a",
    "\u95e8\u6e90\u53bf" = "\u95e8\u6e90\u56de\u65cf\u81ea\u6cbb\u53bf",
    "\u540c\u4ec1\u53bf" = "\u540c\u4ec1\u5e02",
    "\u6cb3\u5357\u53bf" = "\u6cb3\u5357\u8499\u53e4\u65cf\u81ea\u6cbb\u53bf",
    "\u5927\u67f4\u65e6\u884c\u59d4" = "\u5927\u67f4\u65e6\u884c\u653f\u59d4\u5458\u4f1a"
  )
  dplyr::recode(z, !!!aliases, .default = z)
}
write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

panel <- readRDS(panel_path) |> mutate(week_start = as.Date(week_start))
bundle <- readRDS(bundle_path)
primary <- bundle$model_data |> mutate(week_start = as.Date(week_start))
primary_counties <- unique(as.character(primary$county))

study_panel <- panel |>
  filter(county %in% primary_counties, iso_year >= 2017, iso_year <= 2024)
met_panel <- study_panel |>
  filter(met_valid_5day_main %in% TRUE, !is.na(temperature), !is.na(altitude))
lab_observed <- met_panel |>
  filter(lab_available %in% TRUE, tests >= 1, positive >= 0, positive <= tests)
ili_observed <- study_panel |>
  filter(ili_available %in% TRUE, ILI >= 0, visits > 0, ILI <= visits)

altitude_summary <- primary |> distinct(county, city, altitude)
threshold_long <- bind_rows(lapply(thresholds, function(z) {
  met_panel |>
    group_by(county) |>
    summarise(
      threshold = z, met_valid_weeks = n(), cold_weeks = sum(temperature <= z),
      cold_week_pct = 100 * cold_weeks / met_valid_weeks,
      mean_temperature = mean(temperature), .groups = "drop"
    )
}))
met_county <- threshold_long |> filter(threshold == main_threshold)
lab_county <- lab_observed |>
  group_by(county) |>
  summarise(
    lab_observed_weeks = n(), total_tests = sum(tests),
    total_positive = sum(positive), positivity = total_positive / total_tests,
    .groups = "drop"
  )
county_context <- altitude_summary |>
  left_join(met_county, by = "county") |>
  left_join(lab_county, by = "county") |>
  mutate(
    lab_observed_weeks = coalesce(lab_observed_weeks, 0L),
    lab_coverage_pct = 100 * lab_observed_weeks / met_valid_weeks
  )

qinghai_raw <- st_read(shp_path, quiet = TRUE, options = "ENCODING=UTF-8")
province_name <- "\u9752\u6d77\u7701"
qinghai <- qinghai_raw |>
  filter(.data[[names(qinghai_raw)[5]]] == province_name) |>
  mutate(map_county = normalize_county_local(.data[[names(qinghai_raw)[1]]])) |>
  st_make_valid() |>
  group_by(map_county) |>
  summarise(.groups = "drop")

unmatched <- anti_join(
  county_context |> select(county),
  st_drop_geometry(qinghai) |> select(map_county),
  by = c("county" = "map_county")
)
if (nrow(unmatched)) {
  write_csv_utf8(unmatched, file.path(source_dir, "unmatched_counties.csv"))
  stop("County map join failed.")
}

map_data <- qinghai |>
  left_join(county_context, by = c("map_county" = "county")) |>
  mutate(analytic_county = !is.na(altitude)) |>
  st_transform(4326)

threshold_map_data <- qinghai |>
  left_join(
    threshold_long |> left_join(altitude_summary, by = "county"),
    by = c("map_county" = "county")
  ) |>
  st_transform(4326)

study_start <- as.Date("2017-01-02")
study_end <- as.Date("2024-12-23")
week_grid <- tibble(week_start = seq(study_start, study_end, by = "7 days"))
lab_weekly <- lab_observed |>
  group_by(week_start) |>
  summarise(
    laboratory_tests = sum(tests),
    laboratory_positivity_pct = 100 * sum(positive) / sum(tests),
    .groups = "drop"
  )
ili_weekly <- ili_observed |>
  group_by(week_start) |>
  summarise(ILI_rate_pct = 100 * sum(ILI) / sum(visits), .groups = "drop")
weekly <- week_grid |>
  left_join(lab_weekly, by = "week_start") |>
  left_join(ili_weekly, by = "week_start")

palette <- c(
  ink = "#25282B", muted = "#6D7680", line = "#3C4650",
  grid = "#E9EDF0", outside = "#ECEFF1", blue = "#2C6DA4",
  blue_fill = "#DCEAF3", red = "#C85B4B", red_fill = "#F1DCD7",
  teal = "#27847F", teal_fill = "#DCECE9"
)
theme_map <- function(base_size = 7) {
  theme_void(base_size = base_size, base_family = figure_font) +
    theme(
      plot.title = element_text(face = "bold", colour = palette[["ink"]],
                                size = base_size + .8, margin = margin(b = 3)),
      legend.title = element_text(size = base_size - .1),
      legend.text = element_text(size = base_size - .5),
      axis.text.x = element_text(size = base_size - 1.0, colour = palette[["muted"]],
                                 margin = margin(t = 2)),
      axis.text.y = element_text(size = base_size - 1.0, colour = palette[["muted"]],
                                 margin = margin(r = 2)),
      axis.ticks = element_line(linewidth = .22, colour = palette[["muted"]]),
      axis.ticks.length = unit(1.4, "pt"),
      panel.grid.major = element_line(linewidth = .16, colour = "#DCE3E8"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 5, 2, 5)
    )
}
theme_ts <- function(base_size = 7, show_x = FALSE) {
  theme_classic(base_size = base_size, base_family = figure_font) +
    theme(
      axis.line = element_line(linewidth = .3, colour = palette[["line"]]),
      axis.ticks = element_line(linewidth = .25, colour = palette[["line"]]),
      axis.text = element_text(size = base_size - .6, colour = palette[["ink"]]),
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = .25) else element_blank(),
      axis.title = element_blank(),
      panel.grid.major.y = element_line(linewidth = .18, colour = palette[["grid"]]),
      plot.title = element_text(face = "bold", size = base_size + .8),
      plot.margin = margin(1, 5, 1, 5)
    )
}

base_map <- ggplot(map_data) +
  geom_sf(fill = palette[["outside"]], colour = "white", linewidth = .16)
map_axes <- list(
  coord_sf(datum = st_crs(4326), expand = FALSE,
           label_axes = list(bottom = "E", left = "N")),
  scale_x_continuous(
    breaks = c(90, 95, 100),
    labels = function(x) paste0(x, "\u00B0E")
  ),
  scale_y_continuous(
    breaks = c(32, 36, 40),
    labels = function(x) paste0(x, "\u00B0N")
  )
)
p_a <- base_map +
  geom_sf(aes(fill = altitude), colour = "white", linewidth = .16) +
  scale_fill_gradientn(
    colours = c("#B9E3CF", "#5CB7AC", "#3979AE", "#2B2342"),
    name = "Altitude (m)", breaks = c(2000, 3000, 4000),
    na.value = palette[["outside"]],
    guide = guide_colorbar(title.position = "top", barwidth = unit(26, "mm"),
                           barheight = unit(2.2, "mm"))
  ) + map_axes + labs(title = "a  Altitude of county government seats") + theme_map() +
  theme(legend.position = "bottom")

p_b <- base_map +
  geom_sf(aes(fill = cold_week_pct), colour = "white", linewidth = .16) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#BFD7EA", "#5B91BD", "#174A75"),
    name = paste0("Weeks with mean 2-m air temperature\n\u2264 \u221210 ", degree_c, " (%)"),
    breaks = c(0, 20, 40), limits = c(0, 55), oob = squish,
    na.value = palette[["outside"]],
    guide = guide_colorbar(title.position = "top", barwidth = unit(26, "mm"),
                           barheight = unit(2.2, "mm"))
  ) + map_axes + labs(title = "b  Frequency of cold weeks") + theme_map() +
  theme(legend.position = "bottom")

p_c <- base_map +
  geom_sf(aes(fill = lab_coverage_pct), colour = "white", linewidth = .16) +
  scale_fill_gradientn(
    colours = c("#F7F2E8", "#DFC39F", "#B87848", "#754126"),
    name = "Meteorologically valid weeks\nwith laboratory data (%)",
    breaks = c(25, 50, 75, 100), limits = c(0, 100), oob = squish,
    na.value = palette[["outside"]],
    guide = guide_colorbar(title.position = "top", barwidth = unit(26, "mm"),
                           barheight = unit(2.2, "mm"))
  ) + map_axes + labs(title = "c  Laboratory surveillance coverage") + theme_map() +
  theme(legend.position = "bottom")

year_breaks <- c(as.Date("2017-01-02"), as.Date(paste0(2018:2024, "-01-01")))
x_limits <- as.Date(c("2017-01-02", "2024-12-30"))
make_ts <- function(data, y, label, line_colour, fill_colour,
                    title = NULL, show_x = FALSE) {
  ggplot(data, aes(week_start, .data[[y]])) +
    geom_area(fill = fill_colour, alpha = .55, linewidth = 0, na.rm = TRUE) +
    geom_line(colour = line_colour, linewidth = .32, na.rm = TRUE) +
    annotate("text", x = as.Date("2017-01-09"), y = Inf, label = label,
             hjust = 0, vjust = 1.25, size = 2.2, family = figure_font,
             fontface = "bold", colour = palette[["ink"]]) +
    scale_x_date(breaks = year_breaks, date_labels = "%Y", limits = x_limits,
                 expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(.035, .08))) +
    labs(title = title) + theme_ts(show_x = show_x)
}
p_d1 <- make_ts(weekly, "laboratory_tests", "Influenza tests (n)",
                palette[["blue"]], palette[["blue_fill"]],
                "d  Weekly influenza surveillance indicators", FALSE)
p_d2 <- make_ts(weekly, "laboratory_positivity_pct", "Influenza test positivity (%)",
                palette[["red"]], palette[["red_fill"]], NULL, FALSE)
p_d3 <- make_ts(weekly, "ILI_rate_pct", "ILI consultation rate (18 counties, %)",
                palette[["teal"]], palette[["teal_fill"]], NULL, TRUE)

figure1 <- ((p_a | p_b | p_c) /
              (p_d1 / p_d2 / p_d3 + plot_layout(heights = c(1, 1, 1)))) +
  plot_layout(heights = c(.88, 1.05))

save_multi <- function(plot, base, width_mm, height_mm) {
  width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
  grDevices::png(paste0(base, ".png"), width = width_in, height = height_in,
                 units = "in", res = 600, bg = "white", type = "cairo-png")
  print(plot); grDevices::dev.off()
  grDevices::tiff(paste0(base, ".tiff"), width = width_in, height = height_in,
                  units = "in", res = 600, compression = "lzw", bg = "white")
  print(plot); grDevices::dev.off()
  ggsave(paste0(base, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, bg = "white")
  svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in,
                   bg = "white", fix_text_size = FALSE)
  print(plot); grDevices::dev.off()
}
save_multi(
  figure1,
  file.path(figure_dir, "Figure1_AGU_no_title_with_coordinates_20261008"),
  183, 150
)

writeLines(
  c(
    "Figure 1. Ambient-temperature exposure and influenza surveillance across Qinghai's altitude gradient.",
    "County-level surveillance, 2017\u20132024. a, County-government-seat altitude across the 45 study counties. b, Proportion of meteorologically valid study weeks with weekly mean 2-m air temperature at or below \u221210 \u00B0C. c, Proportion of meteorologically valid study weeks with laboratory surveillance data. d, Weekly numbers of influenza tests, influenza test positivity, and the influenza-like illness (ILI) consultation rate; the ILI series was available for 18 counties."
  ),
  file.path(figure_dir, "Figure_1_AGU_legend_20261008.txt"),
  useBytes = TRUE
)

threshold_map_data <- threshold_map_data |>
  mutate(threshold_label = factor(
    threshold, levels = c(-5, -10, -15),
    labels = paste0(c("-5 ", "-10 ", "-15 "), degree_c)
  ))
supp_threshold <- ggplot(threshold_map_data) +
  geom_sf(aes(fill = cold_week_pct), colour = "white", linewidth = .16) +
  facet_wrap(~threshold_label, nrow = 1) +
  coord_sf(datum = NA) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#BFD7EA", "#5B91BD", "#174A75"),
    name = "Weeks at or below threshold (%)", limits = c(0, 70), oob = squish,
    guide = guide_colorbar(title.position = "top", barwidth = unit(30, "mm"),
                           barheight = unit(2.5, "mm"))
  ) +
  labs(
    title = "Sensitivity of county cold-week prevalence to the air-temperature threshold",
    subtitle = "All meteorologically valid county-weeks, 2017-2024"
  ) +
  theme_void(base_family = figure_font, base_size = 8) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(colour = palette[["muted"]], size = 8)
  )
save_multi(
  supp_threshold,
  file.path(supp_dir, "Result1S_T2m_cold_threshold_sensitivity"),
  183, 72
)

write_csv_utf8(county_context, file.path(source_dir, "result1_county_context_main_threshold.csv"))
write_csv_utf8(threshold_long, file.path(source_dir, "result1_cold_threshold_sensitivity.csv"))
write_csv_utf8(weekly, file.path(source_dir, "result1_weekly_surveillance.csv"))
write_csv_utf8(
  tibble(
    check = c("primary county count", "primary county-weeks", "primary tests",
              "primary positives", "main-threshold counties with cold weeks"),
    value = c(n_distinct(primary$county), nrow(primary), sum(primary$tests),
              sum(primary$positive), sum(county_context$cold_weeks > 0)),
    expected = c(45, 9758, 77930, 11072, 45)
  ),
  file.path(qa_dir, "v4_result1_qa.csv")
)

message("Created v4 T2m Result 1 and threshold-sensitivity figure.")
