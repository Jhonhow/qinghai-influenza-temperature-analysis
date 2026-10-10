suppressPackageStartupMessages({
  library(dplyr)
  library(geepack)
})

analysis_id <- "confirmatory_anomaly_shape_20260831"
checkpoint_dir <- file.path(
  "model_outputs", "reanalysis_v4_t2m_20260825", analysis_id
)
source_dir <- file.path(
  "result0715", "source_data", paste0("v4_", analysis_id)
)
bundle <- readRDS(file.path(checkpoint_dir, "confirmatory_anomaly_bundle.rds"))
d <- bundle$confirm_dat |>
  arrange(county, week_index) |>
  mutate(
    county = factor(county),
    week_id = factor(week_id),
    positivity = positive / tests
  )

spline_names <- bundle$spline_names
spline_int_names <- bundle$spline_int_names
clim_names <- grep("^clim_", names(d), value = TRUE)

make_formula <- function(include_interaction) {
  terms <- c(
    spline_names,
    if (include_interaction) spline_int_names else character(),
    clim_names, "county", "week_id"
  )
  as.formula(paste("positivity ~", paste(terms, collapse = " + ")))
}

fit_or_load <- function(path, formula) {
  if (file.exists(path)) return(readRDS(path))
  fit <- geeglm(
    formula,
    family = binomial(link = "logit"),
    data = d,
    weights = tests,
    id = county,
    waves = week_index,
    corstr = "ar1",
    std.err = "san.se"
  )
  saveRDS(fit, path)
  fit
}

fit_reduced <- fit_or_load(
  file.path(checkpoint_dir, "twoway_fixed_effects_gee_no_int.rds"),
  make_formula(FALSE)
)
fit_full <- fit_or_load(
  file.path(checkpoint_dir, "twoway_fixed_effects_gee_alt.rds"),
  make_formula(TRUE)
)

beta <- coef(fit_full)
variance <- fit_full$geese$vbeta
dimnames(variance) <- list(names(beta), names(beta))
keep <- intersect(spline_int_names, names(beta))
b <- beta[keep]
v <- variance[keep, keep, drop = FALSE]
statistic <- as.numeric(t(b) %*% solve(v, b))
p_value <- pchisq(statistic, df = length(keep), lower.tail = FALSE)

result <- tibble(
  analysis = "County and province-week fixed effects; GEE-AR1",
  statistic = statistic,
  df = length(keep),
  p_value = p_value,
  alpha_reduced = fit_reduced$geese$alpha,
  alpha_full = fit_full$geese$alpha,
  clusters = nlevels(d$county),
  observations = nrow(d)
)

write.csv(
  result,
  file.path(source_dir, "confirmatory_twoway_fixed_effects_gee.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

message("Completed two-way fixed-effects GEE sensitivity analysis.")
