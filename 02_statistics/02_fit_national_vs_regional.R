# fit_national_vs_regional.R
#
# Generates the two CSVs consumed by FigureS10.ipynb — the comparison of
# each model's fit when estimated NATIONALLY (one set of coefficients for
# the whole country) versus REGIONALLY (a separate fit per region, on that
# region's data only):
#
#   1. cases_by_region_month_all_models.csv
#      One national fit per model (identical to fit_model.R's fit), with
#      fitted vs. observed cases aggregated to region x month.
#   2. cases_by_region_month_all_models_by_region.csv
#      A separate fit per model PER REGION, refit from scratch on that
#      region's rows only, then aggregated to month within that region.
#      'Year + Month' is skipped here: it has no region-varying terms, so a
#      region-only refit of it is not a meaningful comparison (this matches
#      the model list already present in the existing static file).
#
# Both files share the same column schema:
#   model_name, region, date_first_symptoms, total_population,
#   total_actual_cases, total_pred_cases, actual_incidence, predicted_incidence
#
# This is a plain local run (no SLURM array) — each model gets 1 national
# fit + up to 5 regional fits, so ~22 x 6 = ~130 fepois() calls total.
# Regional subsets are ~5x smaller than the national dataset, so those fits
# are individually faster, but expect this to take a while end to end.

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
data_path        <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"
national_out_path <- "../data/cases_by_region_month_all_models_v2.csv"
regional_out_path <- "../data/cases_by_region_month_all_models_by_region_v2.csv"

# ---------------------------------------------------------------
# 3. Model definitions (identical to fit_model.R)
# ---------------------------------------------------------------
models_with_temp <- list(
  'Year + Month' = c(
    'month',
    'year',
    'city_residency'
  ),
  'Year|region + Month|region' = c(
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha',
    'gdp_per_capita',
    'births',
    'Pr_0priorinf',
    'Pr_1priorinf', 'Pr_2priorinf',
    'year_region',
    'month_region'
  ),
  'Climate(lag 1-5)' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Month|region' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'city_residency'
  ),
  'Year + Month + P(Climate)' = c(
    'mean_2m_air_temp_degree1_lag1','mean_2m_air_temp_degree1_lag2',
    'mean_2m_air_temp_degree1_lag3',
    'mean_2m_air_temp_degree2_lag1','mean_2m_air_temp_degree2_lag2',
    'mean_2m_air_temp_degree2_lag3',
    'mean_2m_air_temp_degree3_lag1','mean_2m_air_temp_degree3_lag2',
    'mean_2m_air_temp_degree3_lag3',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3',
    'city_residency',
    'month',
    'year'
  ),
  'Climate(lag 1) + Year|region + Month|region' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'total_precipitation_lag1',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-3) + Year|region + Month|region' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region + Month|region' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + Socio' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region','year_region','city_residency',
    'urban_area_ha','gdp_per_capita','births'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region','year_region','city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + SeroRepla' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region','year_region','city_residency',
    'serotype_replacement'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + Immunity' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region','year_region','city_residency',
    'Pr_0priorinf','Pr_1priorinf', 'Pr_2priorinf'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births',
    'Pr_0priorinf','Pr_1priorinf', 'Pr_2priorinf',
    'year_region','month_region'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births',
    'year_region','month_region'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births'
  ),
  'Climate(lag 1-5) + Socio' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency','urban_area_ha','gdp_per_capita','births'
  ),
  'Climate(lag 1-5) + PriorCases' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency','immunity_lag1','immunity_lag2','immunity_lag3'
  ),
  'Climate(lag 1-5) + SeroRepla' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency','serotype_replacement'
  ),
  'Climate(lag 1-5) + Immunity' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency','Pr_0priorinf','Pr_1priorinf', 'Pr_2priorinf'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births',
    'Pr_0priorinf','Pr_1priorinf', 'Pr_2priorinf'
  ),
  'Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'month_region',
    'city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births',
    'Pr_0priorinf','Pr_1priorinf', 'Pr_2priorinf'
  )
)

# Models with no region-varying term — a region-only refit of these is not a
# meaningful "regional fit" comparison, so they're skipped in the by-region
# output (matches the model list already present in the existing static file).
skip_for_regional <- c('Year + Month')

# ---------------------------------------------------------------
# 4. Load and prepare data
# ---------------------------------------------------------------
cat("Loading data...\n")
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# ── Municipality-level posterior immunity (posterior mean) ──────────────────────
# Replace the old macro-region-level Pr_0priorinf / Pr_1priorinf with
# municipality-level (city_residency x year) posterior-mean values — see
# 01_model_fit/build_immunity_city_draws.R. This is a point-estimate (no
# bootstrap) script, so the posterior mean is used, same as fit_model.R.
city_health_region   <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)
immunity_mean_region <- read_csv("../data/immunity_posterior_mean_region.csv", show_col_types = FALSE)

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency") %>%
  left_join(immunity_mean_region, by = c("code_health_region", "year"))

temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

dengue_temp <- dengue_temp %>%
  mutate(
    temp_bs_lag11 = temp_bs_fitted_lag1[,1], temp_bs_lag12 = temp_bs_fitted_lag1[,2],
    temp_bs_lag13 = temp_bs_fitted_lag1[,3], temp_bs_lag14 = temp_bs_fitted_lag1[,4],
    temp_bs_lag21 = temp_bs_fitted_lag2[,1], temp_bs_lag22 = temp_bs_fitted_lag2[,2],
    temp_bs_lag23 = temp_bs_fitted_lag2[,3], temp_bs_lag24 = temp_bs_fitted_lag2[,4],
    temp_bs_lag31 = temp_bs_fitted_lag3[,1], temp_bs_lag32 = temp_bs_fitted_lag3[,2],
    temp_bs_lag33 = temp_bs_fitted_lag3[,3], temp_bs_lag34 = temp_bs_fitted_lag3[,4],
    temp_bs_lag41 = temp_bs_fitted_lag4[,1], temp_bs_lag42 = temp_bs_fitted_lag4[,2],
    temp_bs_lag43 = temp_bs_fitted_lag4[,3], temp_bs_lag44 = temp_bs_fitted_lag4[,4],
    temp_bs_lag51 = temp_bs_fitted_lag5[,1], temp_bs_lag52 = temp_bs_fitted_lag5[,2],
    temp_bs_lag53 = temp_bs_fitted_lag5[,3], temp_bs_lag54 = temp_bs_fitted_lag5[,4]
  )

regions_all <- sort(unique(dengue_temp$region))
cat("Data prepared. Regions:", paste(regions_all, collapse = ", "), "\n\n")

# ---------------------------------------------------------------
# 5. Helpers
# ---------------------------------------------------------------
build_formula <- function(vars) {
  fixed_effects <- intersect(vars, c("city_residency", "month", "year",
                                     "month_childs", "year_childs",
                                     "year_region", "month_region"))
  predictors    <- setdiff(vars, fixed_effects)
  if (length(predictors) > 0 && length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "),
                         "|", paste(fixed_effects, collapse = " + "))
  } else if (length(predictors) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "))
  } else if (length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~ 1 |", paste(fixed_effects, collapse = " + "))
  } else {
    formula_str <- "n_cases ~ 1"
  }
  as.formula(formula_str)
}

# Fits `model_formula` on `data`, and aggregates fitted vs. observed cases to
# region x date_first_symptoms. Returns NULL on any failure (bad fit, no
# predictions) so the caller can skip and move on.
fit_and_aggregate <- function(data, model_formula) {
  fit_result <- tryCatch({
    fixest::fepois(model_formula, offset = ~log_pop_offset, data = data, nthreads = 1)
  }, error = function(e) {
    cat("  ERROR fitting:", e$message, "\n")
    NULL
  })
  if (is.null(fit_result)) return(NULL)

  fitted_full <- tryCatch(predict(fit_result, sample = "original"), error = function(e) NULL)
  if (is.null(fitted_full)) return(NULL)

  data %>%
    mutate(.fitted_cases = as.numeric(fitted_full)) %>%
    filter(!is.na(.fitted_cases)) %>%
    group_by(region, date_first_symptoms) %>%
    summarise(
      total_population    = sum(population, na.rm = TRUE),
      total_actual_cases  = sum(n_cases,    na.rm = TRUE),
      total_pred_cases    = sum(.fitted_cases, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      actual_incidence    = total_actual_cases / total_population * 100000,
      predicted_incidence = total_pred_cases   / total_population * 100000
    )
}

# ---------------------------------------------------------------
# 6. Fit every model nationally, and (except skip_for_regional) per region
# ---------------------------------------------------------------
national_results <- list()
regional_results  <- list()

model_names <- names(models_with_temp)

for (model_idx in seq_along(models_with_temp)) {
  model_name <- model_names[model_idx]
  model_vars <- models_with_temp[[model_idx]]
  model_formula <- build_formula(model_vars)

  cat(sprintf("[%2d/%2d] %s\n", model_idx, length(models_with_temp), model_name))

  # --- National fit ---
  cat("  Fitting nationally...\n")
  nat_agg <- fit_and_aggregate(dengue_temp, model_formula)
  if (!is.null(nat_agg)) {
    nat_agg$model_name <- model_name
    national_results[[length(national_results) + 1]] <- nat_agg
  }

  # --- Regional-only fits ---
  if (model_name %in% skip_for_regional) {
    cat("  Skipping regional fits (no region-varying term).\n\n")
    next
  }

  for (region in regions_all) {
    cat("  Fitting region:", region, "...\n")
    region_data <- dengue_temp %>% filter(region == !!region)

    reg_agg <- fit_and_aggregate(region_data, model_formula)
    if (!is.null(reg_agg)) {
      reg_agg$model_name <- model_name
      regional_results[[length(regional_results) + 1]] <- reg_agg
    }
  }
  cat("\n")
}

# ---------------------------------------------------------------
# 7. Combine and save
# ---------------------------------------------------------------
col_order <- c("model_name", "region", "date_first_symptoms", "total_population",
               "total_actual_cases", "total_pred_cases", "actual_incidence",
               "predicted_incidence")

national_table <- bind_rows(national_results) %>% select(all_of(col_order))
regional_table <- bind_rows(regional_results) %>% select(all_of(col_order))

write_csv(national_table, national_out_path)
write_csv(regional_table, regional_out_path)

cat("\n=============================================\n")
cat("Saved national fit table to:", national_out_path, "(", nrow(national_table), "rows )\n")
cat("Saved regional fit table to:", regional_out_path, "(", nrow(regional_table), "rows )\n")
cat("Done!\n")
