# fit_all_models_once.R
#
# Fits each model once on the full dataset (no bootstrap).
# Outputs three CSV files:
#   1. all_models_coefficients.csv       — one row per coefficient per model,
#      including estimate, std error, z-value, p-value, and significance stars.
#   2. all_models_fit_stats.csv          — one row per model with AIC, R², pseudo-R².
#   3. all_models_regional_incidence.csv — observed and fitted monthly incidence
#      (cases per 100k) aggregated by region, year, and month, for every model.

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
data_path  <- "data/model_input_brazil_immunity_city_with_priorinf.csv"
output_dir <- "data/statistics/full_model_fits"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definitions (identical to bootstrap script)
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
    'Pr_1priorinf',
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
    'Pr_0priorinf','Pr_1priorinf'
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
    'Pr_0priorinf','Pr_1priorinf',
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
    'city_residency','Pr_0priorinf','Pr_1priorinf'
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
    'Pr_0priorinf','Pr_1priorinf'
  ),
  'Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'month_region',
    'city_residency',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births',
    'Pr_0priorinf','Pr_1priorinf'
  )
)

# ---------------------------------------------------------------
# 4. Load and prepare data
# ---------------------------------------------------------------
cat("Loading data...\n")
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# Fit B-splines on full data (consistent knots)
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

cat("Data prepared.\n\n")

# ---------------------------------------------------------------
# 5. Helper: build formula (identical logic to bootstrap script)
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

# ---------------------------------------------------------------
# 6. Helper: significance stars from p-value
# ---------------------------------------------------------------
sig_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE       ~ ""
  )
}

# ---------------------------------------------------------------
# 7. Fit all models and collect results
# ---------------------------------------------------------------
all_coefs      <- list()
all_fit_stats  <- list()
all_incidence  <- list()

for (model_idx in seq_along(models_with_temp)) {

  model_name <- names(models_with_temp)[model_idx]
  model_vars <- models_with_temp[[model_idx]]

  cat(sprintf("[%2d/%2d] Fitting: %s\n",
              model_idx, length(models_with_temp), model_name))

  model_formula <- build_formula(model_vars)
  print(model_formula)

  fit_result <- tryCatch({
    fixest::fepois(
      model_formula,
      offset   = ~log_pop_offset,
      data     = dengue_temp,
      nthreads = 1
    )
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    NULL
  })

  print(fit_result)
  
  if (is.null(fit_result)) next

  # --- Model-level fit statistics ---
  fit_stats <- data.frame(
    model_index      = model_idx,
    model_name       = model_name,
    aic              = AIC(fit_result),
    bic              = BIC(fit_result),
    r_squared        = r2(fit_result, type = "cor2"),
    pseudo_r_squared = r2(fit_result, type = "pr2"),
    n_obs            = fit_result$nobs,
    n_fixed_effects  = length(fixef(fit_result)),
    stringsAsFactors = FALSE
  )
  all_fit_stats[[model_idx]] <- fit_stats

  # --- Coefficient table ---
  # coeftable() returns estimate, std error, t/z stat, p-value for non-FE terms
  ct <- coeftable(fit_result)

  if (!is.null(ct) && nrow(ct) > 0) {
    coef_df <- as.data.frame(ct)
    coef_df$term        <- rownames(ct)
    rownames(coef_df)   <- NULL

    # Standardise column names (fixest uses "Std. Error", "z value", "Pr(>|z|)")
    colnames(coef_df) <- make.names(colnames(coef_df))  # safe names
    # Rename to consistent names regardless of fixest version
    col_map <- c(
      "Estimate"      = "estimate",
      "Std..Error"    = "std_error",
      "z.value"       = "z_value",
      "t.value"       = "z_value",      # alias
      "Pr...z.."      = "p_value",
      "Pr...t.."      = "p_value"       # alias
    )
    for (old in names(col_map)) {
      if (old %in% colnames(coef_df)) {
        colnames(coef_df)[colnames(coef_df) == old] <- col_map[[old]]
      }
    }

    coef_df <- coef_df %>%
      mutate(
        model_index  = model_idx,
        model_name   = model_name,
        significance = sig_stars(p_value),
        irr          = exp(estimate)   # incidence rate ratio for Poisson
      ) %>%
      select(model_index, model_name, term,
             estimate, std_error, z_value, p_value, significance, irr)

    all_coefs[[model_idx]] <- coef_df
  }

  # --- Regional monthly incidence: observed vs fitted ---
  # Attach fitted counts to the rows used by this model, then aggregate
  # up to region × year × month, summing cases and population separately
  # before converting to incidence per 100k.
  #
  # We use fitted() which returns values only for the rows that were
  # actually used in estimation (complete cases for this model's variables).

  # predict(..., sample = "original") always returns a vector of length
  # nrow(dengue_temp), padding dropped rows (NAs, zero-FE singletons, etc.)
  # with NA — so no manual index reconstruction is needed.
  fitted_full <- tryCatch(
    predict(fit_result, sample = "original"),
    error = function(e) NULL
  )

  if (!is.null(fitted_full)) {

    incidence_df <- dengue_temp %>%
      mutate(.fitted_cases = as.numeric(fitted_full)) %>%
      filter(!is.na(.fitted_cases)) %>%
      { if (all(c("region", "year", "month", "n_cases", "population") %in% names(.))) . else NULL }

    if (!is.null(incidence_df) && nrow(incidence_df) > 0) {
      incidence_agg <- incidence_df %>%
        group_by(region, year, month) %>%
        summarise(
          observed_cases   = sum(n_cases,         na.rm = TRUE),
          fitted_cases     = sum(.fitted_cases,   na.rm = TRUE),
          total_population = sum(population,      na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          observed_incidence_per100k = observed_cases / total_population * 100000,
          fitted_incidence_per100k   = fitted_cases   / total_population * 100000,
          model_index                = model_idx,
          model_name                 = model_name
        ) %>%
        select(model_index, model_name, region, year, month,
               observed_cases, fitted_cases, total_population,
               observed_incidence_per100k, fitted_incidence_per100k)

      all_incidence[[model_idx]] <- incidence_agg
    }
  }

  cat(sprintf("  AIC=%.1f  BIC=%.1f  R²=%.4f  Pseudo-R²=%.4f  n=%d\n",
              fit_stats$aic, fit_stats$bic, fit_stats$r_squared, fit_stats$pseudo_r_squared,
              fit_stats$n_obs))
}

# ---------------------------------------------------------------
# 8. Combine and save
# ---------------------------------------------------------------
coef_table        <- bind_rows(all_coefs)
fit_stats_table   <- bind_rows(all_fit_stats)
incidence_table   <- bind_rows(all_incidence)

coef_path         <- file.path(output_dir, "all_models_coefficients.csv")
fit_stats_path    <- file.path(output_dir, "all_models_fit_stats.csv")
incidence_path    <- file.path(output_dir, "all_models_regional_incidence.csv")

write_csv(coef_table,      coef_path)
write_csv(fit_stats_table, fit_stats_path)
write_csv(incidence_table, incidence_path)

cat("\n=============================================\n")
cat("Saved coefficient table to:      ", coef_path,      "\n")
cat("Saved fit statistics table to:   ", fit_stats_path, "\n")
cat("Saved regional incidence table to:", incidence_path, "\n")

# ---------------------------------------------------------------
# 9. Print a quick summary to console
# ---------------------------------------------------------------
cat("\n--- Model fit statistics summary ---\n")
fit_stats_table %>%
  select(model_index, model_name, aic, r_squared, pseudo_r_squared, n_obs) %>%
  mutate(across(c(aic, r_squared, pseudo_r_squared), ~round(., 4))) %>%
  print(n = Inf)

cat("\nDone!\n")