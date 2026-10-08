# fit_city_fe_models.R
#
# Fits C15-c(PC,SR,S) and C15-c(PC,SR,S,IM) — Climate(lag 1-5) + PriorCases +
# SeroRepla + Socio [+ Immunity] — using observed temperature, with
# city_residency as the ONLY fixed effect (no year_region/month_region).
# Data prep (spline basis, immunity join) is identical to fit_model.R so the
# two fits are directly comparable to the rest of the pipeline's outputs.
#
# Outputs three CSV files under data/statistics/city_fe_model_fits/:
#   1. coefficients.csv        — one row per coefficient per model.
#   2. fit_stats.csv           — one row per model with AIC, R², pseudo-R².
#   3. city_predicted_incidence.csv — observed and fitted monthly incidence
#      (cases per 100k) per city_residency, for both models.

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
data_path  <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"
output_dir <- "../data/statistics/city_fe_model_fits"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definitions
#    C15-c(PC,SR,S,IM) is identical to the "Climate(lag 1-5) + PriorCases +
#    SeroRepla + Socio + Immunity" entry in fit_model.R (city_residency is
#    already the only FE there). C15-c(PC,SR,S) is the same minus the
#    Immunity (Pr_0/1/2priorinf) terms — not present in fit_model.R, added
#    here.
# ---------------------------------------------------------------
models_with_temp <- list(
  'C15-c(PC,SR,S,IM)' = c(
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

  'C15-c(PC,SR,S)' = c(
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

  'C15-Yr-Mr-c(PC,SR,S,IM)' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'year_region',
    'month_region',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births',
    'Pr_0priorinf','Pr_1priorinf'
  ),

  'C15-Yr-Mr-c(PC,SR,S)' = c(
    'temp_bs_lag11','temp_bs_lag12','temp_bs_lag13','temp_bs_lag14',
    'temp_bs_lag21','temp_bs_lag22','temp_bs_lag23','temp_bs_lag24',
    'temp_bs_lag31','temp_bs_lag32','temp_bs_lag33','temp_bs_lag34',
    'temp_bs_lag41','temp_bs_lag42','temp_bs_lag43','temp_bs_lag44',
    'temp_bs_lag51','temp_bs_lag52','temp_bs_lag53','temp_bs_lag54',
    'total_precipitation_lag1','total_precipitation_lag2',
    'total_precipitation_lag3','total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'year_region',
    'month_region',
    'immunity_lag1','immunity_lag2','immunity_lag3',
    'serotype_replacement',
    'urban_area_ha','gdp_per_capita','births'
  )
)

# ---------------------------------------------------------------
# 4. Load and prepare data (identical to fit_model.R)
# ---------------------------------------------------------------
cat("Loading data...\n")
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

city_health_region   <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)
immunity_mean_region <- read_csv("../data/immunity_posterior_mean_region.csv", show_col_types = FALSE)

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency") %>%
  left_join(immunity_mean_region, by = c("code_health_region", "year"))

# Fit B-splines on full data (consistent knots), using OBSERVED temperature —
# same basis as every other model in this pipeline.
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
# 5. Helper: build formula (identical logic to fit_model.R)
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
# 6. Fit both models and collect results
# ---------------------------------------------------------------
all_coefs     <- list()
all_fit_stats <- list()
all_incidence <- list()

for (model_idx in seq_along(models_with_temp)) {

  model_name <- names(models_with_temp)[model_idx]
  model_vars <- models_with_temp[[model_idx]]

  cat(sprintf("[%d/%d] Fitting: %s\n", model_idx, length(models_with_temp), model_name))

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
  ct <- coeftable(fit_result)

  if (!is.null(ct) && nrow(ct) > 0) {
    coef_df <- as.data.frame(ct)
    coef_df$term      <- rownames(ct)
    rownames(coef_df) <- NULL

    colnames(coef_df) <- make.names(colnames(coef_df))
    col_map <- c(
      "Estimate"   = "estimate",
      "Std..Error" = "std_error",
      "z.value"    = "z_value",
      "t.value"    = "z_value",
      "Pr...z.."   = "p_value",
      "Pr...t.."   = "p_value"
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
        irr          = exp(estimate)
      ) %>%
      select(model_index, model_name, term,
             estimate, std_error, z_value, p_value, significance, irr)

    all_coefs[[model_idx]] <- coef_df
  }

  # --- Per-city monthly incidence: observed vs fitted ---
  # predict(..., sample = "original") returns length nrow(dengue_temp),
  # padding rows dropped at estimation time (NAs, singleton FEs, etc.) with NA.
  fitted_full <- tryCatch(
    predict(fit_result, sample = "original"),
    error = function(e) NULL
  )

  if (!is.null(fitted_full)) {

    incidence_df <- dengue_temp %>%
      mutate(.fitted_cases = as.numeric(fitted_full)) %>%
      filter(!is.na(.fitted_cases))

    incidence_agg <- incidence_df %>%
      group_by(city_residency, region, date_first_symptoms, year, month) %>%
      summarise(
        observed_cases   = sum(n_cases,       na.rm = TRUE),
        fitted_cases     = sum(.fitted_cases, na.rm = TRUE),
        population       = sum(population,    na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        observed_incidence_per100k = observed_cases / population * 100000,
        fitted_incidence_per100k   = fitted_cases   / population * 100000,
        model_index                = model_idx,
        model_name                 = model_name
      ) %>%
      select(model_index, model_name, city_residency, region, date_first_symptoms,
             year, month, observed_cases, fitted_cases, population,
             observed_incidence_per100k, fitted_incidence_per100k)

    all_incidence[[model_idx]] <- incidence_agg
  }

  cat(sprintf("  AIC=%.1f  BIC=%.1f  R²=%.4f  Pseudo-R²=%.4f  n=%d\n",
              fit_stats$aic, fit_stats$bic, fit_stats$r_squared, fit_stats$pseudo_r_squared,
              fit_stats$n_obs))
}

# ---------------------------------------------------------------
# 7. Combine and save
# ---------------------------------------------------------------
coef_table      <- bind_rows(all_coefs)
fit_stats_table <- bind_rows(all_fit_stats)
incidence_table <- bind_rows(all_incidence)

coef_path      <- file.path(output_dir, "coefficients.csv")
fit_stats_path <- file.path(output_dir, "fit_stats.csv")
incidence_path <- file.path(output_dir, "city_predicted_incidence.csv")

write_csv(coef_table,      coef_path)
write_csv(fit_stats_table, fit_stats_path)
write_csv(incidence_table, incidence_path)

cat("\n=============================================\n")
cat("Saved coefficient table to:        ", coef_path,      "\n")
cat("Saved fit statistics table to:     ", fit_stats_path, "\n")
cat("Saved per-city incidence table to: ", incidence_path, "\n")

cat("\n--- Model fit statistics summary ---\n")
fit_stats_table %>%
  select(model_index, model_name, aic, r_squared, pseudo_r_squared, n_obs) %>%
  mutate(across(c(aic, r_squared, pseudo_r_squared), ~round(., 4))) %>%
  print(n = Inf)

cat("\nDone!\n")
