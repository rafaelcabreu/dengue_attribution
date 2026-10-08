# 03_predict_ensemble_single_fit_nat.R
#
# Single-fit (no bootstrap) ensemble prediction for Figure3.ipynb's
# "nat_brazil-bs_dengue_ensemble_predictions.csv" input: fits the flagship
# immunity model ONCE on the full historical data (posterior-mean immunity,
# same as 02_fit_model_incidence.R), then predicts onto each of the 525
# HadGEM NAT (counterfactual) ensemble climate files (the same
# data/predict-all/ext-nat-ens*.csv used by
# 03_fit_sensitivity_parallel_full_nat.R's bootstrap pipeline), so the
# resulting ensemble spread reflects ONLY climate uncertainty — no state
# resampling, no immunity-draw resampling.
#
# Output is kept at region x date_first_symptoms x ensemble-member
# granularity (column "file_number", matching what Figure3.ipynb expects and
# then aggregates itself), not collapsed to a single number per city/region.
#
# This is a plain local run (like 02_fit_model_incidence.R /
# 04_extract_model_params_by_region.R) — one fit, then 525 predictions in a
# loop, no SLURM array needed.

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths and settings
# ---------------------------------------------------------------
data_path    <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"
ensemble_dir <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/predict-all"
output_path  <- "../data/nat_brazil-bs_dengue_ensemble_predictions.csv"

tlimit    <- 12
n_bs_cols <- 4

model_name <- "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity"
model_vars <- c(
  'temp_bs_lag11', 'temp_bs_lag12', 'temp_bs_lag13', 'temp_bs_lag14',
  'temp_bs_lag21', 'temp_bs_lag22', 'temp_bs_lag23', 'temp_bs_lag24',
  'temp_bs_lag31', 'temp_bs_lag32', 'temp_bs_lag33', 'temp_bs_lag34',
  'temp_bs_lag41', 'temp_bs_lag42', 'temp_bs_lag43', 'temp_bs_lag44',
  'temp_bs_lag51', 'temp_bs_lag52', 'temp_bs_lag53', 'temp_bs_lag54',
  'total_precipitation_lag1', 'total_precipitation_lag2', 'total_precipitation_lag3',
  'total_precipitation_lag4', 'total_precipitation_lag5',
  'city_residency',
  'immunity_lag1', 'immunity_lag2', 'immunity_lag3',
  'serotype_replacement',
  'urban_area_ha', 'gdp_per_capita', 'births',
  'Pr_0priorinf', 'Pr_1priorinf', 'Pr_2priorinf',
  'year_region', 'month_region'
)

# ---------------------------------------------------------------
# 3. Load and prepare training data
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# ── Municipality-level posterior immunity (posterior mean) ──────────────────────
# Replace the old macro-region-level Pr_0priorinf / Pr_1priorinf with
# municipality-level (city_residency x year) posterior-mean values — see
# 01_model_fit/build_immunity_city_draws.R. Single fixed-effect fit (no
# bootstrap), so the posterior mean is used, same as fit_model.R.
city_health_region   <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)
immunity_mean_region <- read_csv("../data/immunity_posterior_mean_region.csv", show_col_types = FALSE)

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency") %>%
  left_join(immunity_mean_region, by = c("code_health_region", "year"))

# ---------------------------------------------------------------
# 4. Fit B-splines on full training data (knots fixed here)
# ---------------------------------------------------------------
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = n_bs_cols)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = n_bs_cols)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = n_bs_cols)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = n_bs_cols)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = n_bs_cols)

add_spline_columns <- function(df) {
  bs1 <- predict(temp_bs_fitted_lag1, df$mean_2m_air_temp_degree1_lag1)
  bs2 <- predict(temp_bs_fitted_lag2, df$mean_2m_air_temp_degree1_lag2)
  bs3 <- predict(temp_bs_fitted_lag3, df$mean_2m_air_temp_degree1_lag3)
  bs4 <- predict(temp_bs_fitted_lag4, df$mean_2m_air_temp_degree1_lag4)
  bs5 <- predict(temp_bs_fitted_lag5, df$mean_2m_air_temp_degree1_lag5)
  df %>%
    mutate(
      temp_bs_lag11 = bs1[,1], temp_bs_lag12 = bs1[,2],
      temp_bs_lag13 = bs1[,3], temp_bs_lag14 = bs1[,4],
      temp_bs_lag21 = bs2[,1], temp_bs_lag22 = bs2[,2],
      temp_bs_lag23 = bs2[,3], temp_bs_lag24 = bs2[,4],
      temp_bs_lag31 = bs3[,1], temp_bs_lag32 = bs3[,2],
      temp_bs_lag33 = bs3[,3], temp_bs_lag34 = bs3[,4],
      temp_bs_lag41 = bs4[,1], temp_bs_lag42 = bs4[,2],
      temp_bs_lag43 = bs4[,3], temp_bs_lag44 = bs4[,4],
      temp_bs_lag51 = bs5[,1], temp_bs_lag52 = bs5[,2],
      temp_bs_lag53 = bs5[,3], temp_bs_lag54 = bs5[,4]
    )
}

dengue_temp <- add_spline_columns(dengue_temp)
cat("Training data prepared with B-splines.\n\n")

# ---------------------------------------------------------------
# 5. Model formula
# ---------------------------------------------------------------
fixed_effects <- intersect(model_vars, c('month', 'year', 'city_residency',
                                          'month_childs', 'year_childs',
                                          'year_region', 'month_region'))
predictors    <- setdiff(model_vars, fixed_effects)

model_formula <- as.formula(
  paste("n_cases ~", paste(predictors, collapse = " + "),
        "|", paste(fixed_effects, collapse = " + "))
)
cat("Formula:", deparse(model_formula), "\n\n")

# ---------------------------------------------------------------
# 6. Fit model once on full training data
# ---------------------------------------------------------------
fitted_model <- fixest::fepois(
  fml           = model_formula,
  offset        = ~log_pop_offset,
  data          = dengue_temp,
  combine.quick = FALSE
)
cat("Model fitted.\n\n")

# ---------------------------------------------------------------
# 7. Predict across all 525 ensemble members, keeping region x
#    date_first_symptoms x ensemble-member resolution
# ---------------------------------------------------------------
file_list   <- paste0(ensemble_dir, "/ext-nat-ens", sprintf("%03d", 0:524), ".csv")
n_ensembles <- length(file_list)
results_list <- vector("list", n_ensembles)

for (ensemble_idx in seq_along(file_list)) {
  cat("  Ensemble member", ensemble_idx, "of", n_ensembles, "\n")

  results_list[[ensemble_idx]] <- tryCatch({
    new_data <- read_csv(file_list[[ensemble_idx]], show_col_types = FALSE)
    new_data$log_pop_offset <- log(new_data$population / 100000)
    new_data <- add_spline_columns(new_data)

    # Same posterior-mean immunity as the training data — climate is the
    # only thing that differs across ensemble members.
    new_data <- new_data %>%
      select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
      left_join(city_health_region, by = "city_residency") %>%
      left_join(immunity_mean_region, by = c("code_health_region", "year"))

    new_data <- new_data %>%
      filter(
        mean_2m_air_temp_degree1_lag1 > tlimit,
        mean_2m_air_temp_degree1_lag2 > tlimit,
        mean_2m_air_temp_degree1_lag3 > tlimit,
        mean_2m_air_temp_degree1_lag4 > tlimit,
        mean_2m_air_temp_degree1_lag5 > tlimit
      )

    pred_cases <- predict(fitted_model, newdata = new_data)

    new_data %>%
      mutate(pred_cases = pred_cases) %>%
      group_by(region, date_first_symptoms) %>%
      summarise(
        total_population = sum(population,  na.rm = TRUE),
        total_pred_cases = sum(pred_cases,  na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        predicted_incidence = total_pred_cases / total_population * 100000,
        model_type          = model_name,
        file_number         = ensemble_idx - 1   # 0-indexed, matches ext-ens000..524 naming
      ) %>%
      select(model_type, region, date_first_symptoms, total_population,
             total_pred_cases, predicted_incidence, file_number)

  }, error = function(e) {
    cat("    ERROR on ensemble member", ensemble_idx, ":", e$message, "\n")
    NULL
  })
}

final_results <- bind_rows(results_list)
cat("\nSuccessfully processed", n_distinct(final_results$file_number), "of", n_ensembles, "ensemble members.\n")

write_csv(final_results, output_path)
cat("Saved results to:", output_path, "\n")
cat("Done!\n")
