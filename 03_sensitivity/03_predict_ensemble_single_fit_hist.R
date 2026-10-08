# 03_predict_ensemble_single_fit_hist.R
#
# Single-fit (no bootstrap) ensemble prediction for Figure3.ipynb's
# "hist_brazil-bs_dengue_ensemble_predictions.csv" input: fits the flagship
# immunity model ONCE on the full historical data (posterior-mean immunity,
# same as 02_fit_model_incidence.R), then predicts onto a HISTORICAL BASELINE
# climate ensemble — unlike act/nat (one season, 525 members), this spans
# MULTIPLE past years (Figure3 filters to year < 2014), so the loop below is
# nested over year x ensemble-member.
#
# *** PLACEHOLDER ***
# hist_ensemble_dir / hist_years / the file-naming pattern in the loop below
# are all placeholders — fill in the real path and naming convention for the
# historical ensemble climate files once known, then remove this notice.
#
# Output is kept at region x date_first_symptoms x year x ensemble-member
# granularity (column "file_number"), matching what Figure3.ipynb expects
# (it aggregates further itself via groupby(['region', 'year', 'file_number'])).
#
# This is a plain local run (like 02_fit_model_incidence.R /
# 04_extract_model_params_by_region.R) — one fit, then many predictions in a
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
data_path <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"

# *** PLACEHOLDER — fill in the real historical ensemble directory, the years
# it covers, and the exact file-naming pattern (see the loop in section 7). ***
hist_ensemble_dir <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/predict-hist"   # PLACEHOLDER
hist_years        <- 2005:2013                 # PLACEHOLDER — matches Figure3.ipynb's `year < 2014` filter, confirm exact range
n_hist_ensembles  <- 525                        # PLACEHOLDER — confirm member count for the historical ensemble

output_path <- "../data/hist_brazil-bs_dengue_ensemble_predictions.csv"

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
# 7. Predict across every historical year x ensemble member, keeping
#    region x date_first_symptoms x year x ensemble-member resolution
# ---------------------------------------------------------------
results_list <- list()

for (hist_year in hist_years) {
  for (ensemble_idx in seq_len(n_hist_ensembles)) {
    cat("  Year", hist_year, "- ensemble member", ensemble_idx, "of", n_hist_ensembles, "\n")

    # *** PLACEHOLDER file-naming pattern — replace once known ***
    file_path <- sprintf("%s/ext-hist-%d-ens%03d.csv",
                         hist_ensemble_dir, hist_year, ensemble_idx - 1)

    result <- tryCatch({
      new_data <- read_csv(file_path, show_col_types = FALSE)
      new_data$log_pop_offset <- log(new_data$population / 100000)
      new_data <- add_spline_columns(new_data)

      # Same posterior-mean immunity as the training data — climate is the
      # only thing that differs across ensemble members/years.
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
          year                = hist_year,
          file_number         = ensemble_idx - 1   # 0-indexed
        ) %>%
        select(model_type, region, date_first_symptoms, year, total_population,
               total_pred_cases, predicted_incidence, file_number)

    }, error = function(e) {
      cat("    ERROR on year", hist_year, "ensemble member", ensemble_idx, ":", e$message, "\n")
      NULL
    })

    if (!is.null(result)) results_list[[length(results_list) + 1]] <- result
  }
}

final_results <- bind_rows(results_list)
cat("\nSuccessfully processed", nrow(distinct(final_results, year, file_number)),
    "of", length(hist_years) * n_hist_ensembles, "year x ensemble combinations.\n")

write_csv(final_results, output_path)
cat("Saved results to:", output_path, "\n")
cat("Done!\n")
