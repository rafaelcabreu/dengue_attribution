# 03_fit_sensitivity_parallel_full_nat_by_region.R
#
# BY-REGION variant of 03_fit_sensitivity_parallel_full_nat.R: instead of one
# national model applied to every region's ensemble rows, this fits a
# SEPARATE model per region — using only that region's rows of the
# state-resampled bootstrap data — and predicts each region's ensemble rows
# with that region's own coefficients. Mirrors the national-vs-regional
# pattern already used elsewhere in this pipeline (02_fit_national_vs_regional.R,
# 04_extract_model_params_by_region.R): same model formula, refit separately
# per region rather than adding region as a covariate/interaction.
#
# The region loop happens INSIDE each task (not as an extra SLURM array
# dimension), so the array is the same size as the national script's: one
# bootstrap iteration x one model per task. Each task now fits ~5 region
# models instead of 1, so expect longer per-task runtime.
#
# The SLURM array index encodes BOTH bootstrap_idx and model_idx.
# With n_bootstrap=100 and n_models=18 (this script currently has one fewer
# model than the _act variant — 'Climate(lag 1-5) + PriorCases + SeroRepla +
# Socio' isn't in this list; kept as-is, not re-added), total tasks = 1800:
#   task 1..100    -> model 1, bootstrap 1..100
#   task 101..200  -> model 2, bootstrap 1..100
#   ...
#
# This mapping is handled automatically below.
#
# Submit with:
#   sbatch 03_submit_bootstrap_nat_by_region.sbatch
#
# After all tasks complete, combine outputs with:
#   Rscript 03_aggregate_results_nat_by_region.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript model_slurm_parallel.R <task_id>")

task_id     <- as.integer(args[1])
n_bootstrap <- 100
n_models    <- 18

# Derive model_idx and bootstrap_idx from flat task_id
model_idx     <- ceiling(task_id / n_bootstrap)
bootstrap_idx <- task_id - (model_idx - 1) * n_bootstrap

cat(sprintf("Task %d -> model_idx=%d, bootstrap_idx=%d\n",
            task_id, model_idx, bootstrap_idx))

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths and settings
# ---------------------------------------------------------------
data_path      <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"
bootstrap_path <- "../data/bootstrap_state_samples.csv"
ensemble_dir   <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/predict-all"
output_dir     <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/sensitivity-nat/bootstrap_chunks_by_region"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tlimit    <- 12
n_bs_cols <- 4

# ---------------------------------------------------------------
# 3. Model definitions
# ---------------------------------------------------------------
models_with_temp <- list(
  'Climate(lag 1-5)' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Month|region' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'city_residency'
  ),
  'Year + Month + P(Climate)' = c(
    'mean_2m_air_temp_degree1_lag1',
    'mean_2m_air_temp_degree1_lag2',
    'mean_2m_air_temp_degree1_lag3',
    'mean_2m_air_temp_degree2_lag1',
    'mean_2m_air_temp_degree2_lag2',
    'mean_2m_air_temp_degree2_lag3',
    'mean_2m_air_temp_degree3_lag1',
    'mean_2m_air_temp_degree3_lag2',
    'mean_2m_air_temp_degree3_lag3',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'city_residency',
    'month',
    'year'
  ),
  'Climate(lag 1) + Year|region + Month|region' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'total_precipitation_lag1',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-3) + Year|region + Month|region' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region + Month|region' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + Socio' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'year_region',
    'city_residency',
    'urban_area_ha',
    'gdp_per_capita',
    'births'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'year_region',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + SeroRepla' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'year_region',
    'city_residency',
    'serotype_replacement'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + Immunity' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'month_region',
    'year_region',
    'city_residency',
    'Pr_0priorinf',
    'Pr_1priorinf', 'Pr_2priorinf'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
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
  'Climate(lag 1-5) + Socio' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'urban_area_ha',
    'gdp_per_capita',
    'births'
  ),
  'Climate(lag 1-5) + PriorCases' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3'
  ),
  'Climate(lag 1-5) + SeroRepla' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'serotype_replacement'
  ),
  'Climate(lag 1-5) + Immunity' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'Pr_0priorinf',
    'Pr_1priorinf', 'Pr_2priorinf'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha',
    'gdp_per_capita',
    'births',
    'Pr_0priorinf',
    'Pr_1priorinf', 'Pr_2priorinf'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',
    'temp_bs_lag21',
    'temp_bs_lag22',
    'temp_bs_lag23',
    'temp_bs_lag24',
    'temp_bs_lag31',
    'temp_bs_lag32',
    'temp_bs_lag33',
    'temp_bs_lag34',
    'temp_bs_lag41',
    'temp_bs_lag42',
    'temp_bs_lag43',
    'temp_bs_lag44',
    'temp_bs_lag51',
    'temp_bs_lag52',
    'temp_bs_lag53',
    'temp_bs_lag54',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha',
    'gdp_per_capita',
    'births'
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
  )
)


if (model_idx < 1 || model_idx > length(models_with_temp)) {
  stop(paste("model_idx", model_idx, "out of range — must be 1 to", length(models_with_temp)))
}

model_name <- names(models_with_temp)[model_idx]
model_vars <- models_with_temp[[model_idx]]
cat("Model:", model_name, "\n\n")

# ---------------------------------------------------------------
# 4. Load and prepare training data
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# ── Municipality-level posterior immunity draws ─────────────────────────────────
# Replace the old macro-region-level Pr_0priorinf / Pr_1priorinf with
# municipality-level (city_residency x year) values that carry posterior
# uncertainty — see 01_model_fit/build_immunity_city_draws.R. This task's
# bootstrap_idx is paired with immunity draw sample_id (cycling through
# n_immunity_samples), and that same draw is held fixed across every ensemble
# member below — immunity isn't part of the climate counterfactual, only the
# state resampling and the immunity draw vary together per bootstrap_idx.
city_health_region <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency")

immunity_draws_sampled <- read_csv("../data/immunity_draws_sampled1000.csv", show_col_types = FALSE)
n_immunity_samples      <- max(immunity_draws_sampled$sample_id)
this_sample_id          <- ((bootstrap_idx - 1) %% n_immunity_samples) + 1

immunity_this_draw <- immunity_draws_sampled %>%
  filter(sample_id == this_sample_id) %>%
  select(code_health_region, year, Pr_0priorinf, Pr_1priorinf, Pr_2priorinf)

dengue_temp <- dengue_temp %>%
  left_join(immunity_this_draw, by = c("code_health_region", "year"))

bootstrap_states <- read_csv2(bootstrap_path, show_col_types = FALSE) %>%
  mutate(states_list = map(resampled_states, ~str_split(.x, ",")[[1]]))

cat("Loaded", nrow(bootstrap_states), "bootstrap state rows\n")

# ---------------------------------------------------------------
# 5. Fit B-splines on full training data (knots fixed here,
#    consistent across all bootstrap/ensemble predictions)
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
cat("Base data prepared with B-splines!\n\n")

# ---------------------------------------------------------------
# 6. Model formula
# ---------------------------------------------------------------
create_model_formula <- function(vars) {
  fixed_effects <- intersect(vars, c('month', 'year', 'city_residency',
                                     'month_childs', 'year_childs',
                                     'year_region', 'month_region'))
  predictors <- setdiff(vars, fixed_effects)
  formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "))
  if (length(fixed_effects) > 0) {
    formula_str <- paste(formula_str, "|", paste(fixed_effects, collapse = " + "))
  }
  as.formula(formula_str)
}

model_formula <- create_model_formula(model_vars)
cat("Model formula:", deparse(model_formula), "\n\n")

# ---------------------------------------------------------------
# 7. Build bootstrap sample for this task's bootstrap_idx
# ---------------------------------------------------------------
resampled_states <- bootstrap_states$states_list[[bootstrap_idx]]
cat("Resampled states:", paste(resampled_states, collapse = ", "), "\n")

ensemble_bootstrap_ids <- data.frame(
  state_residency = resampled_states,
  stringsAsFactors = FALSE
) %>%
  mutate(bootstrap_id = 1:n())

bootstrap_training_data <- dengue_temp %>%
  left_join(ensemble_bootstrap_ids,
            by = "state_residency",
            relationship = "many-to-many") %>%
  filter(!is.na(bootstrap_id))

cat("Bootstrap training data:", nrow(bootstrap_training_data), "rows\n")

# ---------------------------------------------------------------
# 8. Fit one model PER REGION on this task's bootstrap sample
#    (region-only refit — separate coefficients per region, same model
#    formula, mirrors 02_fit_national_vs_regional.R's regional fits)
# ---------------------------------------------------------------
regions_all <- sort(unique(dengue_temp$region))
cat("Regions:", paste(regions_all, collapse = ", "), "\n")

region_models <- list()
for (region_name in regions_all) {
  region_train <- bootstrap_training_data %>% filter(region == region_name)

  region_models[[region_name]] <- tryCatch({
    fixest::fepois(
      fml           = model_formula,
      offset        = ~log_pop_offset,
      data          = region_train,
      combine.quick = FALSE
    )
  }, error = function(e) {
    cat("  ERROR fitting region", region_name, ":", e$message, "\n")
    NULL
  })

  if (!is.null(region_models[[region_name]])) {
    cat("  Fitted region:", region_name, "(", nrow(region_train), "rows )\n")
  }
}
cat("Region models fitted!\n\n")

# ---------------------------------------------------------------
# 9. Apply each region's model to that region's rows, for all ensemble members
# ---------------------------------------------------------------
file_list    <- paste0(ensemble_dir, "/ext-nat-ens", sprintf("%03d", 0:524), ".csv")
results_list <- vector("list", length(file_list))

for(ensemble_idx in seq_along(file_list)) {
  file_path <- file_list[ensemble_idx]
  cat("  Ensemble member", ensemble_idx, "of", length(file_list), "\n")
  
  tryCatch({
    new_data <- read_csv(file_path, show_col_types = FALSE)
    new_data$log_pop_offset <- log(new_data$population / 100000)
    new_data <- add_spline_columns(new_data)

    # Same fixed immunity draw as the training data (this_sample_id) — only
    # climate is counterfactual across ensemble members, immunity is held
    # at the observed posterior draw for this bootstrap_idx.
    new_data <- new_data %>%
      select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
      left_join(city_health_region, by = "city_residency") %>%
      left_join(immunity_this_draw, by = c("code_health_region", "year"))

    new_data <- new_data %>%
      left_join(ensemble_bootstrap_ids,
                by = "state_residency",
                relationship = "many-to-many") %>%
      filter(!is.na(bootstrap_id))
    
    new_data <- new_data %>%
      filter(
        mean_2m_air_temp_degree1_lag1 > tlimit &
        mean_2m_air_temp_degree1_lag2 > tlimit &
        mean_2m_air_temp_degree1_lag3 > tlimit &
        mean_2m_air_temp_degree1_lag4 > tlimit &
        mean_2m_air_temp_degree1_lag5 > tlimit
      )
    
    # Predict each region's rows with that region's own model — a region
    # whose model failed to fit contributes no predictions (dropped below),
    # not a fallback to some other region's coefficients.
    new_data$pred_cases <- NA_real_
    for (region_name in names(region_models)) {
      model_r <- region_models[[region_name]]
      if (is.null(model_r)) next
      idx <- which(new_data$region == region_name)
      if (length(idx) == 0) next
      new_data$pred_cases[idx] <- tryCatch(
        predict(model_r, newdata = new_data[idx, ]),
        error = function(e) NA_real_
      )
    }
    new_data <- new_data %>% filter(!is.na(pred_cases))

    aggregated_data <- new_data %>%
      mutate(
        pred_incidence = pred_cases / population * 100000
      ) %>%
      group_by(region, date_first_symptoms) %>%
      summarise(
        total_population   = sum(population,                       na.rm = TRUE),
        total_pred_cases   = sum(pred_cases,                       na.rm = TRUE),
        total_actual_cases = sum(dengue_inc * population / 100000, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        predicted_incidence = total_pred_cases   / total_population * 100000,
        actual_incidence    = total_actual_cases / total_population * 100000,
        model_type          = model_name,
        bootstrap_iteration = bootstrap_idx,
        ensemble_member     = ensemble_idx
      ) %>%
      select(model_type, region, date_first_symptoms, total_population, total_pred_cases,
             predicted_incidence, actual_incidence, bootstrap_iteration, ensemble_member) %>%
      arrange(region, date_first_symptoms)
    
    results_list[[ensemble_idx]] <- aggregated_data
    
  }, error = function(e) {
    cat("    ERROR on ensemble member", ensemble_idx, ":", e$message, "\n")
    results_list[[ensemble_idx]] <<- data.frame(
      model_type          = model_name,
      region              = NA,
      date_first_symptoms = as.Date(NA),
      total_population    = NA,
      total_pred_cases    = NA,
      predicted_incidence = NA,
      actual_incidence    = NA,
      bootstrap_iteration = bootstrap_idx,
      ensemble_member     = ensemble_idx
    )
  })
}

# ---------------------------------------------------------------
# 10. Save this task's results
#     Filename encodes model index and bootstrap index so outputs
#     from different models/bootstraps never collide.
# ---------------------------------------------------------------
final_results <- bind_rows(results_list) %>% filter(!is.na(region))

out_file <- sprintf("predictions_model%02d_bootstrap%04d.csv", model_idx, bootstrap_idx)
out_path <- file.path(output_dir, out_file)
write_csv(final_results, out_path)
cat("Saved results to:", out_path, "\n")
cat("Done!\n")
