# rsq_aic_bootstrap_slurm_parallel.R
#
# Parallelised R-squared/AIC bootstrap across 6 active models x 1000 iterations,
# using pre-sampled bootstrap states from a CSV file.
#
# Each task runs a CHUNK of bootstrap iterations for one model:
#   n_models=6, n_chunks=10, chunk_size=100 -> 60 total tasks
#   task 1..10  -> model 1, bootstraps 1..100, 101..200, ..., 901..1000
#   task 11..20 -> model 2, bootstraps 1..100, ...
#   ...
#   task 51..60 -> model 6, bootstraps 1..100, ...
#
# Each task writes:
#   bootstrap_chunks_rsq/rsq_model<MM>_chunk<CCC>.csv
#
# After all tasks complete, run:
#   Rscript aggregate_rsq_aic_bootstrap.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript rsq_aic_bootstrap_slurm_parallel.R <task_id>")

task_id    <- as.integer(args[1])
n_chunks   <- 10    # chunks per model
chunk_size <- 100   # bootstrap iterations per chunk; n_chunks * chunk_size = n_boot (1000)

# Derive model_idx and chunk_idx from flat task_id
model_idx <- ceiling(task_id / n_chunks)
chunk_idx <- task_id - (model_idx - 1) * n_chunks

# Bootstrap index range for this chunk
boot_start <- (chunk_idx - 1) * chunk_size + 1
boot_end   <- chunk_idx * chunk_size

cat(sprintf("Task %d -> model_idx=%d, chunk_idx=%d, bootstraps %d-%d\n",
            task_id, model_idx, chunk_idx, boot_start, boot_end))

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
data_path      <- "data/model_input_brazil_immunity_city_with_priorinf.csv"
bootstrap_path <- "data/bootstrap_state_samples.csv"
#output_dir     <- "~/cpdn_nonnerc/aaim/dengue/statistics/bootstrap_chunks_rsq"
output_dir     <- "data/statistics/bootstrap_chunks_rsq"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definitions
#    Only the active (uncommented) models from the original script.
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
    'Pr_1priorinf'
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
    'Pr_1priorinf',
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
    'Pr_1priorinf'
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
    'Pr_1priorinf'
  ),
  'Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'month_region',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha',
    'gdp_per_capita',
    'births',
    'Pr_0priorinf',
    'Pr_1priorinf'
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
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio' = c(
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
    'year_region',
    'month_region'
  )
)

n_models <- length(models_with_temp)
if (model_idx < 1 || model_idx > n_models) {
  stop(paste("model_idx", model_idx, "out of range — must be 1 to", n_models))
}

model_name <- names(models_with_temp)[model_idx]
model_vars <- models_with_temp[[model_idx]]
cat("Model:", model_name, "\n\n")

# ---------------------------------------------------------------
# 4. Load and prepare data
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# Fit B-splines on full data (consistent knots across all tasks)
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

cat("Data prepared with B-splines!\n\n")

# ---------------------------------------------------------------
# 5. Load pre-sampled bootstrap states
#    Each task only needs the rows for its own bootstrap range,
#    but loading the full file is cheap and keeps the code simple.
# ---------------------------------------------------------------
bootstrap_states <- read_csv2(
  bootstrap_path,
  col_names = TRUE,
  locale    = locale(encoding = "UTF-8"),
  n_max     = 1000,
  show_col_types = FALSE
) %>%
  separate_rows(resampled_states, sep = ",") %>%
  mutate(resampled_states = str_trim(resampled_states)) %>%
  filter(!is.na(resampled_states), resampled_states != "") %>%
  group_by(bootstrap_iteration) %>%
  mutate(state_boot_id = row_number()) %>%
  ungroup()

n_boot_available <- max(bootstrap_states$bootstrap_iteration)
cat("Pre-sampled bootstrap iterations available:", n_boot_available, "\n")

if (boot_end > n_boot_available) {
  cat(sprintf("WARNING: boot_end (%d) > available iterations (%d); capping.\n",
              boot_end, n_boot_available))
  boot_end <- n_boot_available
}

# ---------------------------------------------------------------
# 6. Helper functions (identical logic to original serial script)
# ---------------------------------------------------------------
boot_strat_presampled <- function(df_full, bootstrap_iter, bootstrap_data,
                                  state_id_var = "state_residency") {
  states_for_iter <- bootstrap_data %>%
    filter(bootstrap_iteration == bootstrap_iter) %>%
    select(state_boot_id, state_name = resampled_states)

  if (!state_id_var %in% colnames(df_full)) return(data.frame())

  join_by <- setNames("state_name", state_id_var)

  df_full %>%
    inner_join(states_for_iter, by = join_by, relationship = "many-to-many") %>%
    filter(!is.na(state_boot_id))
}

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

boot_fit_model_rsq <- function(df_full, bootstrap_iter, bootstrap_data,
                               model_vars, model_formula,
                               state_id_var = "state_residency") {
  boot_data <- boot_strat_presampled(df_full, bootstrap_iter,
                                     bootstrap_data, state_id_var)

  if (nrow(boot_data) == 0)                                  return(NULL)
  if (length(setdiff(model_vars, names(boot_data))) > 0)     return(NULL)
  if (!"n_cases"     %in% names(boot_data))                  return(NULL)
  if (!"population"  %in% names(boot_data))                  return(NULL)

  tryCatch({
    model <- fixest::fepois(model_formula,
                            offset   = ~log_pop_offset,
                            data     = boot_data,
                            nthreads = 1)
    data.frame(
      r_squared        = r2(model, type = "cor2"),
      pseudo_r_squared = r2(model, type = "pr2"),
      aic              = AIC(model),
      n_obs            = model$nobs,
      bootstrap_iteration = bootstrap_iter
    )
  }, error = function(e) NULL)
}

# ---------------------------------------------------------------
# 7. Run bootstrap iterations for this chunk
# ---------------------------------------------------------------
model_formula <- build_formula(model_vars)
cat("Model formula:", deparse(model_formula), "\n\n")
cat(sprintf("Running bootstrap iterations %d to %d...\n", boot_start, boot_end))

chunk_results <- list()

for (i in boot_start:boot_end) {
  if ((i - boot_start + 1) %% 10 == 0) {
    cat("  Completed", i - boot_start + 1, "of", boot_end - boot_start + 1,
        "iterations in this chunk\n")
  }

  result <- tryCatch({
    boot_fit_model_rsq(df_full        = dengue_temp,
                       bootstrap_iter = i,
                       bootstrap_data = bootstrap_states,
                       model_vars     = model_vars,
                       model_formula  = model_formula,
                       state_id_var   = "state_residency")
  }, error = function(e) {
    cat("  Error in iteration", i, ":", e$message, "\n")
    NULL
  })

  if (!is.null(result)) {
    result$model_name <- model_name
    chunk_results[[length(chunk_results) + 1]] <- result
  }
}

cat("Successful iterations in this chunk:", length(chunk_results),
    "of", boot_end - boot_start + 1, "\n")

# ---------------------------------------------------------------
# 8. Save this chunk's results
# ---------------------------------------------------------------
if (length(chunk_results) > 0) {
  chunk_df <- bind_rows(chunk_results)
  out_file <- sprintf("rsq_model%02d_chunk%03d.csv", model_idx, chunk_idx)
  out_path <- file.path(output_dir, out_file)
  write_csv(chunk_df, out_path)
  cat("Saved results to:", out_path, "\n")
} else {
  cat("WARNING: No successful iterations — no output file written.\n")
}

cat("Done!\n")