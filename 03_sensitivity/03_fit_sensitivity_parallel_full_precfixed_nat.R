# model_slurm_parallel.R
#
# Single-fit (all regions together) version with two models,
# designed to run one bootstrap iteration x one model per SLURM array task.
#
# The SLURM array index encodes BOTH bootstrap_idx and model_idx.
# With n_bootstrap=100 and n_models=2, total tasks = 200:
#   task 1..100   -> model 1, bootstrap 1..100
#   task 101..200 -> model 2, bootstrap 1..100
#
# This mapping is handled automatically below.
#
# Submit with:
#   sbatch submit_bootstrap.sh
#
# After all tasks complete, combine outputs with:
#   Rscript aggregate_results.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript model_slurm_parallel.R <task_id>")

task_id     <- as.integer(args[1])
n_bootstrap <- 100
n_models    <- 1

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
data_path      <- "data/model_input_brazil_immunity_city_with_priorinf.csv"
bootstrap_path <- "data/bootstrap_state_samples.csv"
ensemble_dir   <- "data/predict-fixed-prec"
output_dir     <- "data/sensitivity-fixedprec-nat/bootstrap_chunks"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tlimit    <- 12
n_bs_cols <- 4

# ---------------------------------------------------------------
# 3. Model definitions
# ---------------------------------------------------------------
models_with_temp <- list(
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
# 8. Fit model on bootstrap sample
# ---------------------------------------------------------------
bootstrap_model <- fixest::fepois(
  fml           = model_formula,
  offset        = ~log_pop_offset,
  data          = bootstrap_training_data,
  combine.quick = FALSE
)
cat("Bootstrap model fitted successfully!\n\n")

# ---------------------------------------------------------------
# 9. Apply bootstrap model to all ensemble members
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
    
    pred_cases <- predict(bootstrap_model, newdata = new_data)
    
    aggregated_data <- new_data %>%
      mutate(
        pred_cases     = pred_cases,
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