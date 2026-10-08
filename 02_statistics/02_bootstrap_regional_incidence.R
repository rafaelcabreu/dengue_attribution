# bootstrap_regional_incidence.R
#
# State block bootstrap for the two models shown in Figure 2's regional
# time-series panel:
#   - "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity"
#     (fixed effects model, i.e. C15-Yr-Mr-c(PC,SR,S,IM))
#   - "Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity"
#     (no fixed effects model, i.e. C15-c(PC,SR,S,IM))
#
# For each bootstrap iteration, the model is refit on state-block-resampled
# data (same presampled draws used by 02_bootstrap_parallel_full_sqcorr.R),
# and the fitted values are aggregated to region x date_first_symptoms on
# THAT SAME resampled sample (not repredicted onto the full original data).
#
# This mirrors how the existing R^2/AIC block bootstrap evaluates fit
# quality on the resampled sample itself, and avoids the fact that a
# state-block bootstrap can drop entire states (and therefore all of a
# state's cities) from a given draw: those cities' city_residency fixed
# effect would not exist in that draw's fit, so predicting them from the
# full original data would silently propagate NA/undefined fitted counts.
# Aggregating cases/population/incidence within the same resampled sample
# sidesteps that: pred_regional_incidence is a ratio (cases / population)
# so it stays comparable across draws despite states appearing 0-3+ times.
#
# Run this AFTER unzipping data/bootstrap_state_samples.csv.zip to
# data/bootstrap_state_samples.csv.
#
# Array mapping (n_models=2, n_chunks=10, chunk_size=100):
#   Total tasks = n_models * n_chunks = 20
#   task 1..10  -> model 1 (fixed effects),    chunks 1..10 (boots 1-100, ..., 901-1000)
#   task 11..20 -> model 2 (no fixed effects),  chunks 1..10
#
# Each task writes:
#   data/statistics/bootstrap_chunks_regional_incidence/regional_incidence_model<MM>_chunk<CCC>.csv
#
# After all tasks complete, run:
#   Rscript 02_aggregate_regional_incidence_bootstrap.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 02_bootstrap_regional_incidence.R <task_id>")

task_id    <- as.integer(args[1])
n_chunks   <- 10     # chunks per model
chunk_size <- 100    # bootstrap iterations per chunk — n_chunks * chunk_size = n_boot (1000)
n_boot     <- n_chunks * chunk_size

model_idx <- ceiling(task_id / n_chunks)
chunk_idx <- task_id - (model_idx - 1) * n_chunks

boot_start <- (chunk_idx - 1) * chunk_size + 1
boot_end   <- chunk_idx * chunk_size

cat(sprintf("Task %d -> model_idx=%d, chunk_idx=%d, bootstraps %d-%d\n",
            task_id, model_idx, chunk_idx, boot_start, boot_end))

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
data_path      <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"
bootstrap_path <- "../data/bootstrap_state_samples.csv"
output_dir     <- "../data/statistics/bootstrap_chunks_regional_incidence"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definitions — must match Figure2.ipynb / fit_model.R exactly
# ---------------------------------------------------------------
models_with_temp <- list(
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'temp_bs_lag11', 'temp_bs_lag12', 'temp_bs_lag13', 'temp_bs_lag14',
    'temp_bs_lag21', 'temp_bs_lag22', 'temp_bs_lag23', 'temp_bs_lag24',
    'temp_bs_lag31', 'temp_bs_lag32', 'temp_bs_lag33', 'temp_bs_lag34',
    'temp_bs_lag41', 'temp_bs_lag42', 'temp_bs_lag43', 'temp_bs_lag44',
    'temp_bs_lag51', 'temp_bs_lag52', 'temp_bs_lag53', 'temp_bs_lag54',
    'total_precipitation_lag1', 'total_precipitation_lag2',
    'total_precipitation_lag3', 'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1', 'immunity_lag2', 'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha', 'gdp_per_capita', 'births',
    'Pr_0priorinf', 'Pr_1priorinf', 'Pr_2priorinf',
    'year_region', 'month_region'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity' = c(
    'temp_bs_lag11', 'temp_bs_lag12', 'temp_bs_lag13', 'temp_bs_lag14',
    'temp_bs_lag21', 'temp_bs_lag22', 'temp_bs_lag23', 'temp_bs_lag24',
    'temp_bs_lag31', 'temp_bs_lag32', 'temp_bs_lag33', 'temp_bs_lag34',
    'temp_bs_lag41', 'temp_bs_lag42', 'temp_bs_lag43', 'temp_bs_lag44',
    'temp_bs_lag51', 'temp_bs_lag52', 'temp_bs_lag53', 'temp_bs_lag54',
    'total_precipitation_lag1', 'total_precipitation_lag2',
    'total_precipitation_lag3', 'total_precipitation_lag4',
    'total_precipitation_lag5',
    'city_residency',
    'immunity_lag1', 'immunity_lag2', 'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha', 'gdp_per_capita', 'births',
    'Pr_0priorinf', 'Pr_1priorinf', 'Pr_2priorinf'
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

# ── Municipality-level posterior immunity draws ─────────────────────────────────
# Replace the old macro-region-level Pr_0priorinf / Pr_1priorinf with
# municipality-level (city_residency x year) values that carry posterior
# uncertainty — see 01_model_fit/build_immunity_city_draws.R. Each bootstrap
# iteration below is paired with one posterior draw (cycling through
# n_immunity_samples), the same pairing scheme used in the R²/AIC bootstrap.
city_health_region <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency")

immunity_draws_sampled <- read_csv("../data/immunity_draws_sampled1000.csv", show_col_types = FALSE)
n_immunity_samples      <- max(immunity_draws_sampled$sample_id)

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
# 6. Helper functions
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

# Fit on the resampled draw and aggregate fitted counts to region x date
# WITHIN THAT SAME DRAW (see header comment for why we don't repredict
# onto the full, unresampled dataset).
boot_fit_regional_incidence <- function(df_full, bootstrap_iter, bootstrap_data,
                                        model_vars, model_formula,
                                        state_id_var = "state_residency",
                                        immunity_data = NULL, n_immunity_samples = NULL) {
  boot_data <- boot_strat_presampled(df_full, bootstrap_iter,
                                     bootstrap_data, state_id_var)

  if (nrow(boot_data) == 0)                                  return(NULL)

  # Pair this bootstrap iteration with one posterior immunity draw (cycling
  # through the n_immunity_samples available), so state-resampling and
  # immunity uncertainty vary together across the same set of iterations.
  if (!is.null(immunity_data)) {
    this_sample_id <- ((bootstrap_iter - 1) %% n_immunity_samples) + 1
    boot_data <- boot_data %>%
      left_join(
        immunity_data %>% filter(sample_id == this_sample_id) %>%
          select(code_health_region, year, Pr_0priorinf, Pr_1priorinf, Pr_2priorinf),
        by = c("code_health_region", "year")
      )
  }

  if (length(setdiff(model_vars, names(boot_data))) > 0)     return(NULL)
  if (!"n_cases"     %in% names(boot_data))                  return(NULL)
  if (!"population"  %in% names(boot_data))                  return(NULL)

  fit <- tryCatch({
    fixest::fepois(model_formula,
                   offset   = ~log_pop_offset,
                   data     = boot_data,
                   nthreads = 1)
  }, error = function(e) NULL)

  if (is.null(fit)) return(NULL)

  fitted_vals <- tryCatch(predict(fit, sample = "original"), error = function(e) NULL)
  if (is.null(fitted_vals)) return(NULL)

  boot_data %>%
    mutate(.fitted_cases = as.numeric(fitted_vals)) %>%
    filter(!is.na(.fitted_cases)) %>%
    group_by(region, date_first_symptoms) %>%
    summarise(
      total_actual_cases = sum(n_cases,        na.rm = TRUE),
      total_pred_cases   = sum(.fitted_cases,  na.rm = TRUE),
      total_population   = sum(population,     na.rm = TRUE),
      n_observations     = n(),
      n_cities            = n_distinct(city_residency),
      .groups = "drop"
    ) %>%
    mutate(
      actual_regional_incidence = total_actual_cases / total_population * 100000,
      pred_regional_incidence   = total_pred_cases   / total_population * 100000,
      bootstrap_iteration       = bootstrap_iter
    )
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
    boot_fit_regional_incidence(df_full             = dengue_temp,
                                bootstrap_iter       = i,
                                bootstrap_data       = bootstrap_states,
                                model_vars           = model_vars,
                                model_formula        = model_formula,
                                state_id_var         = "state_residency",
                                immunity_data        = immunity_draws_sampled,
                                n_immunity_samples    = n_immunity_samples)
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
  out_file <- sprintf("regional_incidence_model%02d_chunk%03d.csv", model_idx, chunk_idx)
  out_path <- file.path(output_dir, out_file)
  write_csv(chunk_df, out_path)
  cat("Saved results to:", out_path, "\n")
} else {
  cat("WARNING: No successful iterations — no output file written.\n")
}

cat("Done!\n")
