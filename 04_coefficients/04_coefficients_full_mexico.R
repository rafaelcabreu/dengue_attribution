# coefficients_full_mexico.R
#
# Mexico-data counterpart to 04_coefficients_full.R — full 1000-iteration
# state-block bootstrap (same resampling mechanism: boot_strat_newID() /
# boot_fit_model(), seeded per global iteration index for reproducibility),
# chunked for a SLURM array exactly like the original script.
#
# Two models, both the same formula
# 'Climate(lag 1-5) + Year|region + Month|region', differing only in which
# countries' rows are kept for training:
#   model 1: country %in% c('Brazil', 'Mexico')
#   model 2: country == 'Mexico'
#
# Array mapping (n_models=2, n_chunks=10, chunk_size=100):
#   Total tasks = n_models * n_chunks = 20
#   task 1..10  -> model 1 (Brazil + Mexico), chunks 1..10 (boots 1-100, ..., 901-1000)
#   task 11..20 -> model 2 (Mexico only),     chunks 1..10
#
# Each task writes:
#   data/coefficients/bootstrap_chunks_mexico/coefs_model<MM>_chunk<CCC>.csv
#
# After all tasks complete, run:
#   Rscript 04_aggregate_coefficients_mexico.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 04_coefficients_full_mexico.R <task_id>")

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
data_path  <- "data/model_input_mexico-peru-brazil_immunity.csv"
output_dir <- "data/coefficients/bootstrap_chunks_mexico"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definitions
#    Same formula for both — only the training-data country filter differs.
# ---------------------------------------------------------------
model_formula_vars <- c(
  'temp_bs_lag11', 'temp_bs_lag12', 'temp_bs_lag13', 'temp_bs_lag14',
  'temp_bs_lag21', 'temp_bs_lag22', 'temp_bs_lag23', 'temp_bs_lag24',
  'temp_bs_lag31', 'temp_bs_lag32', 'temp_bs_lag33', 'temp_bs_lag34',
  'temp_bs_lag41', 'temp_bs_lag42', 'temp_bs_lag43', 'temp_bs_lag44',
  'temp_bs_lag51', 'temp_bs_lag52', 'temp_bs_lag53', 'temp_bs_lag54',
  'total_precipitation_lag1', 'total_precipitation_lag2',
  'total_precipitation_lag3', 'total_precipitation_lag4',
  'total_precipitation_lag5',
  'month_region', 'year_region', 'city_residency'
)

models <- list(
  'Climate(lag 1-5) + Year|region + Month|region (Brazil + Mexico)' = c('Brazil', 'Mexico'),
  'Climate(lag 1-5) + Year|region + Month|region (Mexico only)'     = c('Mexico')
)

n_models <- length(models)
if (model_idx < 1 || model_idx > n_models) {
  stop(paste("model_idx", model_idx, "out of range — must be 1 to", n_models))
}

model_name      <- names(models)[model_idx]
model_countries <- models[[model_idx]]
cat("Model:", model_name, "\n")
cat("Training countries:", paste(model_countries, collapse = ", "), "\n\n")

# ---------------------------------------------------------------
# 4. Load and prepare data (this model's country subset only)
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# The combined Mexico-Peru-Brazil file has `region` but not the region x
# year / region x month interaction columns the Brazil-only file already
# has precomputed — build them here.
dengue_temp <- dengue_temp %>%
  mutate(
    year_region  = paste(year, region, sep = "_"),
    month_region = paste(month, region, sep = "_")
  ) %>%
  filter(country %in% model_countries)

cat("Training rows after country filter:", nrow(dengue_temp), "\n")

# Fit B-splines on this model's full (filtered) training data — consistent
# knots across all bootstrap iterations for this model.
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

# Unique states for bootstrap sampling (within this model's country subset)
df_states <- dengue_temp %>%
  select(state_residency) %>%
  distinct() %>%
  filter(!is.na(state_residency))

cat("Number of unique states:", nrow(df_states), "\n")

# ---------------------------------------------------------------
# 5. Helper functions (identical logic to 04_coefficients_full.R)
# ---------------------------------------------------------------
boot_strat_newID <- function(df_ids, df_full,
                             state_id_var = "state_residency",
                             seed = 1234) {
  set.seed(seed)
  ids <- df_ids %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    mutate(state_boot_id = 1:n())

  df_full %>%
    left_join(ids, by = state_id_var, relationship = "many-to-many") %>%
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

boot_fit_model <- function(df_ids, df_full, model_vars,
                           state_id_var = "state_residency",
                           seed = 1234) {
  boot_data     <- boot_strat_newID(df_ids, df_full, state_id_var, seed)
  model_formula <- build_formula(model_vars)
  model <- fixest::fepois(model_formula,
                          offset   = ~log_pop_offset,
                          data     = boot_data,
                          nthreads = 1)
  coef(model)
}

# ---------------------------------------------------------------
# 6. Run bootstrap iterations for this chunk
# ---------------------------------------------------------------
model_formula <- build_formula(model_formula_vars)
cat("Model formula:", deparse(model_formula), "\n\n")
cat(sprintf("Running bootstrap iterations %d to %d...\n", boot_start, boot_end))

chunk_results <- list()

for (i in boot_start:boot_end) {
  if ((i - boot_start + 1) %% 10 == 0) {
    cat("  Completed", i - boot_start + 1, "of", chunk_size, "iterations in this chunk\n")
  }

  result <- tryCatch({
    boot_fit_model(df_ids         = df_states,
                   df_full        = dengue_temp,
                   model_vars     = model_formula_vars,
                   state_id_var   = "state_residency",
                   seed           = i)          # seed = global iteration index for reproducibility
  }, error = function(e) {
    cat("  Error in iteration", i, ":", e$message, "\n")
    NULL
  })

  if (!is.null(result)) {
    result_df <- as.data.frame(t(result))
    result_df$boot_iteration <- i
    result_df$model_name     <- model_name
    chunk_results[[length(chunk_results) + 1]] <- result_df
  }
}

cat("Successful iterations in this chunk:", length(chunk_results), "of", chunk_size, "\n")

# ---------------------------------------------------------------
# 7. Save this chunk's results
# ---------------------------------------------------------------
if (length(chunk_results) > 0) {
  chunk_df <- bind_rows(chunk_results)
  out_file <- sprintf("coefs_model%02d_chunk%03d.csv", model_idx, chunk_idx)
  out_path <- file.path(output_dir, out_file)
  write_csv(chunk_df, out_path)
  cat("Saved results to:", out_path, "\n")
} else {
  cat("WARNING: No successful iterations — no output file written.\n")
}

cat("Done!\n")
