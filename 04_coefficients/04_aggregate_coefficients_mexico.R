# aggregate_coefficients_mexico.R
#
# Run this AFTER all SLURM array tasks from 04_coefficients_full_mexico.R
# have completed. Mirrors 04_aggregate_coefficients_full.R: reads all
# per-chunk CSVs from bootstrap_chunks_mexico/, checks for missing tasks,
# and writes the same output files as the original serial script for each
# model:
#   <model_name>_coef_state_blockboot1000.csv
#   <model_name>_coef_state_blockboot_summary1000.csv
# Plus combined files across both models:
#   all_models_bootstrap_results_mexico1000.csv
#   all_models_original_coefficients_mexico.csv

library(tidyverse)
library(fixest)
library(splines)

data_path  <- "data/model_input_mexico-peru-brazil_immunity.csv"
chunks_dir <- "data/coefficients/bootstrap_chunks_mexico"
n_chunks   <- 10
chunk_size <- 100
n_boot     <- n_chunks * chunk_size   # 1000

# ---------------------------------------------------------------
# 1. Model definitions (must match 04_coefficients_full_mexico.R exactly)
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

build_formula <- function(vars) {
  fixed_effects <- intersect(vars, c("city_residency", "month", "year",
                                     "month_childs", "year_childs",
                                     "year_region", "month_region"))
  predictors <- setdiff(vars, fixed_effects)
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
# 2. Load full data once (needed for the "original", non-bootstrapped fits)
# ---------------------------------------------------------------
cat("Loading data for original model fits...\n")
dengue_temp_full <- read_csv(data_path, show_col_types = FALSE)
dengue_temp_full$log_pop_offset <- log(dengue_temp_full$population / 100000)
dengue_temp_full <- dengue_temp_full %>%
  mutate(
    year_region  = paste(year, region, sep = "_"),
    month_region = paste(month, region, sep = "_")
  )

# ---------------------------------------------------------------
# 3. Check which tasks completed
# ---------------------------------------------------------------
cat("\nChecking for completed chunk files...\n")
for (model_idx in 1:n_models) {
  model_name <- names(models)[model_idx]
  missing <- c()
  for (chunk_idx in 1:n_chunks) {
    fname <- sprintf("coefs_model%02d_chunk%03d.csv", model_idx, chunk_idx)
    if (!file.exists(file.path(chunks_dir, fname))) missing <- c(missing, chunk_idx)
  }
  if (length(missing) == 0) {
    cat(sprintf("  Model %02d (%s): all %d chunk files found.\n",
                model_idx, model_name, n_chunks))
  } else {
    cat(sprintf("  Model %02d (%s): MISSING chunks: %s\n",
                model_idx, model_name, paste(missing, collapse = ", ")))
  }
}
cat("\n")

# ---------------------------------------------------------------
# 4. Aggregate per model, fit original model, write outputs
# ---------------------------------------------------------------
all_bootstrap_results <- list()
all_original_results  <- list()

for (model_idx in 1:n_models) {
  model_name      <- names(models)[model_idx]
  model_countries <- models[[model_idx]]
  cat("======= Aggregating model", model_idx, ":", model_name, "=======\n")

  pattern     <- sprintf("coefs_model%02d_chunk\\d{3}\\.csv", model_idx)
  found_files <- sort(list.files(chunks_dir, pattern = pattern, full.names = TRUE))

  if (length(found_files) == 0) {
    cat("  No chunk files found — skipping.\n\n")
    next
  }

  cat("  Reading", length(found_files), "chunk files...\n")
  boot_results <- map_dfr(found_files, ~read_csv(.x, show_col_types = FALSE))
  boot_results$model_name <- model_name

  cat("  Total bootstrap iterations recovered:", nrow(boot_results), "of", n_boot, "\n")

  write_csv(boot_results,
            paste0('data/coefficients/', model_name, "_coef_state_blockboot", n_boot, ".csv"))

  boot_summary <- boot_results %>%
    select(-model_name, -boot_iteration) %>%
    summarise(across(everything(), list(
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x,   na.rm = TRUE),
      q025 = ~quantile(.x, 0.025, na.rm = TRUE),
      q975 = ~quantile(.x, 0.975, na.rm = TRUE)
    ))) %>%
    mutate(model_name = model_name)

  write_csv(boot_summary,
            paste0('data/coefficients/', model_name, "_coef_state_blockboot_summary", n_boot, ".csv"))

  all_bootstrap_results[[model_name]] <- boot_results

  # Fit original (non-bootstrapped) model on this model's country subset
  cat("  Fitting original model...\n")
  tryCatch({
    dengue_temp <- dengue_temp_full %>% filter(country %in% model_countries)

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

    original_model <- fixest::fepois(build_formula(model_formula_vars),
                                     offset   = ~log_pop_offset,
                                     data     = dengue_temp,
                                     nthreads = 1)
    original_coefs <- as.data.frame(t(coef(original_model)))
    original_coefs$model_name <- model_name
    all_original_results[[model_name]] <- original_coefs
    cat("  Original model fitted successfully.\n")
  }, error = function(e) {
    cat("  ERROR fitting original model:", e$message, "\n")
  })

  cat("  Done.\n\n")
}

# ---------------------------------------------------------------
# 5. Write combined outputs across both models
# ---------------------------------------------------------------
if (length(all_bootstrap_results) > 0) {
  combined_bootstrap <- bind_rows(all_bootstrap_results)
  write_csv(combined_bootstrap,
            paste0("all_models_bootstrap_results_mexico", n_boot, ".csv"))
  cat("Combined bootstrap results saved.\n")
}

if (length(all_original_results) > 0) {
  combined_original <- bind_rows(all_original_results)
  write_csv(combined_original, "all_models_original_coefficients_mexico.csv")
  cat("Combined original coefficients saved.\n")
}

cat("\nAggregation complete!\n")
