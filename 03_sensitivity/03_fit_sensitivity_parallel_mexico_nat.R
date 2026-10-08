# fit_sensitivity_parallel_mexico_nat.R
#
# Mexico-data counterpart to 03_fit_sensitivity_parallel_full_nat.R. Identical
# to 03_fit_sensitivity_parallel_mexico_act.R except it applies the fitted
# model to the NATURAL (counterfactual, no-climate-change) ensemble members
# instead of the actual-climate ones. See that script's header for the full
# rationale (no bootstrap; two models differing only by training-data country
# filter).
#
# ASSUMED ensemble file naming (unverified — data/predict-all/ was empty and
# predict-all.tar.xz was truncated in the environment this was written in;
# confirm/fix `ensemble_file_prefix` below once the real files are in place):
#   data/predict-all/ext-mex-nat-ens000.csv ... ext-mex-nat-ens524.csv
#
# Array mapping: task_id IS model_idx directly (no bootstrap dimension).
#   task 1 -> model 1 (Brazil + Mexico training data)
#   task 2 -> model 2 (Mexico only training data)
#
# Each task writes:
#   data/sensitivity-mexico-nat/predictions_model<MM>.csv

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID / model_idx)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 03_fit_sensitivity_parallel_mexico_nat.R <model_idx>")

model_idx <- as.integer(args[1])

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
data_path           <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue_attribution/data/model_input_mexico-peru-brazil_immunity.csv"
ensemble_dir        <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/predict-mexico-peru"
ensemble_file_prefix <- "ext-nat-ens"   # ASSUMED — see header comment
output_dir          <- "/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/data/sensitivity-mexico-nat"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tlimit    <- 12
n_bs_cols <- 4

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
  list(model_name = 'Climate(lag 1-5) + Year|region + Month|region (Brazil + Mexico)',
       countries  = c('Brazil', 'Mexico')),
  list(model_name = 'Climate(lag 1-5) + Year|region + Month|region (Mexico only)',
       countries  = c('Mexico'))
)

if (model_idx < 1 || model_idx > length(models)) {
  stop(paste("model_idx", model_idx, "out of range — must be 1 to", length(models)))
}

model_name       <- models[[model_idx]]$model_name
model_countries  <- models[[model_idx]]$countries
cat("Model:", model_name, "\n")
cat("Training countries:", paste(model_countries, collapse = ", "), "\n\n")

# ---------------------------------------------------------------
# 4. Load and prepare training data
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
  )

dengue_temp <- dengue_temp %>% filter(country %in% model_countries)
cat("Training rows after country filter:", nrow(dengue_temp), "\n")

# ---------------------------------------------------------------
# 5. Fit B-splines on the (filtered) training data
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
cat("Training data prepared with B-splines!\n\n")

# ---------------------------------------------------------------
# 6. Fit the model ONCE on the observed training data (no bootstrap)
# ---------------------------------------------------------------
model_formula <- as.formula(
  paste("n_cases ~",
        paste(setdiff(model_formula_vars, c("month_region", "year_region", "city_residency")), collapse = " + "),
        "|",
        paste(intersect(model_formula_vars, c("month_region", "year_region", "city_residency")), collapse = " + "))
)
cat("Model formula:", deparse(model_formula), "\n\n")

fitted_model <- fixest::fepois(
  fml           = model_formula,
  offset        = ~log_pop_offset,
  data          = dengue_temp,
  combine.quick = FALSE
)
cat("Model fitted successfully!\n\n")

# ---------------------------------------------------------------
# 7. Apply the fitted model to all 525 natural-climate ensemble members
# ---------------------------------------------------------------
file_list    <- paste0(ensemble_dir, "/", ensemble_file_prefix, sprintf("%03d", 0:524), ".csv")
results_list <- vector("list", length(file_list))

for (ensemble_idx in seq_along(file_list)) {
  file_path <- file_list[ensemble_idx]
  cat("  Ensemble member", ensemble_idx, "of", length(file_list), "\n")

  tryCatch({
    new_data <- read_csv(file_path, show_col_types = FALSE)

    if ("country" %in% names(new_data)) {
      new_data <- new_data %>% filter(country %in% model_countries)
    }

    new_data$log_pop_offset <- log(new_data$population / 100000)
    new_data <- new_data %>%
      mutate(
        year_region  = paste(year, region, sep = "_"),
        month_region = paste(month, region, sep = "_")
      )
    new_data <- add_spline_columns(new_data)

    new_data <- new_data %>%
      filter(
        mean_2m_air_temp_degree1_lag1 > tlimit &
        mean_2m_air_temp_degree1_lag2 > tlimit &
        mean_2m_air_temp_degree1_lag3 > tlimit &
        mean_2m_air_temp_degree1_lag4 > tlimit &
        mean_2m_air_temp_degree1_lag5 > tlimit
      )

    pred_cases <- predict(fitted_model, newdata = new_data)

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
        ensemble_member     = ensemble_idx
      ) %>%
      select(model_type, region, date_first_symptoms, total_population, total_pred_cases,
             predicted_incidence, actual_incidence, ensemble_member) %>%
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
      ensemble_member     = ensemble_idx
    )
  })
}

# ---------------------------------------------------------------
# 8. Save this model's results
# ---------------------------------------------------------------
final_results <- bind_rows(results_list) %>% filter(!is.na(region))

out_file <- sprintf("predictions_model%02d.csv", model_idx)
out_path <- file.path(output_dir, out_file)
write_csv(final_results, out_path)
cat("Saved results to:", out_path, "\n")
cat("Done!\n")
