# bootstrap_capital_incidence.R
#
# State block bootstrap for the fixed-effects model shown in FigureS12.ipynb's
# municipality time-series panel (observed vs. fitted incidence, one capital
# city per region):
#   "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity"
#   (i.e. C15-Yr-Mr-c(PC,SR,S,IM))
#
# Sibling script to 02_bootstrap_regional_incidence.R — same state-block
# resampling procedure, same "aggregate within the same resampled draw"
# logic (see that script's header for the full rationale), but aggregated
# to city_residency x date_first_symptoms instead of region x
# date_first_symptoms, and filtered down to a handful of capital cities
# before writing output (the model is still fit on the full resampled
# national data; only the output rows are filtered).
#
# NOTE: a state-block bootstrap resamples STATES, not individual cities.
# For any one of these single-city capitals, its home state is absent from
# roughly 1/e (~37%) of draws (classic bootstrap "not selected" rate for
# resampling with replacement) — those draws contribute no row for that
# city. This is expected, not a bug: quantiles below are computed over
# whatever number of valid draws each city-month actually has (reported as
# n_valid_iterations), which is still ~630 of 1000 draws on average.
#
# Run this AFTER unzipping data/bootstrap_state_samples.csv.zip to
# data/bootstrap_state_samples.csv.
#
# Array mapping (n_chunks=10, chunk_size=100):
#   Total tasks = n_chunks = 10
#   task 1..10 -> chunks 1..10 (boots 1-100, 101-200, ..., 901-1000)
#
# Each task writes:
#   data/statistics/bootstrap_chunks_capital_incidence/capital_incidence_chunk<CCC>.csv
#
# After all tasks complete, run:
#   Rscript 02_aggregate_capital_incidence_bootstrap.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 02_bootstrap_capital_incidence.R <task_id>")

chunk_idx  <- as.integer(args[1])
n_chunks   <- 10     # chunks total
chunk_size <- 100    # bootstrap iterations per chunk — n_chunks * chunk_size = n_boot (1000)
n_boot     <- n_chunks * chunk_size

boot_start <- (chunk_idx - 1) * chunk_size + 1
boot_end   <- chunk_idx * chunk_size

cat(sprintf("Task %d -> chunk_idx=%d, bootstraps %d-%d\n",
            chunk_idx, chunk_idx, boot_start, boot_end))

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
data_path      <- "data/model_input_brazil_immunity_city_with_priorinf.csv"
bootstrap_path <- "data/bootstrap_state_samples.csv"
output_dir     <- "data/statistics/bootstrap_chunks_capital_incidence"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definition — must match FigureS12.ipynb / fit_model.R exactly
# ---------------------------------------------------------------
model_name <- 'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity'
model_vars <- c(
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
  'Pr_0priorinf', 'Pr_1priorinf',
  'year_region', 'month_region'
)

# One capital city per region — must match FigureS12.ipynb's `capitals` dict
target_cities <- c('AM_MANAUS', 'BA_SALVADOR', 'DF_BRASILIA', 'SP_SAO PAULO', 'RS_PORTO ALEGRE')

cat("Model:", model_name, "\n\n")

# ---------------------------------------------------------------
# 4. Load and prepare data
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

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

# Fit on the resampled draw and aggregate fitted counts to city x date
# WITHIN THAT SAME DRAW (see header comment for why — mirrors
# 02_bootstrap_regional_incidence.R's boot_fit_regional_incidence()).
# Only the target capital cities are kept in the returned rows; the model
# itself is still fit on the full resampled national data.
boot_fit_capital_incidence <- function(df_full, bootstrap_iter, bootstrap_data,
                                       model_vars, model_formula, target_cities,
                                       state_id_var = "state_residency") {
  boot_data <- boot_strat_presampled(df_full, bootstrap_iter,
                                     bootstrap_data, state_id_var)

  if (nrow(boot_data) == 0)                                  return(NULL)
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
    filter(city_residency %in% target_cities) %>%
    group_by(city_residency, date_first_symptoms) %>%
    summarise(
      total_actual_cases = sum(n_cases,        na.rm = TRUE),
      total_pred_cases   = sum(.fitted_cases,  na.rm = TRUE),
      total_population   = sum(population,     na.rm = TRUE),
      n_observations     = n(),
      .groups = "drop"
    ) %>%
    mutate(
      actual_capital_incidence = total_actual_cases / total_population * 100000,
      pred_capital_incidence   = total_pred_cases   / total_population * 100000,
      bootstrap_iteration      = bootstrap_iter
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
    boot_fit_capital_incidence(df_full        = dengue_temp,
                               bootstrap_iter = i,
                               bootstrap_data = bootstrap_states,
                               model_vars     = model_vars,
                               model_formula  = model_formula,
                               target_cities  = target_cities,
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
  out_file <- sprintf("capital_incidence_chunk%03d.csv", chunk_idx)
  out_path <- file.path(output_dir, out_file)
  write_csv(chunk_df, out_path)
  cat("Saved results to:", out_path, "\n")
} else {
  cat("WARNING: No successful iterations — no output file written.\n")
}

cat("Done!\n")
