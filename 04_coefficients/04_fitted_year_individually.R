# coef_bootstrap_slurm_by_year.R
#
# Bootstrap sensitivity analysis for a single model, stratified by year.
# Each SLURM array task runs 1000 bootstrap iterations for ONE year of data.
#
# Array mapping:
#   task 1 -> year 1 (first year found in data)
#   task 2 -> year 2
#   ...
#   task N -> year N
#
# The year-to-task mapping is determined by the sorted unique years in the data.
# Check which years are present and set --array=1-N in your submit script accordingly.
#
# Each task writes:
#   bootstrap_by_year/coefs_year<YYYY>.csv
#
# After all tasks complete, run:
#   Rscript aggregate_coef_bootstrap_by_year.R

# ---------------------------------------------------------------
# 0. Parse command-line argument (SLURM_ARRAY_TASK_ID)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript coef_bootstrap_slurm_by_year.R <task_id>")

task_id <- as.integer(args[1])
n_boot  <- 1000

cat(sprintf("Task %d — will run %d bootstrap iterations\n", task_id, n_boot))

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
data_path  <- "data/model_input_brazil_immunity_city_with_priorinf.csv"
output_dir <- "data/coefficients/bootstrap_by_year"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Single model definition
# ---------------------------------------------------------------
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
  'Pr_0priorinf', 'Pr_1priorinf',
  'year_region', 'month_region'
)

# ---------------------------------------------------------------
# 4. Load and prepare full data
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# Fit splines on FULL data so knots are consistent across all year-tasks
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

precip_ns_fitted_lag1 <- ns(dengue_temp$total_precipitation_lag1, df = 3)
precip_ns_fitted_lag2 <- ns(dengue_temp$total_precipitation_lag2, df = 3)
precip_ns_fitted_lag3 <- ns(dengue_temp$total_precipitation_lag3, df = 3)
precip_ns_fitted_lag4 <- ns(dengue_temp$total_precipitation_lag4, df = 3)
precip_ns_fitted_lag5 <- ns(dengue_temp$total_precipitation_lag5, df = 3)

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
    temp_bs_lag53 = temp_bs_fitted_lag5[,3], temp_bs_lag54 = temp_bs_fitted_lag5[,4],
    precip_ns_lag11 = precip_ns_fitted_lag1[,1], precip_ns_lag12 = precip_ns_fitted_lag1[,2], precip_ns_lag13 = precip_ns_fitted_lag1[,3],
    precip_ns_lag21 = precip_ns_fitted_lag2[,1], precip_ns_lag22 = precip_ns_fitted_lag2[,2], precip_ns_lag23 = precip_ns_fitted_lag2[,3],
    precip_ns_lag31 = precip_ns_fitted_lag3[,1], precip_ns_lag32 = precip_ns_fitted_lag3[,2], precip_ns_lag33 = precip_ns_fitted_lag3[,3],
    precip_ns_lag41 = precip_ns_fitted_lag4[,1], precip_ns_lag42 = precip_ns_fitted_lag4[,2], precip_ns_lag43 = precip_ns_fitted_lag4[,3],
    precip_ns_lag51 = precip_ns_fitted_lag5[,1], precip_ns_lag52 = precip_ns_fitted_lag5[,2], precip_ns_lag53 = precip_ns_fitted_lag5[,3]
  )

cat("Data prepared with splines.\n")

# ---------------------------------------------------------------
# 5. Map task_id to a calendar year
# ---------------------------------------------------------------
all_years <- sort(unique(dengue_temp$year))
n_years   <- length(all_years)

cat(sprintf("Years found in data (%d): %s\n", n_years, paste(all_years, collapse = ", ")))
cat(sprintf("Set SLURM --array=1-%d in your submit script.\n\n", n_years))

if (task_id < 1 || task_id > n_years) {
  stop(sprintf("task_id %d out of range — data contains %d years (array should be 1-%d)",
               task_id, n_years, n_years))
}

target_year <- all_years[task_id]
cat(sprintf("Task %d -> fitting on year %d\n\n", task_id, target_year))

# Subset data to the target year only
dengue_year <- dengue_temp %>% filter(year == target_year)
cat(sprintf("Rows for year %d: %d\n", target_year, nrow(dengue_year)))

# Unique states within this year's data for stratified bootstrap
df_states_year <- dengue_year %>%
  select(state_residency) %>%
  distinct() %>%
  filter(!is.na(state_residency))

cat(sprintf("Unique states in year %d: %d\n\n", target_year, nrow(df_states_year)))

# ---------------------------------------------------------------
# 6. Helper functions
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
# 7. Run 1000 bootstrap iterations for this year
# ---------------------------------------------------------------
# Note: year_region and month_region are fixed effects. Within a single year,
# year_region is constant per city so it will be absorbed by city_residency FE.
# If this causes a convergence warning, consider dropping year_region from the
# formula for single-year fits — it carries no within-year variation.
model_formula <- build_formula(model_vars)
cat("Model formula:", deparse(model_formula), "\n\n")
cat(sprintf("Running %d bootstrap iterations for year %d...\n", n_boot, target_year))

year_results <- list()

for (i in seq_len(n_boot)) {
  if (i %% 100 == 0) cat(sprintf("  Iteration %d / %d\n", i, n_boot))

  result <- tryCatch({
    boot_fit_model(df_ids       = df_states_year,
                   df_full      = dengue_year,
                   model_vars   = model_vars,
                   state_id_var = "state_residency",
                   seed         = i)   # seed = iteration index for reproducibility
  }, error = function(e) {
    cat(sprintf("  Error in iteration %d: %s\n", i, e$message))
    NULL
  })

  if (!is.null(result)) {
    result_df <- as.data.frame(t(result))
    result_df$boot_iteration <- i
    result_df$model_name     <- model_name
    result_df$year           <- target_year
    year_results[[length(year_results) + 1]] <- result_df
  }
}

cat(sprintf("Successful iterations: %d / %d\n", length(year_results), n_boot))

# ---------------------------------------------------------------
# 8. Save results
# ---------------------------------------------------------------
if (length(year_results) > 0) {
  year_df  <- bind_rows(year_results)
  out_file <- sprintf("coefs_year%04d.csv", target_year)
  out_path <- file.path(output_dir, out_file)
  write_csv(year_df, out_path)
  cat("Saved results to:", out_path, "\n")
} else {
  cat("WARNING: No successful iterations — no output file written.\n")
}

cat("Done!\n")