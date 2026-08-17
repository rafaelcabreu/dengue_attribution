# aggregate_rsq_aic_bootstrap.R
#
# Run this AFTER all SLURM array tasks have completed.
# Reads all per-chunk CSVs from bootstrap_chunks_rsq/, checks for
# missing tasks, and writes the same output files as the original
# serial script for each model:
#   <model_name>_rsquared_aic_state_blockboot_presampled1000.csv
#   <model_name>_rsquared_aic_state_blockboot_presampled_summary1000.csv
# Plus a combined file across all models:
#   all_models_bootstrap_rsquared_aic_presampled1000.csv

library(tidyverse)

output_dir  <- "data/statistics"
chunks_dir <- "data/statistics/bootstrap_chunks_rsq"
n_chunks   <- 10
chunk_size <- 100
n_boot     <- n_chunks * chunk_size   # 1000

# ---------------------------------------------------------------
# 1. Model definitions (must match the parallel script exactly)
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

# ---------------------------------------------------------------
# 2. Check which tasks completed
# ---------------------------------------------------------------
cat("Checking for completed chunk files...\n")
for (model_idx in 1:n_models) {
  model_name <- names(models_with_temp)[model_idx]
  missing <- c()
  for (chunk_idx in 1:n_chunks) {
    fname <- sprintf("rsq_model%02d_chunk%03d.csv", model_idx, chunk_idx)
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
# 3. Aggregate per model and write outputs
# ---------------------------------------------------------------
all_bootstrap_results <- list()

for (model_idx in 1:n_models) {
  model_name <- names(models_with_temp)[model_idx]
  cat("======= Aggregating model", model_idx, ":", model_name, "=======\n")

  pattern     <- sprintf("rsq_model%02d_chunk\\d{3}\\.csv", model_idx)
  found_files <- sort(list.files(chunks_dir, pattern = pattern, full.names = TRUE))

  if (length(found_files) == 0) {
    cat("  No chunk files found — skipping.\n\n")
    next
  }

  cat("  Reading", length(found_files), "chunk files...\n")
  boot_results <- map_dfr(found_files, ~read_csv(.x, show_col_types = FALSE))
  boot_results$model_name <- model_name

  cat("  Total bootstrap iterations recovered:", nrow(boot_results), "of", n_boot, "\n")

  # Save combined results — same filename as original serial script
  write_csv(boot_results,
            file.path(output_dir, paste0(model_name,
                   "_rsquared_aic_state_blockboot_presampled", n_boot, ".csv")))

  # Bootstrap summary statistics — same structure as original serial script
  boot_summary <- boot_results %>%
    select(r_squared, pseudo_r_squared, aic, n_obs) %>%
    summarise(
      r_squared_mean         = mean(r_squared,        na.rm = TRUE),
      r_squared_sd           = sd(r_squared,          na.rm = TRUE),
      r_squared_q025         = quantile(r_squared,    0.025, na.rm = TRUE),
      r_squared_q975         = quantile(r_squared,    0.975, na.rm = TRUE),
      pseudo_r_squared_mean  = mean(pseudo_r_squared, na.rm = TRUE),
      pseudo_r_squared_sd    = sd(pseudo_r_squared,   na.rm = TRUE),
      pseudo_r_squared_q025  = quantile(pseudo_r_squared, 0.025, na.rm = TRUE),
      pseudo_r_squared_q975  = quantile(pseudo_r_squared, 0.975, na.rm = TRUE),
      aic_mean               = mean(aic,              na.rm = TRUE),
      aic_sd                 = sd(aic,                na.rm = TRUE),
      aic_q025               = quantile(aic,          0.025, na.rm = TRUE),
      aic_q975               = quantile(aic,          0.975, na.rm = TRUE),
      mean_n_obs             = mean(n_obs,            na.rm = TRUE),
      n_valid_iterations     = sum(!is.na(r_squared))
    ) %>%
    mutate(model_name = model_name)

  write_csv(boot_summary,
            file.path(output_dir, paste0(model_name,
                   "_rsquared_aic_state_blockboot_presampled_summary", n_boot, ".csv")))

  cat("  Mean R-squared:        ", round(boot_summary$r_squared_mean,       4), "\n")
  cat("  Mean Pseudo R-squared: ", round(boot_summary$pseudo_r_squared_mean, 4), "\n")
  cat("  Mean AIC:              ", round(boot_summary$aic_mean,              2), "\n\n")

  all_bootstrap_results[[model_name]] <- boot_results
}

# ---------------------------------------------------------------
# 4. Write combined output across all models
# ---------------------------------------------------------------
if (length(all_bootstrap_results) > 0) {
  combined_bootstrap <- bind_rows(all_bootstrap_results)
  out_file <- paste0("all_models_bootstrap_rsquared_aic_presampled", n_boot, ".csv")
  write_csv(combined_bootstrap, file.path(output_dir, out_file))
  cat("Combined results saved to:", out_file, "\n")
}

cat("\nAggregation complete!\n")