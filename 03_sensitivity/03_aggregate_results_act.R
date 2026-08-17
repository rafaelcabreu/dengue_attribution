# aggregate_results_regional.R
#
# Run this AFTER all SLURM array tasks have completed.
# Reads all per-task CSVs from sensitivity-act/bootstrap_chunks_regional/,
# checks for missing tasks, and writes the same four summary output
# files per model as the original serial script.

library(tidyverse)

output_dir  <- "data/sensitivity-bspsprec-act"
chunks_dir  <- file.path(output_dir, "bootstrap_chunks")
n_bootstrap <- 100
n_ensemble  <- 525   # ensemble members 0:524

models_with_temp <- list(
  # 'Climate(lag 1-5)',
  # 'Climate(lag 1-5) + Year|region',
  # 'Climate(lag 1-5) + Month|region',
  # 'Year + Month + P(Climate)',
  # 'Climate(lag 1) + Year|region + Month|region',
  # 'Climate(lag 1-3) + Year|region + Month|region',
  'Climate(lag 1-5) + Year|region + Month|region'
  # 'Climate(lag 1-5) + Year|region + Month|region + Socio',
  # 'Climate(lag 1-5) + Year|region + Month|region + PriorCases',
  # 'Climate(lag 1-5) + Year|region + Month|region + SeroRepla',
  # 'Climate(lag 1-5) + Year|region + Month|region + Immunity',
  # 'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity',
  # 'Climate(lag 1-5) + Socio',
  # 'Climate(lag 1-5) + PriorCases',
  # 'Climate(lag 1-5) + SeroRepla',
  # 'Climate(lag 1-5) + Immunity',
  # 'Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity',
  # 'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio'
)
n_models <- length(models_with_temp)

# ---------------------------------------------------------------
# 1. Check which tasks completed
# ---------------------------------------------------------------
cat("Checking for completed chunk files...\n")
for (model_idx in 1:n_models) {
  model_name <- models_with_temp[[model_idx]]
  missing <- c()
  for (bootstrap_idx in 1:n_bootstrap) {
    fname <- sprintf("predictions_model%02d_bootstrap%04d.csv", model_idx, bootstrap_idx)
    if (!file.exists(file.path(chunks_dir, fname))) missing <- c(missing, bootstrap_idx)
  }
  if (length(missing) == 0) {
    cat(sprintf("  Model %d (%s): all %d bootstrap files found.\n",
                model_idx, model_name, n_bootstrap))
  } else {
    cat(sprintf("  Model %d (%s): MISSING bootstrap iterations: %s\n",
                model_idx, model_name, paste(missing, collapse = ", ")))
  }
}
cat("\n")

# ---------------------------------------------------------------
# 2. Load, combine, summarise, and save per model
# ---------------------------------------------------------------
for (model_idx in 1:n_models) {
  model_name <- models_with_temp[[model_idx]]
  cat("======= Aggregating model", model_idx, ":", model_name, "=======\n")
  
  # Find all chunk files for this model
  pattern     <- sprintf("predictions_model%02d_bootstrap\\d{4}\\.csv", model_idx)
  found_files <- sort(list.files(chunks_dir, pattern = pattern, full.names = TRUE))
  
  if (length(found_files) == 0) {
    cat("  No files found — skipping.\n\n")
    next
  }
  
  cat("  Reading", length(found_files), "chunk files...\n")
  final_model_results <- map_dfr(found_files, ~read_csv(.x, show_col_types = FALSE)) %>%
    filter(!is.na(region))
  
  cat("  Total rows:", nrow(final_model_results), "\n")
  
  # Per bootstrap x ensemble summary
  summary_stats_model <- final_model_results %>%
    group_by(model_type, bootstrap_iteration, ensemble_member) %>%
    summarise(
      rows_per_combination = n(),
      total_cases          = sum(total_pred_cases, na.rm = TRUE),
      n_regions            = n_distinct(region),
      .groups = 'drop'
    )
  
  # Per bootstrap summary
  bootstrap_summary_model <- final_model_results %>%
    group_by(model_type, bootstrap_iteration) %>%
    summarise(
      ensemble_members_completed = n_distinct(ensemble_member),
      total_rows                 = n(),
      successful_rate            = ensemble_members_completed / n_ensemble,
      .groups = 'drop'
    )
  
  # Overall model summary
  model_summary_model <- tibble(
    model_type            = model_name,
    total_combinations    = n_distinct(paste(final_model_results$bootstrap_iteration,
                                             final_model_results$ensemble_member)),
    expected_combinations = n_bootstrap * n_ensemble,
    success_rate          = total_combinations / expected_combinations,
    total_rows            = nrow(final_model_results)
  )
  
  # Write outputs — filenames match the original serial script with _regional suffix
  write_csv(final_model_results,
            file.path(output_dir, paste0("predictions_regional_", model_name, ".csv")))
  write_csv(summary_stats_model,
            file.path(output_dir, paste0("summary_stats_regional_", model_name, ".csv")))
  write_csv(bootstrap_summary_model,
            file.path(output_dir, paste0("bootstrap_processing_summary_regional_", model_name, ".csv")))
  write_csv(model_summary_model,
            file.path(output_dir, paste0("model_type_summary_regional_", model_name, ".csv")))
  
  cat("  Saved all outputs for model:", model_name, "\n")
  print(model_summary_model)
  cat("\n")
  
  rm(final_model_results, summary_stats_model, bootstrap_summary_model, model_summary_model)
  gc()
}

cat("Aggregation complete!\n")