#!/usr/bin/env Rscript
# combine_results.R
# Run after all SLURM array tasks finish to merge chunk CSVs into final outputs.
# Usage: Rscript combine_results.R [results_dir]

library(tidyverse)

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "data/results/dengue_validation"

cat("Reading results from:", results_dir, "\n")

result_files <- list.files(results_dir, pattern = "^validation_model.*\\.csv$",
                           full.names = TRUE)

if (length(result_files) == 0) stop("No result files found in: ", results_dir)

cat("Found", length(result_files), "chunk files. Combining...\n")

all_results <- map_dfr(result_files, read_csv, show_col_types = FALSE)

# Save combined file
combined_out <- file.path(results_dir, "all_models_validation_metrics_5year_windows.csv")
write.csv(all_results, combined_out, row.names = FALSE)
cat("Saved:", combined_out, "\n")

# Per-model files
for (mn in unique(all_results$model)) {
  model_results <- filter(all_results, model == mn)
  # Use a filesystem-safe version of the model name
  safe_name <- gsub("[^A-Za-z0-9_]", "_", mn)
  filename  <- file.path(results_dir,
                         paste0("model_validation_", safe_name, "_5year_windows.csv"))
  write.csv(model_results, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
}

# Summary table
summary_results <- all_results %>%
  filter(!is.na(mae)) %>%
  group_by(model) %>%
  summarise(
    n_successful     = n(),
    mean_mae         = mean(mae,         na.rm = TRUE),
    mean_rmse        = mean(rmse,        na.rm = TRUE),
    mean_correlation = mean(correlation, na.rm = TRUE),
    mean_pseudo_r2   = mean(pseudo_r2,   na.rm = TRUE),
    sd_mae           = sd(mae,           na.rm = TRUE),
    sd_rmse          = sd(rmse,          na.rm = TRUE),
    sd_correlation   = sd(correlation,   na.rm = TRUE),
    sd_pseudo_r2     = sd(pseudo_r2,     na.rm = TRUE),
    .groups = 'drop'
  )

summary_out <- file.path(results_dir, "model_validation_summary_5year_windows.csv")
write.csv(summary_results, summary_out, row.names = FALSE)
cat("\nSummary saved:", summary_out, "\n")
print(summary_results)