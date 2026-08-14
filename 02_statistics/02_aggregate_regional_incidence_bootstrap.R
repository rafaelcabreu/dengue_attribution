# aggregate_regional_incidence_bootstrap.R
#
# Run this AFTER all SLURM array tasks from 02_bootstrap_regional_incidence.R
# have completed. Reads all per-chunk CSVs from
# bootstrap_chunks_regional_incidence/, and for each model writes:
#   - <model>_regional_incidence_state_blockboot1000.csv   (all draws, long format)
#   - dengue_regional_aggregation_boot_summary.csv          (fixed effects model, 95% CI)
#   - dengue_regional_aggregation_nofixed_boot_summary.csv  (no fixed effects model, 95% CI)
#
# The summary files add pred_regional_incidence_q025 / _q500 / _q975 per
# region x date_first_symptoms, for direct use as confidence bands in
# Figure2.ipynb (merge onto dengue_regional_aggregation.csv /
# dengue_regional_aggregation_nofixed.csv by region + date_first_symptoms).

library(tidyverse)

output_dir <- "data/statistics"
chunks_dir <- "data/statistics/bootstrap_chunks_regional_incidence"
n_chunks   <- 10
chunk_size <- 100
n_boot     <- n_chunks * chunk_size   # 1000

model_names <- c(
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity',
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity'
)
model_output_slugs <- c(
  "dengue_regional_aggregation_boot_summary.csv",
  "dengue_regional_aggregation_nofixed_boot_summary.csv"
)

n_models <- length(model_names)

cat("Checking for completed chunk files...\n")
for (model_idx in 1:n_models) {
  missing <- c()
  for (chunk_idx in 1:n_chunks) {
    fname <- sprintf("regional_incidence_model%02d_chunk%03d.csv", model_idx, chunk_idx)
    if (!file.exists(file.path(chunks_dir, fname))) missing <- c(missing, chunk_idx)
  }
  if (length(missing) == 0) {
    cat(sprintf("  Model %d (%s): all %d chunk files found.\n",
                model_idx, model_names[model_idx], n_chunks))
  } else {
    cat(sprintf("  Model %d (%s): MISSING chunks: %s\n",
                model_idx, model_names[model_idx], paste(missing, collapse = ", ")))
  }
}
cat("\n")

for (model_idx in 1:n_models) {
  model_name <- model_names[model_idx]
  cat("======= Aggregating model", model_idx, ":", model_name, "=======\n")

  pattern     <- sprintf("regional_incidence_model%02d_chunk\\d{3}\\.csv", model_idx)
  found_files <- sort(list.files(chunks_dir, pattern = pattern, full.names = TRUE))

  if (length(found_files) == 0) {
    cat("  No chunk files found — skipping.\n\n")
    next
  }

  cat("  Reading", length(found_files), "chunk files...\n")
  boot_results <- map_dfr(found_files, ~read_csv(.x, show_col_types = FALSE))

  n_iter_recovered <- n_distinct(boot_results$bootstrap_iteration)
  cat("  Total bootstrap iterations recovered:", n_iter_recovered, "of", n_boot, "\n")

  write_csv(boot_results,
            file.path(output_dir, paste0(model_name,
                   "_regional_incidence_state_blockboot", n_boot, ".csv")))

  boot_summary <- boot_results %>%
    group_by(region, date_first_symptoms) %>%
    summarise(
      pred_regional_incidence_q025 = quantile(pred_regional_incidence, 0.025, na.rm = TRUE),
      pred_regional_incidence_q500 = quantile(pred_regional_incidence, 0.5,   na.rm = TRUE),
      pred_regional_incidence_q975 = quantile(pred_regional_incidence, 0.975, na.rm = TRUE),
      n_valid_iterations           = sum(!is.na(pred_regional_incidence)),
      .groups = "drop"
    ) %>%
    mutate(model_name = model_name)

  write_csv(boot_summary, file.path(output_dir, model_output_slugs[model_idx]))

  cat("  Saved summary:", file.path(output_dir, model_output_slugs[model_idx]), "\n")
  cat("  Median n_valid_iterations per region-date:",
      median(boot_summary$n_valid_iterations), "\n\n")
}

cat("Aggregation complete!\n")
