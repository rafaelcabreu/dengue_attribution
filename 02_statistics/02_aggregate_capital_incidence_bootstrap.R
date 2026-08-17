# aggregate_capital_incidence_bootstrap.R
#
# Run this AFTER all SLURM array tasks from 02_bootstrap_capital_incidence.R
# have completed. Reads all per-chunk CSVs from
# bootstrap_chunks_capital_incidence/, and writes:
#   - capital_incidence_state_blockboot1000.csv    (all draws, long format)
#   - dengue_capital_incidence_boot_summary.csv    (95% CI per city x date)
#
# The summary file adds pred_capital_incidence_q025 / _q500 / _q975 per
# city_residency x date_first_symptoms, for direct use as a confidence band
# in FigureS12.ipynb's municipality time-series panel (merge onto the
# point-estimate data by city_residency + date_first_symptoms).

library(tidyverse)

output_dir <- "data/statistics"
chunks_dir <- "data/statistics/bootstrap_chunks_capital_incidence"
n_chunks   <- 10
chunk_size <- 100
n_boot     <- n_chunks * chunk_size   # 1000

cat("Checking for completed chunk files...\n")
missing <- c()
for (chunk_idx in 1:n_chunks) {
  fname <- sprintf("capital_incidence_chunk%03d.csv", chunk_idx)
  if (!file.exists(file.path(chunks_dir, fname))) missing <- c(missing, chunk_idx)
}
if (length(missing) == 0) {
  cat(sprintf("  All %d chunk files found.\n", n_chunks))
} else {
  cat(sprintf("  MISSING chunks: %s\n", paste(missing, collapse = ", ")))
}
cat("\n")

found_files <- sort(list.files(chunks_dir, pattern = "capital_incidence_chunk\\d{3}\\.csv",
                               full.names = TRUE))

if (length(found_files) == 0) stop("No chunk files found — nothing to aggregate.")

cat("Reading", length(found_files), "chunk files...\n")
boot_results <- map_dfr(found_files, ~read_csv(.x, show_col_types = FALSE))

n_iter_recovered <- n_distinct(boot_results$bootstrap_iteration)
cat("Total bootstrap iterations recovered:", n_iter_recovered, "of", n_boot, "\n")

write_csv(boot_results,
          file.path(output_dir, "capital_incidence_state_blockboot1000.csv"))

boot_summary <- boot_results %>%
  group_by(city_residency, date_first_symptoms) %>%
  summarise(
    pred_capital_incidence_q025 = quantile(pred_capital_incidence, 0.025, na.rm = TRUE),
    pred_capital_incidence_q500 = quantile(pred_capital_incidence, 0.5,   na.rm = TRUE),
    pred_capital_incidence_q975 = quantile(pred_capital_incidence, 0.975, na.rm = TRUE),
    n_valid_iterations          = sum(!is.na(pred_capital_incidence)),
    .groups = "drop"
  ) %>%
  mutate(model_name = 'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity')

write_csv(boot_summary, file.path(output_dir, "dengue_capital_incidence_boot_summary.csv"))

cat("Saved summary:", file.path(output_dir, "dengue_capital_incidence_boot_summary.csv"), "\n")
cat("Median n_valid_iterations per city-date:", median(boot_summary$n_valid_iterations), "\n")
cat("(expect roughly 1000 * (1 - 1/e) =~ 630, since each capital's home\n")
cat(" state is absent from about 1/e of state-block bootstrap draws)\n\n")

cat("Aggregation complete!\n")
