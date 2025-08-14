# load packages
library(tidyverse)
library(magrittr)
library(fixest)

# Load TRAINING data and fit the BASE model (only once)
cat("Loading training data and fitting base model...\n")
dengue_temp <- read_csv("model_input_brazil_immunity_city.csv", show_col_types = FALSE)

# Fit your BASE model
main <- fixest::fepois(dengue_inc ~
  mean_2m_air_temp_degree1_lag1 + mean_2m_air_temp_degree2_lag1 + mean_2m_air_temp_degree3_lag1 +
  mean_2m_air_temp_degree1_lag2 + mean_2m_air_temp_degree2_lag2 + mean_2m_air_temp_degree3_lag2 +
  mean_2m_air_temp_degree1_lag3 + mean_2m_air_temp_degree2_lag3 + mean_2m_air_temp_degree3_lag3 +
  mean_2m_air_temp_degree1_lag4 + mean_2m_air_temp_degree2_lag4 + mean_2m_air_temp_degree3_lag4 +
  mean_2m_air_temp_degree1_lag5 + mean_2m_air_temp_degree2_lag5 + mean_2m_air_temp_degree3_lag5 +
  total_precipitation_lag1 + total_precipitation_lag2 + total_precipitation_lag3 + 
  total_precipitation_lag4 + total_precipitation_lag5 |
  city_residency + year + month,
  weights = ~population,
  data = dengue_temp,
  combine.quick = FALSE
)

cat("Base model fitted successfully!\n\n")

# Create list of all files with your specific path pattern
file_list <- paste0("../../cpdn_nonnerc/aaim/dengue/predict-all/ext-ens", sprintf("%03d", 0:524), ".csv")

# Initialize empty list to store results
all_results <- list()

# Process each file
for(i in 1:length(file_list)) {
  file_path <- file_list[i]
  
  cat("Processing file", i, "of", length(file_list), ":", file_path, "\n")
  
  tryCatch({
    # Load NEW dataset for predictions
    new_data <- read_csv(file_path, show_col_types = FALSE)
    
    # Get predictions using the BASE model (predicted incidence)
    pred_incidence <- predict(main, newdata = new_data)
    
    # Convert predicted incidence back to predicted cases
    # Assuming dengue_inc is per 100,000 population
    data_with_predictions <- new_data %>%
      mutate(
        pred_incidence = pred_incidence,
        pred_cases = pred_incidence * population / 100000  # Convert incidence to cases
      )
    
    # Aggregate by region and date_first_symptoms
    aggregated_data <- data_with_predictions %>%
      group_by(region, date_first_symptoms) %>%
      summarise(
        total_population = sum(population, na.rm = TRUE),
        total_pred_cases = sum(pred_cases, na.rm = TRUE),
        total_actual_cases = sum(dengue_inc * population / 100000, na.rm = TRUE),  # Convert actual incidence to cases
        .groups = 'drop'
      ) %>%
      mutate(
        # Convert back to incidence rates (per 100,000)
        predicted_incidence = total_pred_cases / total_population * 100000,
        actual_incidence = total_actual_cases / total_population * 100000,
        file_number = i
      ) %>%
      select(region, date_first_symptoms, total_population, predicted_incidence, actual_incidence, file_number) %>%
      arrange(region, date_first_symptoms)
    
    # Store result
    all_results[[i]] <- aggregated_data
    
    cat("File", i, "processed successfully -", nrow(aggregated_data), "rows\n")
    
  }, error = function(e) {
    cat("ERROR processing file", i, ":", e$message, "\n")
    # Store empty result with file number for tracking
    all_results[[i]] <- data.frame(
      region = NA, date_first_symptoms = as.Date(NA), 
      total_population = NA, predicted_incidence = NA, actual_incidence = NA, 
      file_number = i
    )
  })
}

# Combine all results into single dataframe
cat("\nCombining all results...\n")
final_results <- bind_rows(all_results)

# Remove rows where all values are NA (failed files)
final_results <- final_results %>%
  filter(!is.na(region))

# Display summary
cat("\n=== FINAL RESULTS SUMMARY ===\n")
cat("Total files processed:", length(file_list), "\n")
cat("Total rows in final dataset:", nrow(final_results), "\n")
cat("Successful files:", length(unique(final_results$file_number)), "\n")
cat("Failed files:", length(file_list) - length(unique(final_results$file_number)), "\n")

# Show summary by file
files_summary <- final_results %>%
  group_by(file_number) %>%
  summarise(rows_per_file = n(), .groups = 'drop')

cat("Rows per file summary:\n")
print(summary(files_summary$rows_per_file))

# Save final combined results
write_csv(final_results, "all_dengue_ensemble_predictions.csv")
cat("\nAll results saved to: all_dengue_ensemble_predictions.csv\n")

# Optional: Save summary of processing
processing_summary <- data.frame(
  file_number = 1:length(file_list),
  file_name = file_list,
  processed_successfully = 1:length(file_list) %in% unique(final_results$file_number)
)

write_csv(processing_summary, "processing_summary.csv")
cat("Processing summary saved to: processing_summary.csv\n")
