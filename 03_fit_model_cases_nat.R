# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

# load data
dengue_temp <- read_csv("data/model_input_brazil_lag1_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)
  
# Create B-spline basis functions for each temperature lag
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

# Add individual B-spline columns to the dataframe
dengue_temp <- dengue_temp %>%
    mutate(
      # Lag 1 B-splines
      temp_bs_lag11 = temp_bs_fitted_lag1[,1],
      temp_bs_lag12 = temp_bs_fitted_lag1[,2],
      temp_bs_lag13 = temp_bs_fitted_lag1[,3],
      temp_bs_lag14 = temp_bs_fitted_lag1[,4],
      
      # Lag 2 B-splines
      temp_bs_lag21 = temp_bs_fitted_lag2[,1],
      temp_bs_lag22 = temp_bs_fitted_lag2[,2],
      temp_bs_lag23 = temp_bs_fitted_lag2[,3],
      temp_bs_lag24 = temp_bs_fitted_lag2[,4],
      
      # Lag 3 B-splines
      temp_bs_lag31 = temp_bs_fitted_lag3[,1],
      temp_bs_lag32 = temp_bs_fitted_lag3[,2],
      temp_bs_lag33 = temp_bs_fitted_lag3[,3],
      temp_bs_lag34 = temp_bs_fitted_lag3[,4],
      
      # Lag 4 B-splines
      temp_bs_lag41 = temp_bs_fitted_lag4[,1],
      temp_bs_lag42 = temp_bs_fitted_lag4[,2],
      temp_bs_lag43 = temp_bs_fitted_lag4[,3],
      temp_bs_lag44 = temp_bs_fitted_lag4[,4],
      
      # Lag 5 B-splines
      temp_bs_lag51 = temp_bs_fitted_lag5[,1],
      temp_bs_lag52 = temp_bs_fitted_lag5[,2],
      temp_bs_lag53 = temp_bs_fitted_lag5[,3],
      temp_bs_lag54 = temp_bs_fitted_lag5[,4],
    )


main <- fixest::fepois(n_cases ~
  temp_bs_lag11 + temp_bs_lag12 + temp_bs_lag13 + temp_bs_lag14 +
  temp_bs_lag21 + temp_bs_lag22 + temp_bs_lag23 + temp_bs_lag24 +
  temp_bs_lag31 + temp_bs_lag32 + temp_bs_lag33 + temp_bs_lag34 +
  temp_bs_lag41 + temp_bs_lag42 + temp_bs_lag43 + temp_bs_lag44 +
  temp_bs_lag51 + temp_bs_lag52 + temp_bs_lag53 + temp_bs_lag54 +
  total_precipitation_lag1 + total_precipitation_lag2 + total_precipitation_lag3 + 
  total_precipitation_lag4 + total_precipitation_lag5 |
  city_residency + year + month,
  offset = ~log_pop_offset, # population weight
  data = dengue_temp,
  combine.quick = FALSE
)

cat("Base model fitted successfully!\n\n")

# Create list of all files with your specific path pattern
file_list <- paste0("../../../cpdn_nonnerc/aaim/dengue/predict-all/ext-nat-ens", sprintf("%03d", 0:524), ".csv")

# Initialize empty list to store results
all_results <- list()

# Process each file
for(i in 1:length(file_list)) {
  file_path <- file_list[i]
  
  cat("Processing file", i, "of", length(file_list), ":", file_path, "\n")
  
  tryCatch({
    # Load NEW dataset for predictions
    new_data <- read_csv(file_path, show_col_types = FALSE)
    new_data$log_pop_offset <- log(new_data$population/100000)

    # Add individual B-spline columns to the dataframe
    # Inside your loop, replace the bs() calls with predict() calls:
    temp_bs_lag1 <- predict(temp_bs_fitted_lag1, new_data$mean_2m_air_temp_degree1_lag1)
    temp_bs_lag2 <- predict(temp_bs_fitted_lag2, new_data$mean_2m_air_temp_degree1_lag2)
    temp_bs_lag3 <- predict(temp_bs_fitted_lag3, new_data$mean_2m_air_temp_degree1_lag3)
    temp_bs_lag4 <- predict(temp_bs_fitted_lag4, new_data$mean_2m_air_temp_degree1_lag4)
    temp_bs_lag5 <- predict(temp_bs_fitted_lag5, new_data$mean_2m_air_temp_degree1_lag5)      

    new_data <- new_data %>%
      mutate(
        # Lag 1 B-splines
        temp_bs_lag11 = temp_bs_lag1[,1],
        temp_bs_lag12 = temp_bs_lag1[,2],
        temp_bs_lag13 = temp_bs_lag1[,3],
        temp_bs_lag14 = temp_bs_lag1[,4],
        
        # Lag 2 B-splines
        temp_bs_lag21 = temp_bs_lag2[,1],
        temp_bs_lag22 = temp_bs_lag2[,2],
        temp_bs_lag23 = temp_bs_lag2[,3],
        temp_bs_lag24 = temp_bs_lag2[,4],
        
        # Lag 3 B-splines
        temp_bs_lag31 = temp_bs_lag3[,1],
        temp_bs_lag32 = temp_bs_lag3[,2],
        temp_bs_lag33 = temp_bs_lag3[,3],
        temp_bs_lag34 = temp_bs_lag3[,4],
        
        # Lag 4 B-splines
        temp_bs_lag41 = temp_bs_lag4[,1],
        temp_bs_lag42 = temp_bs_lag4[,2],
        temp_bs_lag43 = temp_bs_lag4[,3],
        temp_bs_lag44 = temp_bs_lag4[,4],
        
        # Lag 5 B-splines
        temp_bs_lag51 = temp_bs_lag5[,1],
        temp_bs_lag52 = temp_bs_lag5[,2],
        temp_bs_lag53 = temp_bs_lag5[,3],
        temp_bs_lag54 = temp_bs_lag5[,4],
      )

    # Get predictions using the BASE model (predicted incidence)
    pred_cases <- predict(main, newdata = new_data)
    
    # Convert predicted incidence back to predicted cases
    # Assuming dengue_inc is per 100,000 population
    data_with_predictions <- new_data %>%
      mutate(
        pred_cases = pred_cases,
        pred_incidence = pred_cases / population * 100000  # Convert incidence to cases
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
write_csv(final_results, "nat_dengue_ensemble_predictions.csv")
cat("\nAll results saved to: nat_dengue_ensemble_predictions.csv\n")

# Optional: Save summary of processing
processing_summary <- data.frame(
  file_number = 1:length(file_list),
  file_name = file_list,
  processed_successfully = 1:length(file_list) %in% unique(final_results$file_number)
)

write_csv(processing_summary, "processing_summary.csv")
cat("Processing summary saved to: processing_summary.csv\n")
