# load packages
library(tidyverse)
library(magrittr)
library(fixest)

dengue_temp <- read_csv("data/model_input_brazil_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# Fit your original model
main <- fixest::fepois(n_cases ~
  mean_2m_air_temp_degree1_lag1 + mean_2m_air_temp_degree2_lag1 + 
  mean_2m_air_temp_degree1_lag2 + mean_2m_air_temp_degree2_lag2 + 
  mean_2m_air_temp_degree1_lag3 + mean_2m_air_temp_degree2_lag3 + 
  mean_2m_air_temp_degree1_lag4 + mean_2m_air_temp_degree2_lag4 + 
  mean_2m_air_temp_degree1_lag5 + mean_2m_air_temp_degree2_lag5 +
  total_precipitation_lag1 + total_precipitation_lag2 + total_precipitation_lag3 + 
  total_precipitation_lag4 + total_precipitation_lag5 |
  city_residency + year + month,
  offset = ~log_pop_offset, # population weight
  data = dengue_temp,
  combine.quick = FALSE
)

print(summary(main))

# Get fitted values
fitted_log <- fitted(main)
fitted_incidence <- exp(fitted_log)

# Get predictions for all data (handles dropped observations)
pred_incidence <- predict(main, newdata = dengue_temp)

# Create dataset with fitted values
data_with_fitted <- dengue_temp %>%
  mutate(
    pred_incidence = pred_incidence,
    pred_cases = pred_incidence,  # Convert incidence to cases
    actual_cases = n_cases        # Convert actual incidence to cases
  )

# File 1: Save fitted values dataset
write_csv(data_with_fitted, "data/dengue_with_fitted_values.csv")
cat("File 1 saved: data/dengue_with_fitted_values.csv\n")
cat("Contains original data plus fitted incidence and cases\n")

# File 2: Create regional aggregation
regional_aggregation <- data_with_fitted %>%
  group_by(region, date_first_symptoms) %>%
  summarise(
    # Sum cases and population
    total_actual_cases = sum(actual_cases, na.rm = TRUE),
    total_pred_cases = sum(pred_cases, na.rm = TRUE),
    total_population = sum(population, na.rm = TRUE),
    
    # Additional useful metrics
    n_observations = n(),
    n_cities = length(unique(city_residency)),
    
    .groups = 'drop'
  ) %>%
  mutate(
    # Calculate regional incidence rates (per 100,000)
    actual_regional_incidence = total_actual_cases / total_population * 100000,
    pred_regional_incidence = total_pred_cases / total_population * 100000,
    
    # Calculate absolute and relative errors
    absolute_error = abs(actual_regional_incidence - pred_regional_incidence),
    relative_error = absolute_error / actual_regional_incidence * 100
  ) %>%
  arrange(region, date_first_symptoms)

# Save regional aggregation
write_csv(regional_aggregation, "data/dengue_regional_aggregation.csv")
cat("File 2 saved: data/dengue_regional_aggregation.csv\n")
cat("Contains regional incidence rates and prediction errors\n")

# Display summary statistics
cat("\n=== FITTED VALUES SUMMARY ===\n")
cat("Total observations with fitted values:", nrow(data_with_fitted), "\n")
cat("Mean actual incidence:", round(mean(data_with_fitted$dengue_inc, na.rm = TRUE), 3), "\n")
cat("Mean predicted incidence:", round(mean(data_with_fitted$pred_incidence, na.rm = TRUE), 3), "\n")

cat("\n=== REGIONAL AGGREGATION SUMMARY ===\n")
cat("Number of region-date combinations:", nrow(regional_aggregation), "\n")
cat("Number of unique regions:", length(unique(regional_aggregation$region)), "\n")
cat("Date range:", min(regional_aggregation$date_first_symptoms), "to", max(regional_aggregation$date_first_symptoms), "\n")

# Regional performance metrics
regional_stats <- regional_aggregation %>%
  filter(!is.na(actual_regional_incidence) & !is.na(pred_regional_incidence)) %>%
  summarise(
    correlation = cor(actual_regional_incidence, pred_regional_incidence, use = "complete.obs"),
    r_squared = cor(actual_regional_incidence, pred_regional_incidence, use = "complete.obs")^2,
    rmse = sqrt(mean((actual_regional_incidence - pred_regional_incidence)^2, na.rm = TRUE)),
    mae = mean(abs(actual_regional_incidence - pred_regional_incidence), na.rm = TRUE),
    mape = mean(abs((actual_regional_incidence - pred_regional_incidence) / actual_regional_incidence) * 100, na.rm = TRUE)
  )

cat("\nRegional Performance Metrics:\n")
cat("Correlation:", round(regional_stats$correlation, 4), "\n")
cat("R-squared:", round(regional_stats$r_squared, 4), "\n")
cat("RMSE:", round(regional_stats$rmse, 3), "\n")
cat("MAE:", round(regional_stats$mae, 3), "\n")
cat("MAPE:", round(regional_stats$mape, 2), "%\n")

# Show sample of regional data
cat("\n=== SAMPLE REGIONAL DATA ===\n")
print(head(regional_aggregation %>% 
            select(region, date_first_symptoms, actual_regional_incidence, pred_regional_incidence, absolute_error), 10))

# Optional: Extract fixed effects
fixedEffects <- fixest::fixef(main, type = "detailed")
cat("\nFixed effects extracted (use print(fixedEffects) to view)\n")

cat("\n=== FILES CREATED ===\n")
cat("1. dengue_with_fitted_values.csv - Individual level with fitted values\n")
cat("2. dengue_regional_aggregation.csv - Regional aggregation with incidence rates\n")