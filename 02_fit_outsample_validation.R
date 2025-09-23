# Load required libraries
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

# Load and prepare initial data
dengue_temp <- read_csv("model_input_brazil_lag1_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# Convert date column to Date type if not already
dengue_temp$date_first_symptoms <- as.Date(dengue_temp$date_first_symptoms)

# Define models with temperature terms
models_with_temp <- list(
  'base' = c(
    'month',
    'year',
    'city_residency'
  ),
  'standard' = c(
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
    'month',
    'year',
    'city_residency'
  ),
  'with_immunity' = c(
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
    'month',
    'year',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3'
  ),
  'lag_45' = c(
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
    'month',
    'year',
    'city_residency'
  ),
  'childs' = c(
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
    'month_childs',
    'year_childs',
    'city_residency'
  ),
  'lag1' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',
    'total_precipitation_lag1',
    'month',
    'year',
    'city_residency'
  ),
  'socioeconomic' = c(
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
    'month',
    'year',
    'city_residency',
    'water_sanitation',
    'pib_2021',
    'percent_urban',
    'urban_pop'
  ),
  'lag1_cases' = c(
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
    'n_cases_lag1',
    'month',
    'year',
    'city_residency'
  )
)

# Function to build formula for each model
build_formula <- function(model_vars) {
  # Separate fixed effects from regular variables
  fixed_effects <- c('city_residency', 'year', 'month', 'year_childs', 'month_childs')
  fe_vars <- intersect(model_vars, fixed_effects)
  regular_vars <- setdiff(model_vars, fixed_effects)
  
  # Build formula parts
  if(length(regular_vars) > 0) {
    rhs <- paste(regular_vars, collapse = " + ")
  } else {
    rhs <- "1"
  }
  
  if(length(fe_vars) > 0) {
    fe_part <- paste(fe_vars, collapse = " + ")
    formula_str <- paste("n_cases ~", rhs, "|", fe_part)
  } else {
    formula_str <- paste("n_cases ~", rhs)
  }
  
  return(as.formula(formula_str))
}

# Function to calculate McFadden's pseudo R-squared
calculate_pseudo_r2 <- function(model, validation_data) {
  tryCatch({
    # Get predictions from the full model
    predictions <- predict(model, newdata = validation_data)
    
    # Calculate log-likelihood for the full model
    observed <- validation_data$n_cases
    ll_full <- sum(dpois(observed, predictions, log = TRUE), na.rm = TRUE)
    
    # Calculate log-likelihood for null model (intercept only)
    # For null model, use the mean of observed values as prediction
    null_prediction <- mean(observed, na.rm = TRUE)
    ll_null <- sum(dpois(observed, null_prediction, log = TRUE), na.rm = TRUE)
    
    # McFadden's pseudo R-squared
    pseudo_r2 <- 1 - (ll_full / ll_null)
    
    return(pseudo_r2)
    
  }, error = function(e) {
    return(NA)
  })
}

# Function to run one iteration for a specific model
run_iteration <- function(iteration, model_name, model_vars) {
  # Set different seed for each iteration
  set.seed(123 + iteration)
  
  # Select 20% of unique cities randomly
  selected_cities <- sample(unique(dengue_temp$city_residency), 
                           round(length(unique(dengue_temp$city_residency)) * 0.2))
  
  # Get available years for selected cities
  selected_cities_data <- dengue_temp[dengue_temp$city_residency %in% selected_cities, ]
  available_years <- sort(unique(year(selected_cities_data$date_first_symptoms)))
  
  # Check if we have at least 10 years of data (5 for training + 5 for validation)
  if(length(available_years) < 10) {
    return(data.frame(model = model_name, iteration = iteration, 
                     mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
                     note = "Insufficient years of data"))
  }
  
  # Find the latest possible starting year that allows for 10 consecutive years
  max_start_year <- available_years[length(available_years) - 9]
  min_start_year <- available_years[1]
  
  # Randomly select a starting year
  possible_start_years <- available_years[available_years >= min_start_year & 
                                        available_years <= max_start_year]
  
  if(length(possible_start_years) == 0) {
    return(data.frame(model = model_name, iteration = iteration, 
                     mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
                     note = "No valid start years"))
  }
  
  start_year <- sample(possible_start_years, 1)
  training_years <- start_year:(start_year + 4)
  validation_years <- (start_year + 5):(start_year + 9)
  
  # Create training data: selected cities for training years + all other cities
  training_selected <- selected_cities_data[
    year(selected_cities_data$date_first_symptoms) %in% training_years, ]
  
  # All data from non-selected cities
  other_cities_data <- dengue_temp[!dengue_temp$city_residency %in% selected_cities, ]
  
  # Combine training data
  training_data <- rbind(other_cities_data, training_selected)
  
  # Create validation data: selected cities for validation years
  validation_data <- selected_cities_data[
    year(selected_cities_data$date_first_symptoms) %in% validation_years, ]
  
  # Skip if no validation data
  if(nrow(validation_data) == 0) {
    return(data.frame(model = model_name, iteration = iteration, 
                     mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
                     note = "No validation data"))
  }
  
  # Create B-spline basis functions for training data (only if needed)
  if(any(grepl("temp_bs_lag", model_vars))) {
    temp_bs_lag1 <- bs(training_data$mean_2m_air_temp_degree1_lag1, df = 4)
    temp_bs_lag2 <- bs(training_data$mean_2m_air_temp_degree1_lag2, df = 4)
    temp_bs_lag3 <- bs(training_data$mean_2m_air_temp_degree1_lag3, df = 4)
    temp_bs_lag4 <- bs(training_data$mean_2m_air_temp_degree1_lag4, df = 4)
    temp_bs_lag5 <- bs(training_data$mean_2m_air_temp_degree1_lag5, df = 4)
    
    # Add B-spline columns to training data
    training_data <- training_data %>%
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
  }
  
  # Build and fit the model
  model_formula <- build_formula(model_vars)
  
  tryCatch({
    main <- fixest::fepois(model_formula,
                          offset = ~log_pop_offset,
                          data = training_data,
                          combine.quick = FALSE)
    
    # Create B-splines for validation data (only if needed)
    if(any(grepl("temp_bs_lag", model_vars))) {
      temp_bs_lag1_validation <- predict(temp_bs_lag1, validation_data$mean_2m_air_temp_degree1_lag1)
      temp_bs_lag2_validation <- predict(temp_bs_lag2, validation_data$mean_2m_air_temp_degree1_lag2)
      temp_bs_lag3_validation <- predict(temp_bs_lag3, validation_data$mean_2m_air_temp_degree1_lag3)
      temp_bs_lag4_validation <- predict(temp_bs_lag4, validation_data$mean_2m_air_temp_degree1_lag4)
      temp_bs_lag5_validation <- predict(temp_bs_lag5, validation_data$mean_2m_air_temp_degree1_lag5)
  
      # Add B-spline columns to validation data
      validation_data <- validation_data %>%
        mutate(
          # Lag 1 B-splines
          temp_bs_lag11 = temp_bs_lag1_validation[,1],
          temp_bs_lag12 = temp_bs_lag1_validation[,2],
          temp_bs_lag13 = temp_bs_lag1_validation[,3],
          temp_bs_lag14 = temp_bs_lag1_validation[,4],
          # Lag 2 B-splines
          temp_bs_lag21 = temp_bs_lag2_validation[,1],
          temp_bs_lag22 = temp_bs_lag2_validation[,2],
          temp_bs_lag23 = temp_bs_lag2_validation[,3],
          temp_bs_lag24 = temp_bs_lag2_validation[,4],
          # Lag 3 B-splines
          temp_bs_lag31 = temp_bs_lag3_validation[,1],
          temp_bs_lag32 = temp_bs_lag3_validation[,2],
          temp_bs_lag33 = temp_bs_lag3_validation[,3],
          temp_bs_lag34 = temp_bs_lag3_validation[,4],
          # Lag 4 B-splines
          temp_bs_lag41 = temp_bs_lag4_validation[,1],
          temp_bs_lag42 = temp_bs_lag4_validation[,2],
          temp_bs_lag43 = temp_bs_lag4_validation[,3],
          temp_bs_lag44 = temp_bs_lag4_validation[,4],
          # Lag 5 B-splines
          temp_bs_lag51 = temp_bs_lag5_validation[,1],
          temp_bs_lag52 = temp_bs_lag5_validation[,2],
          temp_bs_lag53 = temp_bs_lag5_validation[,3],
          temp_bs_lag54 = temp_bs_lag5_validation[,4],
        )
    }
    
    # Make predictions
    predictions <- predict(main, newdata = validation_data)
    
    # Calculate metrics
    observed <- validation_data$n_cases
    predicted <- predictions
    
    mae <- mean(abs(observed - predicted), na.rm = TRUE)
    rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
    correlation <- cor(observed, predicted, use = "complete.obs")
    
    # Calculate pseudo R-squared
    pseudo_r2 <- calculate_pseudo_r2(main, validation_data)
    
    return(data.frame(model = model_name, iteration = iteration, 
                     mae = mae, rmse = rmse, correlation = correlation,
                     pseudo_r2 = pseudo_r2,
                     training_years = paste(training_years, collapse = "-"),
                     validation_years = paste(validation_years, collapse = "-"),
                     note = "Success"))
    
  }, error = function(e) {
    return(data.frame(model = model_name, iteration = iteration, 
                     mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
                     training_years = NA, validation_years = NA,
                     note = paste("Error:", e$message)))
  })
}

# Run 100 iterations for each model and save individual files
all_results <- data.frame()

for(model_name in names(models_with_temp)) {
  cat("Running model:", model_name, "\n")
  model_vars <- models_with_temp[[model_name]]
  
  # Initialize results for this model
  model_results <- data.frame()
  
  for(i in 1:100) {
    if(i %% 20 == 0) cat("  Completed", i, "iterations for", model_name, "\n")
    
    result <- run_iteration(i, model_name, model_vars)
    model_results <- rbind(model_results, result)
    all_results <- rbind(all_results, result)
  }
  
  # Save individual model results to separate file
  filename <- paste0("model_validation_", model_name, "_5year_windows.csv")
  write.csv(model_results, filename, row.names = FALSE)
  cat("Saved results for", model_name, "to", filename, "\n")
  
  # Calculate and display summary for this model
  model_summary <- model_results %>%
    filter(!is.na(mae)) %>%
    summarise(
      n_successful = n(),
      mean_mae = mean(mae, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_correlation = mean(correlation, na.rm = TRUE),
      mean_pseudo_r2 = mean(pseudo_r2, na.rm = TRUE),
      sd_mae = sd(mae, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      sd_pseudo_r2 = sd(pseudo_r2, na.rm = TRUE)
    )
  
  cat("Summary for", model_name, ":\n")
  print(model_summary)
  
  cat("Finished model:", model_name, "\n\n")
}

# Save combined results (original functionality maintained)
write.csv(all_results, "all_models_validation_metrics_5year_windows.csv", row.names = FALSE)
cat("Combined results saved to 'all_models_validation_metrics_5year_windows.csv'\n")

# Summary by model (excluding failed iterations)
cat("Overall summary of results by model:\n")
summary_results <- all_results %>%
  filter(!is.na(mae)) %>%
  group_by(model) %>%
  summarise(
    n_successful = n(),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_correlation = mean(correlation, na.rm = TRUE),
    mean_pseudo_r2 = mean(pseudo_r2, na.rm = TRUE),
    sd_mae = sd(mae, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    sd_correlation = sd(correlation, na.rm = TRUE),
    sd_pseudo_r2 = sd(pseudo_r2, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_results)

# Save summary
write.csv(summary_results, "model_validation_summary_5year_windows.csv", row.names = FALSE)

# Show details of failed iterations if any
failed_iterations <- all_results %>%
  filter(is.na(mae)) %>%
  count(model, note)

if(nrow(failed_iterations) > 0) {
  cat("\nFailed iterations by model and reason:\n")
  print(failed_iterations)
}

# Show year range information
cat("\nData year range information:\n")
actual_years <- sort(unique(year(dengue_temp$date_first_symptoms)))
cat("Actual years in dataset:", min(actual_years), "to", max(actual_years), "\n")
cat("Number of years available:", length(actual_years), "\n")

# Show sample of successful year combinations
successful_combinations <- all_results %>%
  filter(!is.na(mae)) %>%
  select(model, training_years, validation_years) %>%
  distinct() %>%
  head(10)

if(nrow(successful_combinations) > 0) {
  cat("\nSample of successful year combinations (first 10):\n")
  print(successful_combinations)
}