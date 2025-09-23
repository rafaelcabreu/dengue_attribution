# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

# load data
dengue_temp <- read_csv("data/model_input_brazil_lag1_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# Load bootstrap state samples
bootstrap_states <- read_csv2("data/bootstrap_state_samples.csv")
cat("Loaded bootstrap states with", nrow(bootstrap_states), "rows\n")

# Parse the bootstrap states - split the resampled_states column
bootstrap_states <- bootstrap_states %>%
  mutate(
    states_list = map(resampled_states, ~str_split(.x, ",")[[1]])
  )

output_dir <- './sensitivity-nat/'

# Define models with temperature terms
models_with_temp <- list(
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

cat("Base data prepared with B-splines!\n\n")

# Create list of all ensemble files
file_list <- paste0("../../../cpdn_nonnerc/aaim/dengue/predict-all/ext-nat-ens", sprintf("%03d", 0:524), ".csv")

# Number of bootstrap iterations
n_bootstrap <- 100 #nrow(bootstrap_states)

# Initialize empty list to store results
all_results <- list()
result_counter <- 1

# Function to create model formula from variable list
create_model_formula <- function(vars) {
  # Separate fixed effects (factor variables) from regular predictors
  fixed_effects <- intersect(vars, c('month', 'year', 'city_residency', 'month_childs', 'year_childs'))
  predictors <- setdiff(vars, fixed_effects)
  
  # Create formula string
  formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "))
  
  # Add fixed effects if any
  if (length(fixed_effects) > 0) {
    formula_str <- paste(formula_str, "|", paste(fixed_effects, collapse = " + "))
  }
  
  return(as.formula(formula_str))
}

# Process each model type
for(model_name in names(models_with_temp)) {
  
  cat("======= PROCESSING MODEL:", model_name, "=======\n")
  model_vars <- models_with_temp[[model_name]]
  cat("Model variables:", paste(model_vars, collapse = ", "), "\n\n")
  
  # Create model formula
  model_formula <- create_model_formula(model_vars)
  cat("Model formula:", deparse(model_formula), "\n\n")
  
  # Collect only this model's results here
  model_results <- list()
  result_counter <- 1
  
  # Process each bootstrap iteration for this model
  for(bootstrap_idx in 1:n_bootstrap) {
    
    cat("=== MODEL:", model_name, "- BOOTSTRAP ITERATION", bootstrap_idx, "of", n_bootstrap, "===\n")
    
    # Get the resampled states for this bootstrap iteration
    resampled_states <- bootstrap_states$states_list[[bootstrap_idx]]
    cat("Resampled states:", paste(resampled_states, collapse = ", "), "\n")

    # Apply the same stratified bootstrap to ensemble data
    ensemble_bootstrap_ids <- data.frame(
      state_residency = resampled_states,
      stringsAsFactors = FALSE
    ) %>%
      mutate(bootstrap_id = 1:n())
    
    tryCatch({
      # Create bootstrap sample of training data
      bootstrap_training_data <- dengue_temp %>%
        left_join(ensemble_bootstrap_ids,
                  by = "state_residency",
                  relationship = "many-to-many") %>%
        filter(!is.na(bootstrap_id))
          
      cat("Bootstrap training data: ", nrow(bootstrap_training_data), "rows\n")
      
      # Fit model on bootstrap sample
      bootstrap_model <- fixest::fepois(
        fml = model_formula,
        offset = ~log_pop_offset,
        data = bootstrap_training_data,
        combine.quick = FALSE
      )
      
      cat("Bootstrap model fitted successfully!\n")
      
      # Apply this bootstrap model to all ensemble members
      for(ensemble_idx in 1:length(file_list)) {
        file_path <- file_list[ensemble_idx]
        
        cat("  Applying to ensemble member", ensemble_idx, "of", length(file_list), "\n")
        
        tryCatch({
          # Load ensemble data
          new_data <- read_csv(file_path, show_col_types = FALSE)
          new_data$log_pop_offset <- log(new_data$population/100000)

          # --- IMPORTANT FIX: use the correct spline objects for precipitation ---
          # Predict temperature splines
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
          
          # Apply stratified bootstrap to ensemble data
          new_data <- new_data %>%
            left_join(ensemble_bootstrap_ids,
                      by = "state_residency",
                      relationship = "many-to-many") %>%
            filter(!is.na(bootstrap_id))

          tlimit <- 12  
          # Filter data to include only observations where ALL lagged temperatures are above 8°C
          new_data <- new_data %>%
            filter(
              mean_2m_air_temp_degree1_lag1 > tlimit &
              mean_2m_air_temp_degree1_lag2 > tlimit &
              mean_2m_air_temp_degree1_lag3 > tlimit &
              mean_2m_air_temp_degree1_lag4 > tlimit &
              mean_2m_air_temp_degree1_lag5 > tlimit
            )
          
          # Get predictions using the bootstrap model
          pred_cases <- predict(bootstrap_model, newdata = new_data)
          
          # Convert predicted incidence back to predicted cases
          data_with_predictions <- new_data %>%
            mutate(
              pred_cases = pred_cases,
              pred_incidence = pred_cases / population * 100000
            )
          
          # Aggregate by region and date_first_symptoms
          aggregated_data <- data_with_predictions %>%
            group_by(region, date_first_symptoms) %>%
            summarise(
              total_population = sum(population, na.rm = TRUE),
              total_pred_cases = sum(pred_cases, na.rm = TRUE),
              total_actual_cases = sum(dengue_inc * population / 100000, na.rm = TRUE),
              .groups = 'drop'
            ) %>%
            mutate(
              # Convert back to incidence rates (per 100,000)
              predicted_incidence = total_pred_cases / total_population * 100000,
              actual_incidence = total_actual_cases / total_population * 100000,
              model_type = model_name,
              bootstrap_iteration = bootstrap_idx,
              ensemble_member = ensemble_idx
            ) %>%
            select(model_type, region, date_first_symptoms, total_population, total_pred_cases,
                   predicted_incidence, actual_incidence, bootstrap_iteration, ensemble_member) %>%
            arrange(region, date_first_symptoms)
          
          # Store result for this model
          model_results[[result_counter]] <- aggregated_data
          result_counter <- result_counter + 1
          
        }, error = function(e) {
          cat("    ERROR processing ensemble member", ensemble_idx, ":", e$message, "\n")
          
          # Store empty result
          model_results[[result_counter]] <- data.frame(
            model_type = model_name,
            region = NA, 
            date_first_symptoms = as.Date(NA), 
            total_population = NA, 
            total_pred_cases = NA,
            predicted_incidence = NA, 
            actual_incidence = NA, 
            bootstrap_iteration = bootstrap_idx,
            ensemble_member = ensemble_idx
          )
          result_counter <- result_counter + 1
        })
      }
      
      cat("Model", model_name, "- Bootstrap iteration", bootstrap_idx, "completed successfully!\n\n")
      
    }, error = function(e) {
      cat("ERROR in model", model_name, "- bootstrap iteration", bootstrap_idx, ":", e$message, "\n")
      
      # Store empty results for all ensemble members of this failed bootstrap
      for(ensemble_idx in 1:length(file_list)) {
        model_results[[result_counter]] <- data.frame(
          model_type = model_name,
          region = NA, 
          date_first_symptoms = as.Date(NA), 
          total_population = NA, 
          total_pred_cases = NA,
          predicted_incidence = NA, 
          actual_incidence = NA, 
          bootstrap_iteration = bootstrap_idx,
          ensemble_member = ensemble_idx
        )
        result_counter <- result_counter + 1
      }
    })
  }
  
  # ---------- Finish this model: bind, filter, save ----------
  cat("Combining results for model:", model_name, "...\n")
  final_model_results <- bind_rows(model_results) %>% filter(!is.na(region))
  
  # Per-model summaries
  summary_stats_model <- final_model_results %>%
    group_by(model_type, bootstrap_iteration, ensemble_member) %>%
    summarise(
      rows_per_combination = n(),
      total_cases = sum(total_pred_cases, na.rm = TRUE),
      n_regions = n_distinct(region),
      .groups = 'drop'
    )
  
  bootstrap_summary_model <- final_model_results %>%
    group_by(model_type, bootstrap_iteration) %>%
    summarise(
      ensemble_members_completed = n_distinct(ensemble_member),
      total_rows = n(),
      successful_rate = ensemble_members_completed / length(file_list),
      .groups = 'drop'
    ) %>%
    mutate(
      resampled_states = map_chr(bootstrap_iteration, ~paste(bootstrap_states$states_list[[.x]], collapse = ", "))
    )
  
  model_summary_model <- tibble(
    model_type = model_name,
    total_combinations = n_distinct(paste(final_model_results$bootstrap_iteration, final_model_results$ensemble_member)),
    expected_combinations = n_bootstrap * length(file_list),
    success_rate = total_combinations / expected_combinations,
    total_rows = nrow(final_model_results)
  )
  
  # Write per-model files
  pred_path <- file.path(output_dir, paste0("predictions_", model_name, ".csv"))
  write_csv(final_model_results, pred_path)
  cat("Saved predictions to:", pred_path, "\n")
  
  summary_path <- file.path(output_dir, paste0("summary_stats_", model_name, ".csv"))
  write_csv(summary_stats_model, summary_path)
  cat("Saved summary stats to:", summary_path, "\n")
  
  bootsum_path <- file.path(output_dir, paste0("bootstrap_processing_summary_", model_name, ".csv"))
  write_csv(bootstrap_summary_model, bootsum_path)
  cat("Saved bootstrap processing summary to:", bootsum_path, "\n")
  
  modelsum_path <- file.path(output_dir, paste0("model_type_summary_", model_name, ".csv"))
  write_csv(model_summary_model, modelsum_path)
  cat("Saved model summary to:", modelsum_path, "\n")
  
  # Free memory before next model
  rm(model_results, final_model_results, summary_stats_model, bootstrap_summary_model, model_summary_model)
  gc()
  
  cat("Model", model_name, "completed and saved!\n\n")
}
