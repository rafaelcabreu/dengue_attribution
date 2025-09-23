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

# stratified bootstrap function for state-level clustering ----
boot_strat_newID <- function(df_ids, # dataset with unique state IDs
                             df_full, # full dataset
                             state_id_var = "state_residency",
                             seed = 1234){
  set.seed(seed)
  # Sample states with replacement
  ids <- df_ids %>%
    slice_sample(prop = 1, replace = T) %>%
    mutate(state_boot_id = 1:n())

  # Join back to full dataset and create new cluster IDs
  df_out <- df_full %>%
    left_join(ids,
              by = state_id_var,
              relationship = "many-to-many") %>%
    filter(!is.na(state_boot_id))

  return(df_out)
}

# fitting function for bootstrapped samples ----
boot_fit_model <- function(df_ids, # dataset with state IDs
                           df_full, # full dataset
                           model_vars, # variables for this model
                           state_id_var = "state_residency",
                           seed = 1234){

  # Get bootstrapped sample
  boot_data <- boot_strat_newID(df_ids, df_full, state_id_var, seed)

  # Separate fixed effects from regular predictors
  fixed_effects <- c()
  predictors <- c()
  
  for(var in model_vars) {
    if(var %in% c("city_residency", "month", "year", "month_childs", "year_childs")) {
      fixed_effects <- c(fixed_effects, var)
    } else {
      predictors <- c(predictors, var)
    }
  }
  
  # Create formula
  if(length(predictors) > 0 && length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "), "|", paste(fixed_effects, collapse = " + "))
  } else if(length(predictors) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "))
  } else if(length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~ 1 |", paste(fixed_effects, collapse = " + "))
  } else {
    formula_str <- "n_cases ~ 1"
  }
  
  model_formula <- as.formula(formula_str)

  # Fit the model with new cluster structure
  model <- fixest::fepois(model_formula,
                          offset = ~log_pop_offset, # population weight
                          data = boot_data,
                          nthreads = 1)  # Changed to single thread

  return(coef(model))
}

# Create dataset of unique states
df_states <- dengue_temp %>%
  select(state_residency) %>%
  unique() %>%
  filter(!is.na(state_residency))  # Remove any NA values

print(paste("Number of unique states:", nrow(df_states)))
print("Unique states:")
print(df_states$state_residency)

# Set number of bootstrap iterations
n_boot <- 1000  # Adjust as needed

# Storage for all results
all_bootstrap_results <- list()
all_original_results <- list()

# Loop through each model
for(model_name in names(models_with_temp)) {
  cat("\n=== Running bootstrap for model:", model_name, "===\n")
  
  # Get variables for this model
  model_vars <- models_with_temp[[model_name]]
  
  cat("Variables:", paste(model_vars, collapse = ", "), "\n")
  cat("Starting", n_boot, "bootstrap iterations...\n")
  
  # Run block bootstrap for this model - SEQUENTIAL VERSION
  list_boot <- list()
  for(i in 1:n_boot) {
    if(i %% 100 == 0) {  # Progress indicator every 100 iterations
      cat("Completed", i, "iterations...\n")
    }
    
    result <- tryCatch({
      boot_fit_model(df_ids = df_states,
                     df_full = dengue_temp,
                     model_vars = model_vars,
                     state_id_var = "state_residency",
                     seed = i)
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(result)) {
      list_boot[[length(list_boot) + 1]] <- result
    }
  }
  
  cat("Successful bootstrap iterations:", length(list_boot), "\n")
  
  if(length(list_boot) > 0) {
    # Combine results
    boot_results <- do.call("bind_rows", list_boot)
    boot_results$model_name <- model_name
    
    # Store results
    all_bootstrap_results[[model_name]] <- boot_results
    
    # Save individual model results
    write_csv(boot_results,
              paste0(model_name, "_coef_state_blockboot", n_boot, ".csv"))
    
    # Calculate bootstrap statistics
    boot_summary <- boot_results %>%
      select(-model_name) %>%
      summarise_all(list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        q025 = ~quantile(.x, 0.025, na.rm = TRUE),
        q975 = ~quantile(.x, 0.975, na.rm = TRUE)
      )) %>%
      mutate(model_name = model_name)
    
    # Save summary
    write_csv(boot_summary,
              paste0(model_name, "_coef_state_blockboot_summary", n_boot, ".csv"))
    
    cat("Bootstrap completed for", model_name, "\n")
  } else {
    cat("All bootstrap iterations failed for", model_name, "\n")
  }
  
  # Fit original model for comparison
  cat("Fitting original model for comparison...\n")
  
  # Separate fixed effects from regular predictors
  fixed_effects <- c()
  predictors <- c()
  
  for(var in model_vars) {
    if(var %in% c("city_residency", "month", "year", "month_childs", "year_childs")) {
      fixed_effects <- c(fixed_effects, var)
    } else {
      predictors <- c(predictors, var)
    }
  }
  
  # Create formula
  if(length(predictors) > 0 && length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "), "|", paste(fixed_effects, collapse = " + "))
  } else if(length(predictors) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "))
  } else if(length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~ 1 |", paste(fixed_effects, collapse = " + "))
  } else {
    formula_str <- "n_cases ~ 1"
  }
  
  model_formula <- as.formula(formula_str)
  
  original_model <- fixest::fepois(model_formula,
                                   offset = ~log_pop_offset,
                                   data = dengue_temp,
                                   nthreads = 1)  # Single thread for consistency
  
  original_coefs <- data.frame(t(coef(original_model)))
  original_coefs$model_name <- model_name
  all_original_results[[model_name]] <- original_coefs
}

cat("\n=== BOOTSTRAP SUMMARY ===\n")
cat("Models processed:", length(all_bootstrap_results), "\n")
cat("Models with successful bootstrap:", sum(sapply(all_bootstrap_results, function(x) nrow(x) > 0)), "\n")

# Combine all bootstrap results
if(length(all_bootstrap_results) > 0) {
  combined_bootstrap <- do.call("bind_rows", all_bootstrap_results)
  write_csv(combined_bootstrap, paste0("all_models_bootstrap_results", n_boot, ".csv"))
  cat("Combined bootstrap results saved to: all_models_bootstrap_results", n_boot, ".csv\n")
}

# Combine all original results
if(length(all_original_results) > 0) {
  combined_original <- do.call("bind_rows", all_original_results)
  write_csv(combined_original, "all_models_original_coefficients.csv")
  cat("Combined original results saved to: all_models_original_coefficients.csv\n")
}

cat("\nBootstrap analysis completed!\n")