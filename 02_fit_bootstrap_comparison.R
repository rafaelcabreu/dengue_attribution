# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(foreach)
library(doParallel)

# load data
dengue_temp <- read_csv("model_input_brazil_immunity_city.csv")

# Define models with temperature terms
models_with_temp <- list(
  'climate' = c(
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
    'month',
    'city_residency'
  ),
  'standard' = c(
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
    'month',
    'year',
    'city_residency'
  ),
  'with_immunity' = c(
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
    'month',
    'year',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3'
  ),
  'lag_45' = c(
    'mean_2m_air_temp_degree1_lag1',
    'mean_2m_air_temp_degree1_lag2',
    'mean_2m_air_temp_degree1_lag3',
    'mean_2m_air_temp_degree1_lag4',
    'mean_2m_air_temp_degree1_lag5',
    'mean_2m_air_temp_degree2_lag1',
    'mean_2m_air_temp_degree2_lag2',
    'mean_2m_air_temp_degree2_lag3',
    'mean_2m_air_temp_degree2_lag4',
    'mean_2m_air_temp_degree2_lag5',
    'mean_2m_air_temp_degree3_lag1',
    'mean_2m_air_temp_degree3_lag2',
    'mean_2m_air_temp_degree3_lag3',
    'mean_2m_air_temp_degree3_lag4',
    'mean_2m_air_temp_degree3_lag5',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'total_precipitation_lag4',
    'total_precipitation_lag5', 
    'month',
    'year',
    'city_residency'
  )
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
    'mean_2m_air_temp_degree1_lag1',
    'mean_2m_air_temp_degree2_lag1',
    'mean_2m_air_temp_degree3_lag1',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
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
    formula_str <- paste("dengue_inc ~", paste(predictors, collapse = " + "), "|", paste(fixed_effects, collapse = " + "))
  } else if(length(predictors) > 0) {
    formula_str <- paste("dengue_inc ~", paste(predictors, collapse = " + "))
  } else if(length(fixed_effects) > 0) {
    formula_str <- paste("dengue_inc ~ 1 |", paste(fixed_effects, collapse = " + "))
  } else {
    formula_str <- "dengue_inc ~ 1"
  }
  
  model_formula <- as.formula(formula_str)

  # Fit the model with new cluster structure
  model <- fixest::fepois(model_formula,
                          weights = ~population, # population weight
                          data = boot_data,
                          nthreads = 2)

  return(coef(model))
}

# Set up parallel processing
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Create dataset of unique states
df_states <- dengue_temp %>%
  select(state_residency) %>%
  unique() %>%
  filter(!is.na(state_residency))  # Remove any NA values

print(paste("Number of unique states:", nrow(df_states)))
print("Unique states:")
print(df_states$state_residency)

# Set number of bootstrap iterations
n_boot <- 1000 

# Storage for all results
all_bootstrap_results <- list()
all_original_results <- list()

# Loop through each model
for(model_name in names(models_with_temp)) {
  cat("\n=== Running bootstrap for model:", model_name, "===\n")
  
  tryCatch({
    # Get variables for this model
    model_vars <- models_with_temp[[model_name]]
    
    cat("Variables:", paste(model_vars, collapse = ", "), "\n")
    cat("Starting", n_boot, "bootstrap iterations...\n")
    
    # Run block bootstrap for this model
    list_boot <- foreach(i = 1:n_boot,
                         .packages = c("tidyverse", "magrittr", "fixest"),
                         .errorhandling = "remove") %dopar% {
      
      tryCatch({
        boot_fit_model(df_ids = df_states,
                       df_full = dengue_temp,
                       model_vars = model_vars,
                       state_id_var = "state_residency",
                       seed = i)
      }, error = function(e) {
        return(NULL)
      })
    }
    
    # Remove NULL results (failed iterations)
    list_boot <- list_boot[!sapply(list_boot, is.null)]
    
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
      formula_str <- paste("dengue_inc ~", paste(predictors, collapse = " + "), "|", paste(fixed_effects, collapse = " + "))
    } else if(length(predictors) > 0) {
      formula_str <- paste("dengue_inc ~", paste(predictors, collapse = " + "))
    } else if(length(fixed_effects) > 0) {
      formula_str <- paste("dengue_inc ~ 1 |", paste(fixed_effects, collapse = " + "))
    } else {
      formula_str <- "dengue_inc ~ 1"
    }
    
    model_formula <- as.formula(formula_str)
    
    original_model <- fixest::fepois(model_formula,
                                     weights = ~population,
                                     data = dengue_temp)
    
    original_coefs <- data.frame(t(coef(original_model)))
    original_coefs$model_name <- model_name
    all_original_results[[model_name]] <- original_coefs
    
  }, error = function(e) {
    cat("ERROR with model", model_name, ":", e$message, "\n")
  })
}

# Stop parallel processing
stopCluster(cl)

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
