# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(foreach)
library(doParallel)

# load data
dengue_temp <- read_csv("data/model_input_brazil_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# Load pre-sampled bootstrap states (first 10 iterations)
bootstrap_states <- read_csv2("data/bootstrap_state_samples.csv", col_names = TRUE, locale = locale(encoding = "UTF-8"), n_max = 1000)

# Clean up the bootstrap states data
bootstrap_states <- bootstrap_states %>%
  separate_rows(resampled_states, sep = ",") %>%
  mutate(resampled_states = str_trim(resampled_states)) %>%
  filter(!is.na(resampled_states), resampled_states != "") %>%
  group_by(bootstrap_iteration) %>%
  mutate(state_boot_id = row_number()) %>%
  ungroup()

# Check the loaded bootstrap data
cat("Loaded bootstrap samples for", max(bootstrap_states$bootstrap_iteration), "iterations\n")

# Define models with temperature terms
models_with_temp <- list(
  'base' = c(
    'month',
    'year',
    'city_residency'
  ),
  'standard' = c(
    'mean_2m_air_temp_degree1_lag1',
    'mean_2m_air_temp_degree1_lag2',
    'mean_2m_air_temp_degree1_lag3',
    'mean_2m_air_temp_degree2_lag1',
    'mean_2m_air_temp_degree2_lag2',
    'mean_2m_air_temp_degree2_lag3',
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
    'mean_2m_air_temp_degree1_lag1',
    'mean_2m_air_temp_degree2_lag1',
    'total_precipitation_lag1',
    'total_precipitation_lag2',
    'total_precipitation_lag3',
    'month',
    'year',
    'city_residency'
  ),
  'socioeconomic' = c(
    'mean_2m_air_temp_degree1_lag1',
    'mean_2m_air_temp_degree1_lag2',
    'mean_2m_air_temp_degree1_lag3',
    'mean_2m_air_temp_degree2_lag1',
    'mean_2m_air_temp_degree2_lag2',
    'mean_2m_air_temp_degree2_lag3',
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
  )
)

# Function to use pre-sampled states
boot_strat_presampled <- function(df_full, # full dataset
                                  bootstrap_iter, # which bootstrap iteration to use
                                  bootstrap_data, # pre-sampled bootstrap states
                                  state_id_var = "state_residency"){
  
  # Get the states for this bootstrap iteration
  states_for_iter <- bootstrap_data %>%
    filter(bootstrap_iteration == bootstrap_iter) %>%
    select(state_boot_id, state_name = resampled_states)
  
  # Check if the state column exists in df_full
  if(!state_id_var %in% colnames(df_full)) {
    return(data.frame())
  }
  
  # Create a named vector for the join
  join_by <- setNames("state_name", state_id_var)
  
  # Join back to full dataset using the pre-sampled states
  df_out <- df_full %>%
    inner_join(states_for_iter, 
               by = join_by,
               relationship = "many-to-many") %>%
    filter(!is.na(state_boot_id))
  
  return(df_out)
}

# Fitting function for pre-sampled bootstrapped samples
boot_fit_model_rsq_presampled <- function(df_full, # full dataset
                                          bootstrap_iter, # bootstrap iteration
                                          bootstrap_data, # pre-sampled states
                                          model_vars, # variables for this model
                                          state_id_var = "state_residency"){

  # Get bootstrapped sample using pre-sampled states
  boot_data <- boot_strat_presampled(df_full, bootstrap_iter, bootstrap_data, state_id_var)
  
  if(nrow(boot_data) == 0) {
    return(NULL)
  }

  # Check if required columns exist
  missing_vars <- setdiff(model_vars, names(boot_data))
  if(length(missing_vars) > 0) {
    return(NULL)
  }
  
  # Check if n_cases and population exist
  if(!"n_cases" %in% names(boot_data) || !"population" %in% names(boot_data)) {
    return(NULL)
  }

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

  tryCatch({
    # Fit the model with pre-sampled cluster structure
    model <- fixest::fepois(model_formula,
                            offset = ~log_pop_offset,
                            data = boot_data,
                            nthreads = 2)

    # Calculate R-squared using fixest's r2 function
    r_squared <- r2(model, type = "cor2")
    pseudo_r_squared <- r2(model, type = "pr2")  # Pseudo R-squared
    
    return(data.frame(
      r_squared = r_squared,
      pseudo_r_squared = pseudo_r_squared,
      n_obs = model$nobs,
      bootstrap_iteration = bootstrap_iter
    ))
  }, error = function(e) {
    return(NULL)
  })
}

# Set up parallel processing
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Get number of bootstrap iterations from the loaded data
n_boot <- max(bootstrap_states$bootstrap_iteration)
cat("Using", n_boot, "pre-sampled bootstrap iterations\n")

# Storage for all results
all_bootstrap_results <- list()

# Loop through each model
for(model_name in names(models_with_temp)) {
  cat("\n=== Running bootstrap for model:", model_name, "===\n")
  
  tryCatch({
    # Get variables for this model
    model_vars <- models_with_temp[[model_name]]
    
    cat("Variables:", paste(model_vars, collapse = ", "), "\n")
    cat("Starting", n_boot, "bootstrap iterations using pre-sampled states...\n")
    
    # Run block bootstrap for this model using pre-sampled states
    list_boot <- foreach(i = 1:n_boot,
                         .packages = c("tidyverse", "magrittr", "fixest"),
                         .export = c("bootstrap_states"),
                         .errorhandling = "remove") %dopar% {
      
      tryCatch({
        boot_fit_model_rsq_presampled(df_full = dengue_temp,
                                      bootstrap_iter = i,
                                      bootstrap_data = bootstrap_states,
                                      model_vars = model_vars,
                                      state_id_var = "state_residency")
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
                paste0(model_name, "_rsquared_state_blockboot_presampled", n_boot, ".csv"))
      
      # Calculate bootstrap statistics
      boot_summary <- boot_results %>%
        select(r_squared, pseudo_r_squared, n_obs) %>%
        summarise(
          r_squared_mean = mean(r_squared, na.rm = TRUE),
          r_squared_sd = sd(r_squared, na.rm = TRUE),
          r_squared_q025 = quantile(r_squared, 0.025, na.rm = TRUE),
          r_squared_q975 = quantile(r_squared, 0.975, na.rm = TRUE),
          pseudo_r_squared_mean = mean(pseudo_r_squared, na.rm = TRUE),
          pseudo_r_squared_sd = sd(pseudo_r_squared, na.rm = TRUE),
          pseudo_r_squared_q025 = quantile(pseudo_r_squared, 0.025, na.rm = TRUE),
          pseudo_r_squared_q975 = quantile(pseudo_r_squared, 0.975, na.rm = TRUE),
          mean_n_obs = mean(n_obs, na.rm = TRUE),
          n_valid_iterations = sum(!is.na(r_squared))
        ) %>%
        mutate(model_name = model_name)
      
      # Save summary
      write_csv(boot_summary,
                paste0(model_name, "_rsquared_state_blockboot_presampled_summary", n_boot, ".csv"))
      
      cat("Bootstrap completed for", model_name, "\n")
      cat("Mean R-squared:", round(boot_summary$r_squared_mean, 4), "\n")
      cat("Mean Pseudo R-squared:", round(boot_summary$pseudo_r_squared_mean, 4), "\n")
    } else {
      cat("All bootstrap iterations failed for", model_name, "\n")
    }
    
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
  write_csv(combined_bootstrap, paste0("data/all_models_bootstrap_rsquared_presampled", n_boot, ".csv"))
  cat("Combined bootstrap R-squared results saved to: data/all_models_bootstrap_rsquared_presampled", n_boot, ".csv\n")
}

cat("\nBootstrap R-squared analysis with pre-sampled states completed!\n")