# Load necessary packages
library(tidyverse)
library(magrittr)
library(splines)
library(glmmTMB)
library(dplyr)
library(MuMIn)  # For R-squared calculations

# Read data
dengue_temp <- read_csv("data/model_input_brazil_lag1_immunity_city.csv")

# Create population offset
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# Create B-spline basis functions for temperature lags
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

# Add B-spline columns to dataframe
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
    temp_bs_lag54 = temp_bs_fitted_lag5[,4]
  )

# make sure grouping vars are factors
dengue_temp <- dengue_temp %>%
  arrange(city_residency, year, month) %>%
  mutate(
    region = as.factor(region),
    year  = as.factor(year),
    month = as.factor(month(date_first_symptoms)),
    city_residency = as.factor(city_residency)
  )

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

# Fitting function for pre-sampled bootstrapped samples (now with AIC)
boot_fit_model_rsq_presampled <- function(df_full, # full dataset
                                          bootstrap_iter, # bootstrap iteration
                                          bootstrap_data, # pre-sampled states
                                          model_vars, # variables for this model
                                          state_id_var = "state_residency"){

  cat("Bootstrap iteration", bootstrap_iter, "\n")
  # Get bootstrapped sample using pre-sampled states
  boot_data <- boot_strat_presampled(df_full, bootstrap_iter, bootstrap_data, state_id_var)
  
  model <- glmmTMB(
    n_cases ~
      temp_bs_lag11 + temp_bs_lag12 + temp_bs_lag13 + temp_bs_lag14 +
      temp_bs_lag21 + temp_bs_lag22 + temp_bs_lag23 + temp_bs_lag24 +
      temp_bs_lag31 + temp_bs_lag32 + temp_bs_lag33 + temp_bs_lag34 +
      temp_bs_lag41 + temp_bs_lag42 + temp_bs_lag43 + temp_bs_lag44 +
      temp_bs_lag51 + temp_bs_lag52 + temp_bs_lag53 + temp_bs_lag54 +
      total_precipitation_lag1 + total_precipitation_lag2 + total_precipitation_lag3 +
      total_precipitation_lag4 + total_precipitation_lag5 +
      (1 | city_residency) + (1 | year) + ar1(month + 0 | region),
    offset = log_pop_offset,
    family = poisson(link = "log"),
    data = boot_data
  )

  print(summary(model))
  
  null_model <- glmmTMB(n_cases ~ 1, family = poisson(link = "log"), data = boot_data)

  mcfadden_r2 <- 1 - as.numeric(logLik(model) / logLik(null_model))  
  
  fitted_vals   <- fitted(model)
  observed_vals <- model.response(model.frame(model))
  correlation_r2 <- cor(observed_vals, fitted_vals)^2
  overall_municipality_rmse <- sqrt(mean((observed_vals - fitted_vals)^2, na.rm = TRUE))
  
  model_aic <- AIC(model)
  
  base_stats <- data.frame(
    mcfadden_r2    = mcfadden_r2,
    correlation_r2 = correlation_r2,
    aic            = model_aic,
    rmse           = overall_municipality_rmse,
    n_obs          = nobs(model),
    bootstrap_iteration = bootstrap_iter
  )
  
  # Extract fixed effects coefficients
  coefficients <- fixef(model)$cond  # $cond gets the conditional model coefficients
  
  # Convert coefficients to a data frame row (transpose the named vector)
  coef_df <- as.data.frame(t(coefficients))

  stats <- cbind(base_stats, coef_df)

  print(stats)
  
  return(stats)

}

# Get number of bootstrap iterations from the loaded data
n_boot <- 2 #max(bootstrap_states$bootstrap_iteration)
cat("Using", n_boot, "pre-sampled bootstrap iterations\n")

# Storage for all results
all_bootstrap_results <- list()
model_name <- "glmm"

  
# Run block bootstrap for this model using pre-sampled states - SEQUENTIAL VERSION
list_boot <- list()
for(i in 1:n_boot) {
  result <- boot_fit_model_rsq_presampled(df_full = dengue_temp,
                                  bootstrap_iter = i,
                                  bootstrap_data = bootstrap_states,
                                  model_vars = model_vars,
                                  state_id_var = "state_residency")

  
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
            paste0("data/glmm_rsquared_aic_state_blockboot_presampled", n_boot, ".csv"))  
} else {
  cat("All bootstrap iterations failed for", model_name, "\n")
}
