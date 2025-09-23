# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

# Load data
dengue_temp <- read_csv("data/model_input_brazil_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)
  
# Create B-spline basis functions for each temperature lag
# Extract individual columns from each B-spline matrix
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
    temp_bs_lag54 = temp_bs_fitted_lag5[,4]
  )

# Define models
models <- list(
  'Year + Month' = c(
    'month',
    'year',
    'city_residency'
  ),
  'Year + Month + Climate' = c(
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
  'Year + Month + Climate + Immunity' = c(
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
  'Year + Month + Climate(lag1)' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',
    'total_precipitation_lag1',
    'year',
    'month',
    'city_residency'
  ),
  'Year + Month + Climate(lag1-5)' = c(
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
  'Childs et al. (2025)' = c(
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
  'Year + Month + Climate + Socioeconomic' = c(
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
    'urban_pop',
    'percent_urban'
  )
)

# Initialize results storage
model_results <- list()
model_statistics <- data.frame()

# Fit each model
for(model_name in names(models)) {
  cat("\n=== Fitting model:", model_name, "===\n")
  
  tryCatch({
    # Get variables for this model
    model_vars <- models[[model_name]]
    
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
    cat("Formula:", formula_str, "\n")
    
    # Fit model
    fitted_model <- fixest::fepois(
      model_formula,
      offset = ~log_pop_offset,
      data = dengue_temp,
      combine.quick = FALSE
    )
    
    # Store model
    model_results[[model_name]] <- fitted_model
    
    # Extract statistics
    model_stats <- data.frame(
      model_name = model_name,
      n_obs = nobs(fitted_model),
      log_likelihood = as.numeric(logLik(fitted_model)),
      aic = AIC(fitted_model),
      bic = BIC(fitted_model),
      pseudo_r2 = r2(fitted_model, type = "pr2"),
      adj_pseudo_r2 = r2(fitted_model, type = "apr2"),
      squared_cor = r2(fitted_model, type = "cor2"),
      n_variables = length(model_vars),
      stringsAsFactors = FALSE
    )
    
    # Add to results
    model_statistics <- rbind(model_statistics, model_stats)
    
    # Print summary
    cat("Model fitted successfully!\n")
    cat("Log-likelihood:", model_stats$log_likelihood, "\n")
    cat("AIC:", model_stats$aic, "\n")
    cat("BIC:", model_stats$bic, "\n")
    cat("Pseudo R²:", model_stats$pseudo_r2, "\n")
    cat("Squared Correlation:", model_stats$squared_cor, "\n")
    
  }, error = function(e) {
    cat("ERROR fitting model", model_name, ":", e$message, "\n")
    
    # Add failed model to statistics
    model_stats <- data.frame(
      model_name = model_name,
      n_obs = NA,
      log_likelihood = NA,
      aic = NA,
      bic = NA,
      pseudo_r2 = NA,
      adj_pseudo_r2 = NA,
      squared_cor = NA,
      n_variables = length(models[[model_name]]),
      stringsAsFactors = FALSE
    )
    model_statistics <<- rbind(model_statistics, model_stats)
  })
}

# Display final results
cat("\n=== FINAL MODEL COMPARISON ===\n")
print(model_statistics)

# Sort by AIC (lower is better)
model_statistics_sorted <- model_statistics %>%
  arrange(aic) %>%
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE),
    rank_aic = rank(aic, na.last = NA)
  )

cat("\n=== MODELS RANKED BY AIC ===\n")
print(model_statistics_sorted[c("model_name", "aic", "delta_aic", "rank_aic", "pseudo_r2", "n_variables")])

# Save results
write_csv(model_statistics_sorted, "model_comparison_statistics.csv")
cat("\nModel statistics saved to: model_comparison_statistics.csv\n")

# Print best model
best_model_name <- model_statistics_sorted$model_name[1]
cat("\nBest model (lowest AIC):", best_model_name, "\n")

# Print summary of best model
if(best_model_name %in% names(model_results)) {
  cat("\n=== SUMMARY OF BEST MODEL ===\n")
  print(summary(model_results[[best_model_name]]))
}