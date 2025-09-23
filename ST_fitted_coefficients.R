# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

dengue_temp <- read_csv("../model_input_brazil_immunity_city.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# Create B-spline basis functions for temperature
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

# Add B-spline columns to the dataframe
dengue_temp <- dengue_temp %>%
    mutate(
      # Temperature B-splines
      temp_bs_lag11 = temp_bs_fitted_lag1[,1],
      temp_bs_lag12 = temp_bs_fitted_lag1[,2],
      temp_bs_lag13 = temp_bs_fitted_lag1[,3],
      temp_bs_lag14 = temp_bs_fitted_lag1[,4],
      temp_bs_lag21 = temp_bs_fitted_lag2[,1],
      temp_bs_lag22 = temp_bs_fitted_lag2[,2],
      temp_bs_lag23 = temp_bs_fitted_lag2[,3],
      temp_bs_lag24 = temp_bs_fitted_lag2[,4],
      temp_bs_lag31 = temp_bs_fitted_lag3[,1],
      temp_bs_lag32 = temp_bs_fitted_lag3[,2],
      temp_bs_lag33 = temp_bs_fitted_lag3[,3],
      temp_bs_lag34 = temp_bs_fitted_lag3[,4],
      temp_bs_lag41 = temp_bs_fitted_lag4[,1],
      temp_bs_lag42 = temp_bs_fitted_lag4[,2],
      temp_bs_lag43 = temp_bs_fitted_lag4[,3],
      temp_bs_lag44 = temp_bs_fitted_lag4[,4],
      temp_bs_lag51 = temp_bs_fitted_lag5[,1],
      temp_bs_lag52 = temp_bs_fitted_lag5[,2],
      temp_bs_lag53 = temp_bs_fitted_lag5[,3],
      temp_bs_lag54 = temp_bs_fitted_lag5[,4]
    )

# Define models with temperature terms
models_with_temp <- list(
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

# Function to extract coefficients from a model
extract_coefficients <- function(model, model_name) {
  tryCatch({
    # Use coeftable() function which is more reliable for fixest objects
    coef_table <- coeftable(model)
    
    # Check if coefficients exist
    if(is.null(coef_table) || nrow(coef_table) == 0) {
      cat("Warning: No coefficients found for", model_name, "\n")
      return(data.frame(
        model = model_name,
        variable = "No coefficients",
        coefficient = NA,
        std_error = NA,
        z_value = NA,
        p_value = NA,
        significant = "",
        ci_lower = NA,
        ci_upper = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    # Debug: print column names to see structure
    cat("Coefficient table columns for", model_name, ":", colnames(coef_table), "\n")
    
    # Create data frame with coefficient information
    coef_df <- data.frame(
      model = model_name,
      variable = rownames(coef_table),
      stringsAsFactors = FALSE
    )
    
    # Extract columns by name to handle different formats
    col_names <- colnames(coef_table)
    
    # Coefficient estimate
    if("Estimate" %in% col_names) {
      coef_df$coefficient <- coef_table[, "Estimate"]
    } else if(ncol(coef_table) >= 1) {
      coef_df$coefficient <- coef_table[, 1]
    } else {
      coef_df$coefficient <- NA
    }
    
    # Standard error
    if("Std. Error" %in% col_names) {
      coef_df$std_error <- coef_table[, "Std. Error"]
    } else if(ncol(coef_table) >= 2) {
      coef_df$std_error <- coef_table[, 2]
    } else {
      coef_df$std_error <- NA
    }
    
    # Z-value (or t-value)
    if("z value" %in% col_names) {
      coef_df$z_value <- coef_table[, "z value"]
    } else if("t value" %in% col_names) {
      coef_df$z_value <- coef_table[, "t value"]
    } else if(ncol(coef_table) >= 3) {
      coef_df$z_value <- coef_table[, 3]
    } else {
      coef_df$z_value <- NA
    }
    
    # P-value
    if("Pr(>|z|)" %in% col_names) {
      coef_df$p_value <- coef_table[, "Pr(>|z|)"]
    } else if("Pr(>|t|)" %in% col_names) {
      coef_df$p_value <- coef_table[, "Pr(>|t|)"]
    } else if(ncol(coef_table) >= 4) {
      coef_df$p_value <- coef_table[, 4]
    } else {
      coef_df$p_value <- NA
    }
    
    # Add significance stars and confidence intervals
    coef_df <- coef_df %>%
      mutate(
        significant = case_when(
          is.na(p_value) ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          p_value < 0.1 ~ ".",
          TRUE ~ ""
        ),
        # Add confidence intervals
        ci_lower = ifelse(is.na(std_error), NA, coefficient - 1.96 * std_error),
        ci_upper = ifelse(is.na(std_error), NA, coefficient + 1.96 * std_error)
      )
    
    return(coef_df)
    
  }, error = function(e) {
    cat("Error extracting coefficients for", model_name, ":", e$message, "\n")
    
    # Try alternative approach using summary() function
    tryCatch({
      model_summary <- summary(model)
      if(!is.null(model_summary$coefficients)) {
        coef_table <- model_summary$coefficients
        
        return(data.frame(
          model = model_name,
          variable = rownames(coef_table),
          coefficient = coef_table[, 1],
          std_error = if(ncol(coef_table) >= 2) coef_table[, 2] else NA,
          z_value = if(ncol(coef_table) >= 3) coef_table[, 3] else NA,
          p_value = if(ncol(coef_table) >= 4) coef_table[, 4] else NA,
          significant = "",
          ci_lower = NA,
          ci_upper = NA,
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e2) {
      cat("Summary extraction also failed for", model_name, "\n")
    })
    
    # Last resort: try basic coef() function
    tryCatch({
      coefs <- coef(model)
      if(length(coefs) > 0) {
        return(data.frame(
          model = model_name,
          variable = names(coefs),
          coefficient = as.numeric(coefs),
          std_error = NA,
          z_value = NA,
          p_value = NA,
          significant = "",
          ci_lower = NA,
          ci_upper = NA,
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e3) {
      cat("All extraction methods failed for", model_name, "\n")
    })
    
    return(NULL)
  })
}

# Initialize results storage
all_coefficients <- data.frame()

# Fit all models and extract coefficients
cat("=== EXTRACTING COEFFICIENTS FOR ALL MODELS ===\n")
for(model_name in names(models_with_temp)) {
  cat("Fitting model:", model_name, "\n")
  
  model_vars <- models_with_temp[[model_name]]
  
  # Build and fit the model
  model_formula <- build_formula(model_vars)
  
  tryCatch({
    fitted_model <- fixest::fepois(
      model_formula,
      offset = ~log_pop_offset,
      data = dengue_temp,
      combine.quick = FALSE
    )

    print(summary(fitted_model))
    
    # Extract coefficients
    coef_df <- extract_coefficients(fitted_model, model_name)
    if(!is.null(coef_df)) {
      all_coefficients <- rbind(all_coefficients, coef_df)
      cat("  - Successfully extracted", nrow(coef_df), "coefficients\n")
    }
    
  }, error = function(e) {
    cat("  - Error fitting model:", e$message, "\n")
  })
  
  cat("\n")
}

# Save all coefficients to CSV
write_csv(all_coefficients, "model_coefficients_all.csv")
cat("All coefficients saved to: model_coefficients_all.csv\n")

# Display summary of coefficients
cat("\n=== COEFFICIENT SUMMARY ===\n")
coef_summary <- all_coefficients %>%
  group_by(model) %>%
  summarise(
    n_coefficients = n(),
    n_significant_05 = sum(p_value < 0.05, na.rm = TRUE),
    n_significant_01 = sum(p_value < 0.01, na.rm = TRUE),
    n_significant_001 = sum(p_value < 0.001, na.rm = TRUE),
    .groups = 'drop'
  )

print(coef_summary)

# Show most significant coefficients across all models
cat("\n=== MOST SIGNIFICANT COEFFICIENTS (p < 0.001) ===\n")
most_significant <- all_coefficients %>%
  filter(p_value < 0.001) %>%
  arrange(p_value) %>%
  select(model, variable, coefficient, p_value, significant) %>%
  head(20)

if(nrow(most_significant) > 0) {
  print(most_significant)
} else {
  cat("No coefficients with p < 0.001 found\n")
}

cat("\n=== FILE CREATED ===\n")
cat("model_coefficients_all.csv - Contains all model coefficients with statistics\n")