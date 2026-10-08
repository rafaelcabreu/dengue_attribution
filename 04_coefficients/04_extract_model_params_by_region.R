# 04_extract_model_params_by_region.R
#
# Fits the "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla
# + Socio + Immunity" model (C15-Yr-Mr-c(PC,SR,S,IM)) once nationally and once per
# region, and extracts coefficients/vcov/temperature-grid data for Python to
# compute temperature response curves and marginal effects — used by
# FigureS13.ipynb.
#
# Run from the repo root (impacts/dengue_attribution/); it writes its output
# CSVs (model_coefficients_<region>.csv, model_vcov_<region>.csv, etc.) directly
# into the CURRENT directory (no "data/" prefix), matching how FigureS13.ipynb
# reads them.
#
# Moved here from impacts/dengue_revision/extract_model_params_by_region.R.

library(tidyverse)
library(fixest)
library(splines)

cat("=== Extracting Model Parameters by Region for Python ===\n\n")

# Load data ----
cat("Loading data...\n")
dengue_temp <- read_csv("../data/model_input_brazil_immunity_city_with_priorinf.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population/100000)

# ── Municipality-level posterior immunity (posterior mean) ──────────────────────
# Replace the old macro-region-level Pr_0priorinf / Pr_1priorinf with
# municipality-level (city_residency x year) posterior-mean values — see
# 01_model_fit/build_immunity_city_draws.R. This is a point-estimate,
# no-bootstrap fit, so the posterior mean is used, same as fit_model.R.
city_health_region   <- read_csv("../data/code_health_region_to_city_residency.csv")
immunity_mean_region <- read_csv("../data/immunity_posterior_mean_region.csv")

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency") %>%
  left_join(immunity_mean_region, by = c("code_health_region", "year"))

cat("Data loaded:", nrow(dengue_temp), "observations\n")

# Get unique regions
regions <- unique(dengue_temp$region)
cat("Regions found:", paste(regions, collapse = ", "), "\n\n")

# Define region order
region_order <- c('North', 'Northeast', 'Middle-West', 'Southeast', 'South')
regions_to_process <- region_order[region_order %in% regions]

# Add national model to the processing list
regions_to_process <- c('National', regions_to_process)

# Define model variables
model_vars <- c(
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
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha',
    'gdp_per_capita',
    'births',
    'Pr_0priorinf',
    'Pr_1priorinf',
    'Pr_2priorinf',
    'year_region',
    'month_region'
)

# Separate fixed effects from predictors
fixed_effects <- c()
predictors <- c()
for(var in model_vars) {
  if(var %in% c("city_residency", "month", "year", "year_region", "month_region")) {
    fixed_effects <- c(fixed_effects, var)
  } else {
    predictors <- c(predictors, var)
  }
}

# Create formula once (same for all regions)
formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "),
                    "|", paste(fixed_effects, collapse = " + "))
model_formula <- as.formula(formula_str)
cat("Model formula:", formula_str, "\n\n")

# Loop through each region ----
for(region_name in regions_to_process) {
  cat("\n")
  cat("========================================\n")
  cat("PROCESSING REGION:", region_name, "\n")
  cat("========================================\n")

  # Filter data for this region (or use all data for National)
  if(region_name == "National") {
    region_data <- dengue_temp
  } else {
    region_data <- dengue_temp %>% filter(region == region_name)
  }
  cat("N observations:", nrow(region_data), "\n")

  # Create B-spline features for this region ----
  cat("Creating B-spline features...\n")
  temp_bs_fitted_lag1 <- bs(region_data$mean_2m_air_temp_degree1_lag1, df = 4)
  temp_bs_fitted_lag2 <- bs(region_data$mean_2m_air_temp_degree1_lag2, df = 4)
  temp_bs_fitted_lag3 <- bs(region_data$mean_2m_air_temp_degree1_lag3, df = 4)
  temp_bs_fitted_lag4 <- bs(region_data$mean_2m_air_temp_degree1_lag4, df = 4)
  temp_bs_fitted_lag5 <- bs(region_data$mean_2m_air_temp_degree1_lag5, df = 4)

  region_data <- region_data %>%
    mutate(
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

  # Fit model ----
  cat("Fitting model...\n")
  fitted_model <- fixest::fepois(
    model_formula,
    offset = ~log_pop_offset,
    data = region_data
  )

  cat("Model fitted successfully\n")

  # Extract parameters ----
  cat("Extracting model parameters...\n")

  # Coefficients
  coef_vec <- coef(fitted_model)
  coef_df <- data.frame(
    region = region_name,
    coefficient = names(coef_vec),
    value = as.numeric(coef_vec)
  )

  # Variance-covariance matrix
  vcov_mat <- vcov(fitted_model, vcov = ~city_residency)
  vcov_df <- as.data.frame(as.matrix(vcov_mat))
  vcov_df$region <- region_name
  vcov_df$row_name <- rownames(vcov_mat)

  # Temperature data for histogram
  temp_data <- region_data %>%
    select(mean_2m_air_temp_degree1_lag1, population) %>%
    filter(!is.na(mean_2m_air_temp_degree1_lag1)) %>%
    rename(temperature = mean_2m_air_temp_degree1_lag1) %>%
    mutate(region = region_name)

  # Calculate temperature range and mean
  temp_quantiles <- quantile(temp_data$temperature, c(0.01, 0.99), na.rm = TRUE)
  temp_mean <- weighted.mean(temp_data$temperature, temp_data$population, na.rm = TRUE)
  temp_range <- range(temp_data$temperature, na.rm = TRUE)

  # Create metadata
  metadata <- data.frame(
    region = region_name,
    parameter = c("temp_min", "temp_max", "temp_mean", "temp_q01", "temp_q99",
                  "n_obs", "aic", "bic", "pseudo_r2", "df"),
    value = c(temp_range[1], temp_range[2], temp_mean,
             temp_quantiles[1], temp_quantiles[2],
             nobs(fitted_model), AIC(fitted_model), BIC(fitted_model),
             r2(fitted_model, type = "pr2"), 4)
  )

  # Save B-spline information - dense grid for interpolation
  temp_grid <- seq(temp_range[1], temp_range[2], length.out = 1000)
  bs_grid_lag1 <- bs(temp_grid, df = 4,
                     Boundary.knots = attr(temp_bs_fitted_lag1, "Boundary.knots"))

  bs_basis <- data.frame(
    region = region_name,
    temperature = temp_grid,
    basis1 = bs_grid_lag1[,1],
    basis2 = bs_grid_lag1[,2],
    basis3 = bs_grid_lag1[,3],
    basis4 = bs_grid_lag1[,4]
  )

  # Save region-specific parameters ----
  region_safe <- gsub("-", "_", region_name)
  region_safe <- gsub(" ", "_", region_safe)

  cat("Saving parameters to CSV files...\n")
  write_csv(coef_df, paste0("model_coefficients_", region_safe, ".csv"))
  write_csv(vcov_df, paste0("model_vcov_", region_safe, ".csv"))
  write_csv(temp_data, paste0("temperature_data_", region_safe, ".csv"))
  write_csv(metadata, paste0("model_metadata_", region_safe, ".csv"))
  write_csv(bs_basis, paste0("bspline_basis_", region_safe, ".csv"))

  # Print summary
  cat("\n--- Model Summary for", region_name, "---\n")
  cat("Coefficients:", nrow(coef_df), "\n")
  cat("Temperature coefficients:", sum(grepl("temp_bs", coef_df$coefficient)), "\n")
  cat("Temperature range:", round(temp_range[1], 2), "to", round(temp_range[2], 2), "°C\n")
  cat("Mean temperature:", round(temp_mean, 2), "°C\n")
  cat("AIC:", round(AIC(fitted_model), 2), "\n")
  cat("Pseudo R²:", round(r2(fitted_model, type = "pr2"), 4), "\n")
  cat("Files saved:\n")
  cat("  - model_coefficients_", region_safe, ".csv\n", sep = "")
  cat("  - model_vcov_", region_safe, ".csv\n", sep = "")
  cat("  - temperature_data_", region_safe, ".csv\n", sep = "")
  cat("  - model_metadata_", region_safe, ".csv\n", sep = "")
  cat("  - bspline_basis_", region_safe, ".csv\n", sep = "")
}

cat("\n\n")
cat("========================================\n")
cat("PARAMETER EXTRACTION COMPLETE!\n")
cat("========================================\n")
cat("Processed", length(regions_to_process), "regions:\n")
cat(paste(regions_to_process, collapse = ", "), "\n")
cat("\nNow run the Python notebook to create plots.\n")
