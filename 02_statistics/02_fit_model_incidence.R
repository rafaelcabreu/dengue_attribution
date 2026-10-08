# 02_fit_model_incidence.R
#
# Single-fit (no bootstrap) point estimate of the
# "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity"
# model (i.e. C15-Yr-Mr-c(PC,SR,S,IM)), used by:
#   - FigureS12.ipynb (municipality-level observed vs. fitted incidence)
#   - Figure2.ipynb (fixed-effects regional timeseries central line)
#
# Moved here from impacts/dengue_revision/02_fit_model_incidence.R.
#
# Outputs:
#   data/dengue_with_fitted_values.csv  — individual (city-month) level, with
#                                          fitted incidence/cases added
#   data/dengue_regional_aggregation.csv — aggregated to region x date, with
#                                           observed vs. fitted incidence and
#                                           prediction error metrics

# load packages
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

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

# Fit your original model
main <- fixest::fepois(n_cases ~
  temp_bs_lag11 + temp_bs_lag12 + temp_bs_lag13 + temp_bs_lag14 +
  temp_bs_lag21 + temp_bs_lag22 + temp_bs_lag23 + temp_bs_lag24 +
  temp_bs_lag31 + temp_bs_lag32 + temp_bs_lag33 + temp_bs_lag34 +
  temp_bs_lag41 + temp_bs_lag42 + temp_bs_lag43 + temp_bs_lag44 +
  temp_bs_lag51 + temp_bs_lag52 + temp_bs_lag53 + temp_bs_lag54 +
  total_precipitation_lag1 + total_precipitation_lag2 + total_precipitation_lag3 +
  total_precipitation_lag4 + total_precipitation_lag5 + immunity_lag1 + immunity_lag2 +
  immunity_lag3 + serotype_replacement + urban_area_ha + gdp_per_capita + births +
  Pr_0priorinf + Pr_1priorinf + Pr_2priorinf |
  city_residency + year_region + month_region,
  offset = ~log_pop_offset, # population weight
  data = dengue_temp,
  combine.quick = FALSE,
)

print(summary(main))

# ---- McFadden's pseudo-R^2 with proper null (same RE/AR1/offset) ----
null_model <- fixest::fepois(n_cases ~
  1,
  #offset = ~log_pop_offset, # population weight
  data = dengue_temp,
  combine.quick = FALSE
)

print(summary(null_model))

mcfadden_r2 <- 1 - as.numeric(logLik(main) / logLik(null_model))

print(mcfadden_r2)

# Get predictions for all data (handles dropped observations)
pred_cases <- predict(main, newdata = dengue_temp)

# Create dataset with fitted values
data_with_fitted <- dengue_temp %>%
  mutate(
    pred_incidence = pred_cases / population * 100000, # Convert cases to incindence
    pred_cases = pred_cases,
    actual_cases = n_cases        # Convert actual incidence to cases
  )

# 1. Overall municipality-level RMSE (across all observations)
overall_municipality_rmse <- sqrt(mean((data_with_fitted$n_cases - data_with_fitted$pred_cases)^2, na.rm = TRUE))

print(overall_municipality_rmse)

print(r2(main, type = "pr2"))

# File 1: Save fitted values dataset
write_csv(data_with_fitted, "../data/dengue_with_fitted_values.csv")
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
write_csv(regional_aggregation, "../data/dengue_regional_aggregation.csv")
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
