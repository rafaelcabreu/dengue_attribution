# run_validation_array.R
# Called by SLURM array job with arguments: model_index start_iter end_iter results_dir
# Model index is 1-based (1 = first model, 19 = last model)
# Example: Rscript run_validation_array.R 1 1 25 results/dengue_validation

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_validation_array.R <model_index> <start_iter> <end_iter> [results_dir]")
}

model_index <- as.integer(args[1])
start_iter  <- as.integer(args[2])
end_iter    <- as.integer(args[3])
results_dir <- if (length(args) >= 4) args[4] else "results/dengue_validation"

# Load required libraries
library(tidyverse)
library(magrittr)
library(fixest)
library(splines)

# Load and prepare initial data
dengue_temp <- read_csv("../data/model_input_brazil_immunity_city_with_priorinf.csv")
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)
dengue_temp$date_first_symptoms <- as.Date(dengue_temp$date_first_symptoms)

# ── Model definitions ──────────────────────────────────────────────────────────
models_with_temp <- list(
  'Year + Month' = c(
    'month',
    'year',
    'city_residency'
  ),
  'Year|region + Month|region' = c(
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5)' = c(
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
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region' = c(
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
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Month|region' = c(
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
    'month_region',
    'city_residency'
  ),
  'Year + Month + P(Climate)' = c(
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
    'city_residency',
    'month',
    'year'
  ),
  'Climate(lag 1) + Year|region + Month|region' = c(
    'temp_bs_lag11',
    'temp_bs_lag12',
    'temp_bs_lag13',
    'temp_bs_lag14',  
    'total_precipitation_lag1',
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-3) + Year|region + Month|region' = c(
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
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region + Month|region' = c(
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
    'month_region',
    'year_region',
    'city_residency'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + Socio' = c(
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
    'month_region',
    'year_region',
    'city_residency',
    'urban_area_ha',
    'gdp_per_capita',
    'births'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases' = c(
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
    'month_region',
    'year_region',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + SeroRepla' = c(
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
    'month_region',
    'year_region',
    'city_residency',
    'serotype_replacement'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + Immunity' = c(
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
    'month_region',
    'year_region',
    'city_residency',
    'Pr_0priorinf',
    'Pr_1priorinf'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
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
    'year_region',
    'month_region'
  ),
  'Climate(lag 1-5) + Socio' = c(
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
    'urban_area_ha',
    'gdp_per_capita',
    'births'
  ),
  'Climate(lag 1-5) + PriorCases' = c(
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
    'immunity_lag3'
  ),
  'Climate(lag 1-5) + SeroRepla' = c(
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
    'serotype_replacement'
  ),
  'Climate(lag 1-5) + Immunity' = c(
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
    'Pr_0priorinf',
    'Pr_1priorinf'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity' = c(
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
    'Pr_1priorinf'
  ),
  'Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
    'month_region',
    'city_residency',
    'immunity_lag1',
    'immunity_lag2',
    'immunity_lag3',
    'serotype_replacement',
    'urban_area_ha',
    'gdp_per_capita',
    'births',
    'Pr_0priorinf',
    'Pr_1priorinf'
  ),
  'Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity' = c(
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
    'year_region',
    'month_region'
  ),
  'Climate(lag 1-5) + PriorCases + SeroRepla + Socio' = c(
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
    'births'
  ),
  'Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio' = c(
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
    'year_region',
    'month_region'
  )
)

# ── Resolve model by index ─────────────────────────────────────────────────────
n_models <- length(models_with_temp)
if (model_index < 1 || model_index > n_models) {
  stop(paste("model_index must be between 1 and", n_models, "(got", model_index, ")"))
}

model_name <- names(models_with_temp)[model_index]
model_vars <- models_with_temp[[model_index]]

cat("Model index:", model_index, "\n")
cat("Model name: ", model_name, "\n")
cat("Iterations: ", start_iter, "to", end_iter, "\n")
cat("Results dir:", results_dir, "\n\n")

# ── helpers ────────────────────────────────────────────────────────────────────

build_formula <- function(model_vars) {
  # All variables treated as fixed effects (absorbed) by fixest
  fixed_effects <- c('city_residency', 'year', 'month',
                     'year_region', 'month_region',
                     'year_childs', 'month_childs')
  fe_vars      <- intersect(model_vars, fixed_effects)
  regular_vars <- setdiff(model_vars, fixed_effects)

  rhs <- if (length(regular_vars) > 0) paste(regular_vars, collapse = " + ") else "1"

  if (length(fe_vars) > 0) {
    as.formula(paste("n_cases ~", rhs, "|", paste(fe_vars, collapse = " + ")))
  } else {
    as.formula(paste("n_cases ~", rhs))
  }
}

calculate_pseudo_r2 <- function(model, validation_data) {
  tryCatch({
    predictions     <- predict(model, newdata = validation_data)
    observed        <- validation_data$n_cases
    ll_full         <- sum(dpois(observed, predictions, log = TRUE), na.rm = TRUE)
    null_prediction <- mean(observed, na.rm = TRUE)
    ll_null         <- sum(dpois(observed, null_prediction, log = TRUE), na.rm = TRUE)
    1 - (ll_full / ll_null)
  }, error = function(e) NA)
}

run_iteration <- function(iteration, model_name, model_vars) {
  set.seed(123 + iteration)

  selected_cities      <- sample(unique(dengue_temp$city_residency),
                                 round(length(unique(dengue_temp$city_residency)) * 0.2))
  selected_cities_data <- dengue_temp[dengue_temp$city_residency %in% selected_cities, ]
  available_years      <- sort(unique(year(selected_cities_data$date_first_symptoms)))

  if (length(available_years) < 10) {
    return(data.frame(model = model_name, iteration = iteration,
                      mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
                      training_years = NA, validation_years = NA,
                      note = "Insufficient years of data"))
  }

  max_start_year   <- available_years[length(available_years) - 9]
  possible_start   <- available_years[available_years <= max_start_year]
  start_year       <- sample(possible_start, 1)
  training_years   <- start_year:(start_year + 4)
  validation_years <- (start_year + 5):(start_year + 9)

  training_selected <- selected_cities_data[
    year(selected_cities_data$date_first_symptoms) %in% training_years, ]
  other_cities_data <- dengue_temp[!dengue_temp$city_residency %in% selected_cities, ]
  training_data     <- rbind(other_cities_data, training_selected)

  validation_data <- selected_cities_data[
    year(selected_cities_data$date_first_symptoms) %in% validation_years, ]

  if (nrow(validation_data) == 0) {
    return(data.frame(model = model_name, iteration = iteration,
                      mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
                      training_years = NA, validation_years = NA,
                      note = "No validation data"))
  }

  # B-splines (only if needed)
  if (any(grepl("temp_bs_lag", model_vars))) {
    temp_bs_lag1 <- bs(training_data$mean_2m_air_temp_degree1_lag1, df = 4)
    temp_bs_lag2 <- bs(training_data$mean_2m_air_temp_degree1_lag2, df = 4)
    temp_bs_lag3 <- bs(training_data$mean_2m_air_temp_degree1_lag3, df = 4)
    temp_bs_lag4 <- bs(training_data$mean_2m_air_temp_degree1_lag4, df = 4)
    temp_bs_lag5 <- bs(training_data$mean_2m_air_temp_degree1_lag5, df = 4)

    add_splines <- function(df, bs1, bs2, bs3, bs4, bs5) {
      df %>% mutate(
        temp_bs_lag11 = bs1[,1], temp_bs_lag12 = bs1[,2],
        temp_bs_lag13 = bs1[,3], temp_bs_lag14 = bs1[,4],
        temp_bs_lag21 = bs2[,1], temp_bs_lag22 = bs2[,2],
        temp_bs_lag23 = bs2[,3], temp_bs_lag24 = bs2[,4],
        temp_bs_lag31 = bs3[,1], temp_bs_lag32 = bs3[,2],
        temp_bs_lag33 = bs3[,3], temp_bs_lag34 = bs3[,4],
        temp_bs_lag41 = bs4[,1], temp_bs_lag42 = bs4[,2],
        temp_bs_lag43 = bs4[,3], temp_bs_lag44 = bs4[,4],
        temp_bs_lag51 = bs5[,1], temp_bs_lag52 = bs5[,2],
        temp_bs_lag53 = bs5[,3], temp_bs_lag54 = bs5[,4]
      )
    }

    training_data <- add_splines(training_data,
                                 temp_bs_lag1, temp_bs_lag2, temp_bs_lag3,
                                 temp_bs_lag4, temp_bs_lag5)

    v1 <- predict(temp_bs_lag1, validation_data$mean_2m_air_temp_degree1_lag1)
    v2 <- predict(temp_bs_lag2, validation_data$mean_2m_air_temp_degree1_lag2)
    v3 <- predict(temp_bs_lag3, validation_data$mean_2m_air_temp_degree1_lag3)
    v4 <- predict(temp_bs_lag4, validation_data$mean_2m_air_temp_degree1_lag4)
    v5 <- predict(temp_bs_lag5, validation_data$mean_2m_air_temp_degree1_lag5)

    validation_data <- add_splines(validation_data, v1, v2, v3, v4, v5)
  }

  model_formula <- build_formula(model_vars)

  tryCatch({
    main <- fixest::fepois(model_formula,
                           offset = ~log_pop_offset,
                           data   = training_data,
                           combine.quick = FALSE)

    predictions <- predict(main, newdata = validation_data)
    observed    <- validation_data$n_cases

    data.frame(
      model            = model_name,
      iteration        = iteration,
      mae              = mean(abs(observed - predictions), na.rm = TRUE),
      rmse             = sqrt(mean((observed - predictions)^2, na.rm = TRUE)),
      correlation      = cor(observed, predictions, use = "complete.obs"),
      pseudo_r2        = calculate_pseudo_r2(main, validation_data),
      training_years   = paste(training_years, collapse = "-"),
      validation_years = paste(validation_years, collapse = "-"),
      note             = "Success"
    )
  }, error = function(e) {
    data.frame(model = model_name, iteration = iteration,
               mae = NA, rmse = NA, correlation = NA, pseudo_r2 = NA,
               training_years = NA, validation_years = NA,
               note = paste("Error:", e$message))
  })
}

# ── Run the requested iterations ──────────────────────────────────────────────

results <- data.frame()
for (i in start_iter:end_iter) {
  if (i %% 10 == 0) cat("  Iteration", i, "\n")
  results <- rbind(results, run_iteration(i, model_name, model_vars))
}

# Save chunk results — use zero-padded index in filename to avoid special character issues
safe_name <- paste0("model", sprintf("%02d", model_index))
out_file  <- file.path(results_dir,
                       paste0("validation_", safe_name,
                              "_iter", start_iter, "_", end_iter, ".csv"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, out_file, row.names = FALSE)
cat("Saved:", out_file, "\n")