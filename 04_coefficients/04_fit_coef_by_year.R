# 04_fit_coef_by_year.R
#
# Fit the flagship immunity model
# ("Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla +
# Socio + Immunity") once per calendar year, using only that year's data (no
# bootstrap — point estimate only). Loops through all years and saves
# coefficients to a single combined file.
#
# Used by FigureS14.ipynb (the "Year only" per-year temperature-response line,
# read from data/coefficients/coefs_by_year/all_coefs_by_year.csv).
#
# Moved here from impacts/dengue_revision/jasmin/S4_coefficients/coef_fit_by_year.R.

# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)

# ---------------------------------------------------------------
# 2. Paths and settings
# ---------------------------------------------------------------
data_path  <- "../data/model_input_brazil_immunity_city_with_priorinf.csv"
output_dir <- "../data/coefficients/coefs_by_year"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 3. Model definition
# ---------------------------------------------------------------
model_name <- "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity"
model_vars <- c(
  'temp_bs_lag11', 'temp_bs_lag12', 'temp_bs_lag13', 'temp_bs_lag14',
  'temp_bs_lag21', 'temp_bs_lag22', 'temp_bs_lag23', 'temp_bs_lag24',
  'temp_bs_lag31', 'temp_bs_lag32', 'temp_bs_lag33', 'temp_bs_lag34',
  'temp_bs_lag41', 'temp_bs_lag42', 'temp_bs_lag43', 'temp_bs_lag44',
  'temp_bs_lag51', 'temp_bs_lag52', 'temp_bs_lag53', 'temp_bs_lag54',
  'total_precipitation_lag1', 'total_precipitation_lag2', 'total_precipitation_lag3',
  'total_precipitation_lag4', 'total_precipitation_lag5',
  'city_residency',
  'immunity_lag1', 'immunity_lag2', 'immunity_lag3',
  'serotype_replacement',
  'urban_area_ha', 'gdp_per_capita', 'births',
  'Pr_0priorinf', 'Pr_1priorinf', 'Pr_2priorinf',
  'year_region', 'month_region'
)

# ---------------------------------------------------------------
# 4. Load and prepare full data
# ---------------------------------------------------------------
dengue_temp <- read_csv(data_path, show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

# ── Municipality-level posterior immunity (posterior mean) ──────────────────────
# Replace the old macro-region-level Pr_0priorinf / Pr_1priorinf with
# municipality-level (city_residency x year) posterior-mean values — see
# 01_model_fit/build_immunity_city_draws.R. This is a point-estimate,
# no-bootstrap fit, so the posterior mean is used, same as fit_model.R.
city_health_region   <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)
immunity_mean_region <- read_csv("../data/immunity_posterior_mean_region.csv", show_col_types = FALSE)

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf"))) %>%
  left_join(city_health_region, by = "city_residency") %>%
  left_join(immunity_mean_region, by = c("code_health_region", "year"))

# Fit splines on FULL data so knots are consistent across all years
dengue_temp <- dengue_temp %>%
  mutate(
    temp_bs_lag11 = bs(mean_2m_air_temp_degree1_lag1, df = 4)[,1],
    temp_bs_lag12 = bs(mean_2m_air_temp_degree1_lag1, df = 4)[,2],
    temp_bs_lag13 = bs(mean_2m_air_temp_degree1_lag1, df = 4)[,3],
    temp_bs_lag14 = bs(mean_2m_air_temp_degree1_lag1, df = 4)[,4],
    temp_bs_lag21 = bs(mean_2m_air_temp_degree1_lag2, df = 4)[,1],
    temp_bs_lag22 = bs(mean_2m_air_temp_degree1_lag2, df = 4)[,2],
    temp_bs_lag23 = bs(mean_2m_air_temp_degree1_lag2, df = 4)[,3],
    temp_bs_lag24 = bs(mean_2m_air_temp_degree1_lag2, df = 4)[,4],
    temp_bs_lag31 = bs(mean_2m_air_temp_degree1_lag3, df = 4)[,1],
    temp_bs_lag32 = bs(mean_2m_air_temp_degree1_lag3, df = 4)[,2],
    temp_bs_lag33 = bs(mean_2m_air_temp_degree1_lag3, df = 4)[,3],
    temp_bs_lag34 = bs(mean_2m_air_temp_degree1_lag3, df = 4)[,4],
    temp_bs_lag41 = bs(mean_2m_air_temp_degree1_lag4, df = 4)[,1],
    temp_bs_lag42 = bs(mean_2m_air_temp_degree1_lag4, df = 4)[,2],
    temp_bs_lag43 = bs(mean_2m_air_temp_degree1_lag4, df = 4)[,3],
    temp_bs_lag44 = bs(mean_2m_air_temp_degree1_lag4, df = 4)[,4],
    temp_bs_lag51 = bs(mean_2m_air_temp_degree1_lag5, df = 4)[,1],
    temp_bs_lag52 = bs(mean_2m_air_temp_degree1_lag5, df = 4)[,2],
    temp_bs_lag53 = bs(mean_2m_air_temp_degree1_lag5, df = 4)[,3],
    temp_bs_lag54 = bs(mean_2m_air_temp_degree1_lag5, df = 4)[,4],
    precip_ns_lag11 = ns(total_precipitation_lag1, df = 3)[,1],
    precip_ns_lag12 = ns(total_precipitation_lag1, df = 3)[,2],
    precip_ns_lag13 = ns(total_precipitation_lag1, df = 3)[,3],
    precip_ns_lag21 = ns(total_precipitation_lag2, df = 3)[,1],
    precip_ns_lag22 = ns(total_precipitation_lag2, df = 3)[,2],
    precip_ns_lag23 = ns(total_precipitation_lag2, df = 3)[,3],
    precip_ns_lag31 = ns(total_precipitation_lag3, df = 3)[,1],
    precip_ns_lag32 = ns(total_precipitation_lag3, df = 3)[,2],
    precip_ns_lag33 = ns(total_precipitation_lag3, df = 3)[,3],
    precip_ns_lag41 = ns(total_precipitation_lag4, df = 3)[,1],
    precip_ns_lag42 = ns(total_precipitation_lag4, df = 3)[,2],
    precip_ns_lag43 = ns(total_precipitation_lag4, df = 3)[,3],
    precip_ns_lag51 = ns(total_precipitation_lag5, df = 3)[,1],
    precip_ns_lag52 = ns(total_precipitation_lag5, df = 3)[,2],
    precip_ns_lag53 = ns(total_precipitation_lag5, df = 3)[,3]
  )

cat("Data prepared with splines.\n")

# ---------------------------------------------------------------
# 5. Helper: build formula
# ---------------------------------------------------------------
build_formula <- function(vars) {
  fixed_effects <- intersect(vars, c("city_residency", "month", "year",
                                     "month_childs", "year_childs",
                                     "year_region", "month_region"))
  predictors    <- setdiff(vars, fixed_effects)

  if (length(predictors) > 0 && length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "),
                         "|", paste(fixed_effects, collapse = " + "))
  } else if (length(predictors) > 0) {
    formula_str <- paste("n_cases ~", paste(predictors, collapse = " + "))
  } else if (length(fixed_effects) > 0) {
    formula_str <- paste("n_cases ~ 1 |", paste(fixed_effects, collapse = " + "))
  } else {
    formula_str <- "n_cases ~ 1"
  }
  as.formula(formula_str)
}

# ---------------------------------------------------------------
# 6. Loop through years and fit model
# ---------------------------------------------------------------
all_years     <- sort(unique(dengue_temp$year))
model_formula <- build_formula(model_vars)

cat(sprintf("Years found in data (%d): %s\n", length(all_years), paste(all_years, collapse = ", ")))
cat("Model formula:", deparse(model_formula), "\n\n")

all_results <- list()

for (yr in all_years) {
  cat(sprintf("Fitting year %d...\n", yr))

  dengue_year <- dengue_temp %>% filter(year == yr)

  result <- tryCatch({
    model <- fixest::fepois(model_formula,
                            offset   = ~log_pop_offset,
                            data     = dengue_year,
                            nthreads = 1)
    coef(model)
  }, error = function(e) {
    cat(sprintf("  ERROR for year %d: %s\n", yr, e$message))
    NULL
  })

  if (!is.null(result)) {
    result_df            <- as.data.frame(t(result))
    result_df$model_name <- model_name
    result_df$year       <- yr
    all_results[[length(all_results) + 1]] <- result_df
    cat(sprintf("  Done.\n"))
  }
}

# ---------------------------------------------------------------
# 7. Combine and save
# ---------------------------------------------------------------
all_coefs <- bind_rows(all_results)

out_path <- file.path(output_dir, "all_coefs_by_year.csv")
write_csv(all_coefs, out_path)

cat(sprintf("\nSaved coefficients for %d years to: %s\n", nrow(all_coefs), out_path))
