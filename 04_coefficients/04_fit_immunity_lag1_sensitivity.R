# 04_fit_immunity_lag1_sensitivity.R
#
# Sensitivity check: does lagging the posterior immunity variables
# (Pr_0priorinf / Pr_1priorinf / Pr_2priorinf) by ONE YEAR change their
# estimated coefficients, relative to the flagship model's current-year
# immunity ("Climate(lag 1-5) + Year|region + Month|region + PriorCases +
# SeroRepla + Socio + Immunity", the same formula as
# 02_statistics/02_fit_model_incidence.R)?
#
# Fits BOTH versions as single (non-bootstrap) point-estimate models — same
# style as 02_fit_model_incidence.R, not a full 1000-iteration bootstrap,
# since this is a quick "how much do the coefficients move" check, not a
# publication-grade CI — on the IDENTICAL set of rows, so the comparison
# isn't confounded by the extra row loss the lag requires (each city's
# first data year has no prior year to lag from). Writes a side-by-side
# coefficient table, actual (current-year) vs. lag-1-year immunity.
#
# NOTE: "immunity" here means the posterior serological immunity
# (Pr_0priorinf/Pr_1priorinf/Pr_2priorinf) built by
# 01_model_fit/build_immunity_city_draws.R — NOT the pre-existing
# immunity_lag1/immunity_lag2/immunity_lag3 columns already in the base
# CSV (a different, unrelated covariate carried through unchanged in both
# models below).
#
# Run with dengue_attribution/ as the working directory.

library(tidyverse)
library(magrittr)
library(splines)
library(fixest)

# ---------------------------------------------------------------
# 1. Load data and BOTH immunity joins (current-year + lag-1-year)
# ---------------------------------------------------------------
dengue_temp <- read_csv("../data/model_input_brazil_immunity_city_with_priorinf.csv", show_col_types = FALSE)
dengue_temp$log_pop_offset <- log(dengue_temp$population / 100000)

city_health_region   <- read_csv("../data/code_health_region_to_city_residency.csv", show_col_types = FALSE)
immunity_mean_region <- read_csv("../data/immunity_posterior_mean_region.csv", show_col_types = FALSE)

# Immunity lagged by one year: a municipality's posterior-mean immunity from
# YEAR-1, attached to the row for YEAR (shift the year key forward by one).
immunity_lag1yr <- immunity_mean_region %>%
  mutate(year = year + 1) %>%
  rename(
    Pr_0priorinf_lag1 = Pr_0priorinf,
    Pr_1priorinf_lag1 = Pr_1priorinf,
    Pr_2priorinf_lag1 = Pr_2priorinf
  )

dengue_temp <- dengue_temp %>%
  select(-any_of(c("Pr_0priorinf", "Pr_1priorinf", "Pr_2priorinf"))) %>%
  left_join(city_health_region, by = "city_residency") %>%
  left_join(immunity_mean_region, by = c("code_health_region", "year")) %>%
  left_join(immunity_lag1yr, by = c("code_health_region", "year"))

# ---------------------------------------------------------------
# 2. Restrict to rows where BOTH current-year and lag-1-year immunity are
#    available, so the two models are fit on the identical sample. Some
#    loss is expected: each municipality's first data year has no prior
#    year to lag from, plus 12 code_health_region values are entirely
#    absent from immunity_posterior_mean_region.csv for every year (found
#    while testing 04_coefficients/inla_model_tests/test_models_6_14_15.R —
#    affects both joins identically since it's a total gap, not a
#    year-specific one).
# ---------------------------------------------------------------
n_before <- nrow(dengue_temp)
dengue_temp <- dengue_temp %>%
  filter(!is.na(Pr_0priorinf), !is.na(Pr_1priorinf),
         !is.na(Pr_0priorinf_lag1), !is.na(Pr_1priorinf_lag1))
cat(sprintf("Kept %d of %d rows with both current-year and lag-1-year immunity available (%.1f%% dropped).\n",
            nrow(dengue_temp), n_before, 100 * (1 - nrow(dengue_temp) / n_before)))

# ---------------------------------------------------------------
# 3. B-splines on temperature (identical for both fits)
# ---------------------------------------------------------------
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

dengue_temp <- dengue_temp %>%
  mutate(
    temp_bs_lag11 = temp_bs_fitted_lag1[, 1], temp_bs_lag12 = temp_bs_fitted_lag1[, 2],
    temp_bs_lag13 = temp_bs_fitted_lag1[, 3], temp_bs_lag14 = temp_bs_fitted_lag1[, 4],
    temp_bs_lag21 = temp_bs_fitted_lag2[, 1], temp_bs_lag22 = temp_bs_fitted_lag2[, 2],
    temp_bs_lag23 = temp_bs_fitted_lag2[, 3], temp_bs_lag24 = temp_bs_fitted_lag2[, 4],
    temp_bs_lag31 = temp_bs_fitted_lag3[, 1], temp_bs_lag32 = temp_bs_fitted_lag3[, 2],
    temp_bs_lag33 = temp_bs_fitted_lag3[, 3], temp_bs_lag34 = temp_bs_fitted_lag3[, 4],
    temp_bs_lag41 = temp_bs_fitted_lag4[, 1], temp_bs_lag42 = temp_bs_fitted_lag4[, 2],
    temp_bs_lag43 = temp_bs_fitted_lag4[, 3], temp_bs_lag44 = temp_bs_fitted_lag4[, 4],
    temp_bs_lag51 = temp_bs_fitted_lag5[, 1], temp_bs_lag52 = temp_bs_fitted_lag5[, 2],
    temp_bs_lag53 = temp_bs_fitted_lag5[, 3], temp_bs_lag54 = temp_bs_fitted_lag5[, 4]
  )

# ---------------------------------------------------------------
# 4. Fit the flagship model as-is ("actual") and the lag-1-year-immunity
#    version — identical formula, only Pr_0priorinf/Pr_1priorinf/
#    Pr_2priorinf swapped for their _lag1 counterparts.
# ---------------------------------------------------------------
shared_rhs <- paste(
  "temp_bs_lag11 + temp_bs_lag12 + temp_bs_lag13 + temp_bs_lag14 +",
  "temp_bs_lag21 + temp_bs_lag22 + temp_bs_lag23 + temp_bs_lag24 +",
  "temp_bs_lag31 + temp_bs_lag32 + temp_bs_lag33 + temp_bs_lag34 +",
  "temp_bs_lag41 + temp_bs_lag42 + temp_bs_lag43 + temp_bs_lag44 +",
  "temp_bs_lag51 + temp_bs_lag52 + temp_bs_lag53 + temp_bs_lag54 +",
  "total_precipitation_lag1 + total_precipitation_lag2 + total_precipitation_lag3 +",
  "total_precipitation_lag4 + total_precipitation_lag5 +",
  "immunity_lag1 + immunity_lag2 + immunity_lag3 +",
  "serotype_replacement + urban_area_ha + gdp_per_capita + births"
)

formula_actual <- as.formula(paste(
  "n_cases ~", shared_rhs, "+ Pr_0priorinf + Pr_1priorinf + Pr_2priorinf",
  "| city_residency + year_region + month_region"
))

formula_lag1 <- as.formula(paste(
  "n_cases ~", shared_rhs, "+ Pr_0priorinf_lag1 + Pr_1priorinf_lag1 + Pr_2priorinf_lag1",
  "| city_residency + year_region + month_region"
))

cat("Fitting actual (current-year immunity) model...\n")
model_actual <- fixest::fepois(formula_actual, offset = ~log_pop_offset, data = dengue_temp, combine.quick = FALSE)

cat("Fitting lag-1-year-immunity model...\n")
model_lag1 <- fixest::fepois(formula_lag1, offset = ~log_pop_offset, data = dengue_temp, combine.quick = FALSE)

print(summary(model_actual))
print(summary(model_lag1))

# ---------------------------------------------------------------
# 5. Side-by-side coefficient comparison table.
# ---------------------------------------------------------------
extract_coefs <- function(model) {
  ct <- as.data.frame(summary(model)$coeftable)
  ct$term <- rownames(ct)
  rownames(ct) <- NULL
  ct %>% select(term, everything())
}

coefs_actual <- extract_coefs(model_actual) %>% rename_with(~ paste0(.x, "_actual"), -term)
coefs_lag1 <- extract_coefs(model_lag1) %>%
  mutate(term = str_replace(term, "_lag1$", "")) %>%   # align term names (Pr_0priorinf_lag1 -> Pr_0priorinf) for the join
  rename_with(~ paste0(.x, "_lag1"), -term)

comparison <- full_join(coefs_actual, coefs_lag1, by = "term") %>% arrange(term)

immunity_rows <- comparison %>% filter(str_starts(term, "Pr_"))
cat("\n=== Immunity coefficients: actual (current-year) vs. lag-1-year ===\n")
print(immunity_rows)

out_dir <- "04_coefficients/immunity_lag_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(comparison, file.path(out_dir, "coef_comparison_actual_vs_lag1_immunity.csv"))
cat(sprintf("\nFull comparison written to %s\n", file.path(out_dir, "coef_comparison_actual_vs_lag1_immunity.csv")))
