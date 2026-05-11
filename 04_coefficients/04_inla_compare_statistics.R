# run_model.R
# Driven by SLURM_ARRAY_TASK_ID (1-indexed).
# Each array task fits one model specification and extracts posterior samples
# of the temperature B-spline fixed effects.

# ── Packages ───────────────────────────────────────────────────────────────────
library(tidyverse)
library(magrittr)
library(splines)
library(INLA)

# ── Model specifications ───────────────────────────────────────────────────────
models_with_temp <- list(
  # "Year + Month" = c(
  #   "month", "year", "city_residency"
  # ),
  # "Year|region + Month|region" = c(
  #   "month_region", "year_region", "city_residency"
  # ),
  # "Climate(lag 1-5)" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency"
  # ),
  # "Climate(lag 1-5) + Year|region" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","year_region","city_residency"
  # ),
  # "Climate(lag 1-5) + Month|region" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","month_region","city_residency"
  # ),
  # "Year + Month + P(Climate)" = c(
  #   "mean_2m_air_temp_degree1_lag1","mean_2m_air_temp_degree1_lag2",
  #   "mean_2m_air_temp_degree1_lag3","mean_2m_air_temp_degree2_lag1",
  #   "mean_2m_air_temp_degree2_lag2","mean_2m_air_temp_degree2_lag3",
  #   "mean_2m_air_temp_degree3_lag1","mean_2m_air_temp_degree3_lag2",
  #   "mean_2m_air_temp_degree3_lag3","total_precipitation_lag1",
  #   "total_precipitation_lag2","total_precipitation_lag3",
  #   "city_residency","month","year"
  # ),
  # "Climate(lag 1) + Year|region + Month|region" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "total_precipitation_lag1","month_region","year_region","city_residency"
  # ),
  # "Climate(lag 1-3) + Year|region + Month|region" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","month_region","year_region","city_residency"
  # ),
  # "Climate(lag 1-5) + Year|region + Month|region" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","month_region","year_region","city_residency"
  # ),
  # "Climate(lag 1-5) + Year|region + Month|region + Socio" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","month_region","year_region","city_residency",
  #   "urban_area_ha","gdp_per_capita","births"
  # ),
  # "Climate(lag 1-5) + Year|region + Month|region + PriorCases" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","month_region","year_region","city_residency",
  #   "immunity_lag1","immunity_lag2","immunity_lag3"
  # ),
  # "Climate(lag 1-5) + Year|region + Month|region + SeroRepla" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","month_region","year_region","city_residency",
  #   "serotype_replacement"
  # ),
  # "Climate(lag 1-5) + Year|region + Month|region + Immunity" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","month_region","year_region","city_residency",
  #   "Pr_0priorinf","Pr_1priorinf"
  # ),
  # "Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency",
  #   "immunity_lag1","immunity_lag2","immunity_lag3",
  #   "serotype_replacement","urban_area_ha","gdp_per_capita","births",
  #   "Pr_0priorinf","Pr_1priorinf","year_region","month_region"
  # ),
  "Climate(lag 1-5) + Socio" = c(
    "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
    "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
    "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
    "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
    "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
    "total_precipitation_lag1","total_precipitation_lag2",
    "total_precipitation_lag3","total_precipitation_lag4",
    "total_precipitation_lag5","city_residency",
    "urban_area_ha","gdp_per_capita","births"
  )
  # "Climate(lag 1-5) + PriorCases" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency",
  #   "immunity_lag1","immunity_lag2","immunity_lag3"
  # ),
  # "Climate(lag 1-5) + SeroRepla" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency","serotype_replacement"
  # ),
  # "Climate(lag 1-5) + Immunity" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency",
  #   "Pr_0priorinf","Pr_1priorinf"
  # ),
  # "Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency",
  #   "immunity_lag1","immunity_lag2","immunity_lag3",
  #   "serotype_replacement","urban_area_ha","gdp_per_capita","births",
  #   "Pr_0priorinf","Pr_1priorinf"
  # ),
  # "Month|region + PriorCases + SeroRepla + Socio + Immunity" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5","city_residency",
  #   "immunity_lag1","immunity_lag2","immunity_lag3",
  #   "serotype_replacement","urban_area_ha","gdp_per_capita","births",
  #   "Pr_0priorinf","Pr_1priorinf","month_region"
  # ),
  # "Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity (no city)" = c(
  #   "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  #   "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  #   "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  #   "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  #   "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54",
  #   "total_precipitation_lag1","total_precipitation_lag2",
  #   "total_precipitation_lag3","total_precipitation_lag4",
  #   "total_precipitation_lag5",
  #   "immunity_lag1","immunity_lag2","immunity_lag3",
  #   "serotype_replacement","urban_area_ha","gdp_per_capita","births",
  #   "Pr_0priorinf","Pr_1priorinf","year_region","month_region"
  # )
)

# ── Identify this task ─────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else
                                  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("No task ID supplied. Pass it as a command-line argument (Rscript run_model.R <id>) or via SLURM_ARRAY_TASK_ID.")

model_names  <- names(models_with_temp)
n_models     <- length(model_names)
cat(sprintf("Total models: %d  |  Running task: %d\n", n_models, task_id))
if (task_id < 1 || task_id > n_models) {
  stop(sprintf("task_id %d is out of range [1, %d].", task_id, n_models))
}

model_name <- model_names[[task_id]]
model_vars <- models_with_temp[[task_id]]

model_tag <- gsub("[^A-Za-z0-9]+", "_", model_name)
model_tag <- gsub("_+$", "", model_tag)

cat(sprintf("Model %d: %s\n", task_id, model_name))

# ── Output directory ───────────────────────────────────────────────────────────
out_dir <- file.path("results", sprintf("%02d_%s", task_id, model_tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── All possible temp B-spline column names (lags 1-5, basis 1-4) ─────────────
all_temp_bs_names <- paste0(
  "temp_bs_lag",
  rep(1:5, each = 4),
  rep(1:4, times = 5)
)

# All possible precipitation variable names (lags 1-5)
all_precip_names <- paste0("total_precipitation_lag", 1:5)

# Which temp spline / precip terms does this model actually include?
temp_coef_names   <- intersect(all_temp_bs_names, model_vars)
precip_coef_names <- intersect(all_precip_names,  model_vars)

# ── Graceful exit for models with no temperature spline terms ──────────────────
if (length(temp_coef_names) == 0) {
  cat("No temperature B-spline terms in this model — writing empty output.\n")
  write_csv(
    tibble(note = "No temperature B-spline terms in this model specification."),
    file.path(out_dir, sprintf("%02d_%s_fe_samples_temp.csv", task_id, model_tag))
  )
  quit(status = 0)
}

# ── Load data ──────────────────────────────────────────────────────────────────
dengue_temp <- read_csv("data/model_input_brazil_immunity_city_with_priorinf.csv",
                        show_col_types = FALSE)

# ── B-spline basis functions ───────────────────────────────────────────────────
temp_bs_fitted_lag1 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag1, df = 4)
temp_bs_fitted_lag2 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag2, df = 4)
temp_bs_fitted_lag3 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag3, df = 4)
temp_bs_fitted_lag4 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag4, df = 4)
temp_bs_fitted_lag5 <- bs(dengue_temp$mean_2m_air_temp_degree1_lag5, df = 4)

dengue_temp <- dengue_temp %>%
  mutate(
    temp_bs_lag11 = temp_bs_fitted_lag1[, 1],
    temp_bs_lag12 = temp_bs_fitted_lag1[, 2],
    temp_bs_lag13 = temp_bs_fitted_lag1[, 3],
    temp_bs_lag14 = temp_bs_fitted_lag1[, 4],
    temp_bs_lag21 = temp_bs_fitted_lag2[, 1],
    temp_bs_lag22 = temp_bs_fitted_lag2[, 2],
    temp_bs_lag23 = temp_bs_fitted_lag2[, 3],
    temp_bs_lag24 = temp_bs_fitted_lag2[, 4],
    temp_bs_lag31 = temp_bs_fitted_lag3[, 1],
    temp_bs_lag32 = temp_bs_fitted_lag3[, 2],
    temp_bs_lag33 = temp_bs_fitted_lag3[, 3],
    temp_bs_lag34 = temp_bs_fitted_lag3[, 4],
    temp_bs_lag41 = temp_bs_fitted_lag4[, 1],
    temp_bs_lag42 = temp_bs_fitted_lag4[, 2],
    temp_bs_lag43 = temp_bs_fitted_lag4[, 3],
    temp_bs_lag44 = temp_bs_fitted_lag4[, 4],
    temp_bs_lag51 = temp_bs_fitted_lag5[, 1],
    temp_bs_lag52 = temp_bs_fitted_lag5[, 2],
    temp_bs_lag53 = temp_bs_fitted_lag5[, 3],
    temp_bs_lag54 = temp_bs_fitted_lag5[, 4],
    log_pop_offset = log(population / 100000),
    city_id        = as.integer(as.factor(city_residency)),
    year_id        = as.integer(as.factor(year_region)),
    month_id       = as.integer(as.factor(month_region))
  )

# ── Build INLA formula dynamically ────────────────────────────────────────────
random_effect_vars <- c("city_residency", "year_region", "month_region",
                        "year", "month")

re_index_map <- list(
  city_residency = "city_id",
  year_region    = "year_id",
  month_region   = "month_id",
  year           = "year_id",
  month          = "month_id"
)

cyclic_rw1_vars <- c("month_region", "month")

fe_vars <- setdiff(model_vars, random_effect_vars)
re_vars <- intersect(model_vars, random_effect_vars)

if (length(fe_vars) > 0) {
  fe_str <- paste(fe_vars, collapse = " + ")
} else {
  fe_str <- NULL
}

re_str_parts <- sapply(re_vars, function(v) {
  idx_col <- re_index_map[[v]]
  if (v %in% cyclic_rw1_vars) {
    sprintf('f(%s, model = "rw1", cyclic = TRUE)', idx_col)
  } else {
    sprintf('f(%s, model = "iid")', idx_col)
  }
})
re_str <- if (length(re_str_parts) > 0) paste(re_str_parts, collapse = " + ") else NULL

rhs_parts <- c(fe_str, re_str)
rhs_parts <- rhs_parts[!sapply(rhs_parts, is.null)]
rhs       <- paste(rhs_parts, collapse = " + ")

formula_inla <- as.formula(paste("n_cases ~ -1 +", rhs))
cat("Formula:\n"); print(formula_inla)

# ── Fit model ──────────────────────────────────────────────────────────────────
model_inla <- inla(
  formula_inla,
  family  = "poisson",
  data    = dengue_temp,
  offset  = dengue_temp$log_pop_offset,
  control.predictor = list(compute = TRUE, link = 1),
  control.compute   = list(
    dic    = TRUE,
    waic   = TRUE,
    cpo    = TRUE,
    config = TRUE
  ),
  control.inla = list(
    strategy     = "simplified.laplace",
    int.strategy = "ccd"
  ),
  verbose = FALSE
)

# ── Export fixed effects summary ───────────────────────────────────────────────
write_csv(
  model_inla$summary.fixed %>% rownames_to_column("parameter"),
  file.path(out_dir, sprintf("%02d_%s_fixed_effects_summary.csv", task_id, model_tag))
)

# ── Posterior sampling — temperature B-spline + precipitation coefficients ─────
cat("Sampling from posterior...\n")
n_samples <- 2000
samples <- inla.posterior.sample(n_samples, model_inla)

# Verify naming used by INLA (fixed effects get a ':1' suffix in the latent field)
all_latent_names <- rownames(samples[[1]]$latent)
cat("Temperature B-spline terms found in latent field:\n")
print(all_latent_names[grepl("temp_bs", all_latent_names)])
cat("Precipitation terms found in latent field:\n")
print(all_latent_names[grepl("total_precipitation", all_latent_names)])

# Helper to extract named coefficients from a single posterior sample
extract_coefs <- function(s, coef_names) {
  nm <- rownames(s$latent)
  sapply(coef_names, function(coef) {
    idx <- which(nm == paste0(coef, ":1"))
    if (length(idx) == 1) s$latent[idx, 1] else NA_real_
  })
}

# Temperature B-spline samples
fe_samples_temp <- t(sapply(samples, extract_coefs, coef_names = temp_coef_names))
colnames(fe_samples_temp) <- temp_coef_names

write_csv(
  as.data.frame(fe_samples_temp),
  file.path(out_dir, sprintf("%02d_%s_fe_samples_temp.csv", task_id, model_tag))
)
cat(sprintf("Temperature samples written (%d terms).\n", length(temp_coef_names)))

# Precipitation samples (only written if this model includes precipitation terms)
if (length(precip_coef_names) > 0) {
  fe_samples_precip <- t(sapply(samples, extract_coefs, coef_names = precip_coef_names))
  
  # When there is only one precip term, sapply returns a vector and t() gives a
  # 1-row matrix — coerce explicitly so colnames() always works.
  fe_samples_precip <- matrix(fe_samples_precip,
                               nrow  = n_samples,
                               ncol  = length(precip_coef_names),
                               dimnames = list(NULL, precip_coef_names))

  write_csv(
    as.data.frame(fe_samples_precip),
    file.path(out_dir, sprintf("%02d_%s_fe_samples_precip.csv", task_id, model_tag))
  )
  cat(sprintf("Precipitation samples written (%d terms).\n", length(precip_coef_names)))
}

cat(sprintf("Done. Outputs written to: %s\n", out_dir))