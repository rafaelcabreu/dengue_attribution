# test_model06.R
# Runs model: Climate(lag 1-5) + Year|region + Month|region
# Temperature entered as B-spline basis (4 basis functions x 5 lags)
# Precipitation lags 1-5

library(tidyverse)
library(magrittr)
library(splines)
library(INLA)

# ── Fitting strategies (tried in order) ───────────────────────────────────────
strategies <- list(
  list(label        = "adaptive + grid (original)",
       strategy     = "adaptive",
       int.strategy = "grid",
       num.threads  = NULL),
  list(label        = "simplified.laplace + ccd",
       strategy     = "simplified.laplace",
       int.strategy = "ccd",
       num.threads  = NULL),
  list(label        = "simplified.laplace + ccd + single thread",
       strategy     = "simplified.laplace",
       int.strategy = "ccd",
       num.threads  = 1)
)

# ── Load and prepare data ──────────────────────────────────────────────────────
cat("Loading data...\n")
dengue_temp <- read_csv("data/model_input_mexico_immunity_city.csv",
                        show_col_types = FALSE)

# ── B-spline basis functions ───────────────────────────────────────────────────
# Build splines for every lag that appears in *any* model so the column set is
# stable; we only use the subset relevant to this model in the formula.
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
    year_id        = as.integer(as.factor(year)),
    month_id       = as.integer(as.factor(month))
  )

# ── Sanity checks ──────────────────────────────────────────────────────────────
cat("Checking random effect indices...\n")
stopifnot(!anyNA(dengue_temp$city_id),  min(dengue_temp$city_id)  >= 1)
stopifnot(!anyNA(dengue_temp$year_id),  min(dengue_temp$year_id)  >= 1)
stopifnot(!anyNA(dengue_temp$month_id), min(dengue_temp$month_id) >= 1)
cat(sprintf("  city_id:  %d levels\n", n_distinct(dengue_temp$city_id)))
cat(sprintf("  year_id:  %d levels\n", n_distinct(dengue_temp$year_id)))
cat(sprintf("  month_id: %d levels\n", n_distinct(dengue_temp$month_id)))

# ── Define variable names ──────────────────────────────────────────────────────
temp_bs_names <- c(
  "temp_bs_lag11","temp_bs_lag12","temp_bs_lag13","temp_bs_lag14",
  "temp_bs_lag21","temp_bs_lag22","temp_bs_lag23","temp_bs_lag24",
  "temp_bs_lag31","temp_bs_lag32","temp_bs_lag33","temp_bs_lag34",
  "temp_bs_lag41","temp_bs_lag42","temp_bs_lag43","temp_bs_lag44",
  "temp_bs_lag51","temp_bs_lag52","temp_bs_lag53","temp_bs_lag54"
)

precip_names <- c(
  "total_precipitation_lag1","total_precipitation_lag2",
  "total_precipitation_lag3","total_precipitation_lag4",
  "total_precipitation_lag5"
)

# All continuous fixed effects to scale
vars_to_scale <- intersect(
  c(temp_bs_names, precip_names),
  names(dengue_temp)
)

# ── Scale continuous fixed effects ─────────────────────────────────────────────
cat("Scaling continuous fixed effects...\n")
scaling_params <- list()
dengue_temp <- dengue_temp %>%
  mutate(across(
    all_of(vars_to_scale),
    ~ {
      m <- mean(.x, na.rm = TRUE)
      s <- sd(.x,   na.rm = TRUE)
      scaling_params[[cur_column()]] <<- list(mean = m, sd = s)
      if (s > 0) (.x - m) / s else .x - m
    }
  ))

cat("Scaling params:\n")
for (v in vars_to_scale) {
  cat(sprintf("  %-35s  mean=%.3g  sd=%.3g\n", v,
              scaling_params[[v]]$mean, scaling_params[[v]]$sd))
}

# ── Formula ────────────────────────────────────────────────────────────────────
formula_m06 <- n_cases ~ -1 +
  # Temperature B-spline basis, lag 1
  temp_bs_lag11 + temp_bs_lag12 + temp_bs_lag13 + temp_bs_lag14 +
  # Temperature B-spline basis, lag 2
  temp_bs_lag21 + temp_bs_lag22 + temp_bs_lag23 + temp_bs_lag24 +
  # Temperature B-spline basis, lag 3
  temp_bs_lag31 + temp_bs_lag32 + temp_bs_lag33 + temp_bs_lag34 +
  # Temperature B-spline basis, lag 4
  temp_bs_lag41 + temp_bs_lag42 + temp_bs_lag43 + temp_bs_lag44 +
  # Temperature B-spline basis, lag 5
  temp_bs_lag51 + temp_bs_lag52 + temp_bs_lag53 + temp_bs_lag54 +
  # Precipitation lags 1-5
  total_precipitation_lag1 + total_precipitation_lag2 +
  total_precipitation_lag3 + total_precipitation_lag4 +
  total_precipitation_lag5 +
  # Random effects
  f(city_id,  model = "iid") +
  f(year_id,  model = "iid") +
  f(month_id, model = "rw1", cyclic = TRUE)

cat("Formula:\n"); print(formula_m06)

# ── Fit with progressive fallback ─────────────────────────────────────────────
model_inla <- NULL
for (strat in strategies) {
  cat(sprintf("\nTrying strategy: %s\n", strat$label))
  result <- tryCatch({
    inla_args <- list(
      formula           = formula_m06,
      family            = "poisson",
      data              = dengue_temp,
      offset            = dengue_temp$log_pop_offset,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
      control.inla      = list(strategy     = strat$strategy,
                               int.strategy = strat$int.strategy),
      verbose           = FALSE
    )
    if (!is.null(strat$num.threads)) inla_args$num.threads <- strat$num.threads
    do.call(inla, inla_args)
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", conditionMessage(e)))
    NULL
  })

  if (!is.null(result)) {
    cat(sprintf("  SUCCESS with strategy: %s\n", strat$label))
    model_inla <- result
    break
  }
}

if (is.null(model_inla)) stop("All strategies failed for model 06.")

# ── Results summary ────────────────────────────────────────────────────────────
cat(sprintf("\nDIC:  %.2f\n", model_inla$dic$dic))
cat(sprintf("WAIC: %.2f\n", model_inla$waic$waic))

cat("\nFixed effects summary:\n")
print(model_inla$summary.fixed)

# ── Save outputs ───────────────────────────────────────────────────────────────
out_dir <- "results/mexico"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  model_inla$summary.fixed %>% rownames_to_column("parameter"),
  file.path(out_dir, "06_fixed_effects_summary.csv")
)
saveRDS(scaling_params, file.path(out_dir, "06_scaling_params.rds"))

# ── Helper: extract named coefficients into a correctly-shaped matrix ──────────
# Handles the single-term edge case where sapply() returns a vector instead of
# a matrix, causing colnames<- to fail.
extract_coef_matrix <- function(samples, coef_names, n_samples = 2000) {
  mat <- t(sapply(samples, function(s) {
    nm <- rownames(s$latent)
    sapply(coef_names, function(coef) {
      idx <- which(nm == paste0(coef, ":1"))
      if (length(idx) == 1) s$latent[idx, 1] else NA_real_
    })
  }))
  matrix(mat,
         nrow     = n_samples,
         ncol     = length(coef_names),
         dimnames = list(NULL, coef_names))
}

# ── Helper: rescale posterior samples back to original predictor units ─────────
# β_original = β_scaled / sd(x)
# Recovers the log-rate change per one-unit increase in the original predictor.
# Terms absent from scaling_params pass through unchanged.
rescale_coef_matrix <- function(coef_mat, coef_names, scaling_params) {
  for (j in seq_along(coef_names)) {
    nm <- coef_names[j]
    if (!is.null(scaling_params[[nm]]) && scaling_params[[nm]]$sd > 0) {
      coef_mat[, j] <- coef_mat[, j] / scaling_params[[nm]]$sd
    }
  }
  coef_mat
}

# ── Helper: build and write a scaled-vs-rescaled summary table ────────────────
write_summary_comparison <- function(coef_names, samples_scaled, samples_original,
                                     scaling_params, path) {
  tibble(
    parameter        = coef_names,
    sd_original      = sapply(coef_names, function(v) scaling_params[[v]]$sd),
    mean_scaled      = colMeans(samples_scaled,   na.rm = TRUE),
    mean_original    = colMeans(samples_original, na.rm = TRUE),
    sd_scaled        = apply(samples_scaled,   2, sd, na.rm = TRUE),
    sd_original_coef = apply(samples_original, 2, sd, na.rm = TRUE)
  ) %>% write_csv(path)
}

# ── Posterior samples ──────────────────────────────────────────────────────────
# Subset to only those variables actually present in the data
temp_coef_names   <- intersect(temp_bs_names, vars_to_scale)
precip_coef_names <- intersect(precip_names,  vars_to_scale)

cat("\nSampling 2000 posterior draws...\n")
n_samples <- 2000
samples   <- inla.posterior.sample(n_samples, model_inla)

all_latent_names <- rownames(samples[[1]]$latent)
cat("Temperature B-spline terms found in latent field:\n")
print(all_latent_names[grepl("temp_bs", all_latent_names)])
cat("Precipitation terms found in latent field:\n")
print(all_latent_names[grepl("precipitation", all_latent_names)])

# ── Temperature posterior samples ─────────────────────────────────────────────
found <- sapply(temp_coef_names, function(cn) paste0(cn, ":1") %in% all_latent_names)
if (any(!found))
  warning(sprintf("Temp terms not found in latent field: %s",
                  paste(temp_coef_names[!found], collapse = ", ")))

fe_samples_temp_scaled   <- extract_coef_matrix(samples, temp_coef_names, n_samples)
fe_samples_temp_original <- rescale_coef_matrix(fe_samples_temp_scaled,
                                                 temp_coef_names, scaling_params)

write_csv(as.data.frame(fe_samples_temp_scaled),
          file.path(out_dir, "06_fe_samples_temp_scaled.csv"))
write_csv(as.data.frame(fe_samples_temp_original),
          file.path(out_dir, "06_fe_samples_temp_original_scale.csv"))
write_summary_comparison(
  temp_coef_names,
  fe_samples_temp_scaled, fe_samples_temp_original,
  scaling_params,
  file.path(out_dir, "06_fe_summary_temp_scaled_vs_original.csv")
)
cat(sprintf("  Temperature samples written (%d terms, scaled + rescaled).\n",
            length(temp_coef_names)))

# ── Precipitation posterior samples ───────────────────────────────────────────
found <- sapply(precip_coef_names, function(cn) paste0(cn, ":1") %in% all_latent_names)
if (any(!found))
  warning(sprintf("Precip terms not found in latent field: %s",
                  paste(precip_coef_names[!found], collapse = ", ")))

fe_samples_precip_scaled   <- extract_coef_matrix(samples, precip_coef_names, n_samples)
fe_samples_precip_original <- rescale_coef_matrix(fe_samples_precip_scaled,
                                                   precip_coef_names, scaling_params)

write_csv(as.data.frame(fe_samples_precip_scaled),
          file.path(out_dir, "06_fe_samples_precip_scaled.csv"))
write_csv(as.data.frame(fe_samples_precip_original),
          file.path(out_dir, "06_fe_samples_precip_original_scale.csv"))
write_summary_comparison(
  precip_coef_names,
  fe_samples_precip_scaled, fe_samples_precip_original,
  scaling_params,
  file.path(out_dir, "06_fe_summary_precip_scaled_vs_original.csv")
)
cat(sprintf("  Precipitation samples written (%d terms, scaled + rescaled).\n",
            length(precip_coef_names)))

cat(sprintf("\nDone. All outputs written to: %s\n", out_dir))