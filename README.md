# Dengue Attribution

Analysis pipeline for attributing dengue fever incidence in Brazil, Mexico, and Peru to anthropogenic climate change. The pipeline fits a suite of Poisson GLMs with B-spline temperature responses, runs a parallelised bootstrap for uncertainty quantification, and produces publication figures.

Code repository: [http://github.com/rafaelcabreu/dengue_attribution](http://github.com/rafaelcabreu/dengue_attribution)

---

## System requirements

### Software dependencies

The recommended approach is Docker, which bundles all dependencies in a reproducible environment. Alternatively, R and Python can be installed locally.

**Docker (recommended):**
- Docker Engine ≥ 20.10 or Docker Desktop ≥ 4.0
- The image is built for `linux/amd64`. On Apple Silicon (M1/M2/M3), Docker Desktop will run it via Rosetta 2 emulation.

**R (≥ 4.5.1) — tested on 4.5.1:**

| Package | Version |
|---|---|
| `tidyverse` | 2.0.0 |
| `fixest` | 0.12.x |
| `splines` | (base R) |
| `glmmTMB` | 1.1.x |
| `MuMIn` | 1.47.x |
| `INLA` | 24.x (installed from r-inla.org) |
| `magrittr` | 2.0.x |

**Python (≥ 3.9) — tested on 3.9:**

| Package | Version |
|---|---|
| `geopandas` | 1.0.1 |
| `matplotlib` | 3.9.2 |
| `numpy` | 1.26.4 |
| `pandas` | 2.2.2 |
| `scipy` | 1.13.1 |
| `seaborn` | 0.13.2 |
| `shapely` | 2.0.6 |
| `statsmodels` | 0.14.3 |
| `climattr` | (installed from GitHub) |

See `requirements.txt` for the full pinned list.

### Operating systems tested

- Linux (Ubuntu 22.04, x86_64) — primary development environment
- macOS (Sonoma 14, x86_64 and Apple Silicon via Docker emulation)
- Windows 10/11 via Docker Desktop (Linux containers)

### Hardware requirements

- RAM: ≥ 8 GB (≥ 16 GB recommended for INLA steps)
- Disk: ≥ 15 GB (Docker image ~5 GB + data archive ~10 GB)
- No GPU required
- For the full bootstrap and sensitivity analysis (thousands of iterations): an HPC cluster with SLURM is required. Single-task local runs are feasible for testing.

---

## Installation guide

### Option A — Docker (recommended)

Docker handles all R and Python dependencies automatically, including INLA.

```bash
# 1. Clone the repository
git clone <repo-url> dengue_attribution
cd dengue_attribution

# 2. Build the Docker image
docker compose build
```

> **Typical install time:** 30–60 minutes on a standard desktop (first build only — INLA compilation dominates). Subsequent builds take 2–5 minutes using the layer cache.

Verify the environment after building:

```bash
docker compose run --rm jupyter python3 tests/test_python_packages.py
docker compose run --rm jupyter Rscript tests/test_r_packages.R
```

### Option B — Local installation

**R packages:**

```r
install.packages(c("tidyverse", "fixest", "splines", "glmmTMB", "MuMIn", "magrittr"))

# INLA requires a separate repository:
install.packages("INLA",
  repos = c(INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)
```

**Python packages:**

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

> **Typical install time:** 15–30 minutes (INLA download and compilation is the slowest step).

---

## Demo

The demo runs step 1 of the pipeline (model fitting on the full dataset) using the provided data. It requires no HPC cluster and completes in ~5 minutes.

### 1. Download the data

The dataset is archived on Mendeley Data:

> **Mendeley Data:** [https://doi.org/10.17632/dyx8r83ktx.1](https://doi.org/10.17632/dyx8r83ktx.1)

Place all extracted files inside the `data/` folder at the root of this repository. All scripts expect data to be located there.

```bash
# If you downloaded a zip archive:
unzip dengue_attribution_data.zip -d data/

# If you downloaded a tar archive:
tar -xf dengue_attribution_data.tar.gz -C data/
```

### 2. Run the demo

```bash
# With Docker (recommended):
docker compose run --rm jupyter Rscript 01_model_fit/fit_model.R

# Or locally:
Rscript 01_model_fit/fit_model.R
```

### 3. Expected output

The script writes three CSV files to `data/statistics/full_model_fits/`:

| File | Contents |
|---|---|
| `all_models_coefficients.csv` | Estimate, SE, z-value, p-value per coefficient per model |
| `all_models_fit_stats.csv` | AIC, BIC, R², pseudo-R² per model |
| `all_models_regional_incidence.csv` | Observed and fitted incidence per 100k by region/year/month |

You should see model fitting progress printed to the console, ending with a confirmation that all three files have been written.

### 4. Expected run time

~5 minutes on a standard desktop (2020-era laptop with 4 cores and 16 GB RAM).

---

## Instructions for use

### Repository structure

```
.
├── 01_model_fit/               # Single-fit: all models on full dataset
├── 02_statistics/              # Parallelised bootstrap for R²/AIC
├── 03_sensitivity/             # Sensitivity analysis (actual vs. counterfactual climate)
├── 04_coefficients/            # Coefficient bootstrap and INLA comparison
├── Figure*.ipynb               # Publication figures (Python/Jupyter)
├── data/                       # Input and output data (not versioned)
├── utils/                      # Shared utilities (e.g. region mapping JSON)
├── tests/                      # Environment verification scripts
├── Dockerfile                  # Reproducible environment
├── docker-compose.yml          # Jupyter service
└── requirements.txt            # Python dependencies
```

### Running the full pipeline

#### Step 1 — Model fitting (`01_model_fit/`)

Fits each model once on the full dataset using `fixest::fepois` (Poisson GLMM with city/year/month fixed effects). Temperature is encoded as B-splines with 4 degrees of freedom across up to 5 monthly lags.

```bash
# Brazil
Rscript 01_model_fit/fit_model.R

# Mexico
Rscript 01_model_fit/fitted_model_mexico.R

# Inside Docker:
docker compose run --rm jupyter Rscript 01_model_fit/fit_model.R
```

**Outputs** (written to `data/statistics/full_model_fits/`):

| File | Contents |
|---|---|
| `all_models_coefficients.csv` | Estimate, SE, z-value, p-value per coefficient per model |
| `all_models_fit_stats.csv` | AIC, BIC, R², pseudo-R² per model |
| `all_models_regional_incidence.csv` | Observed and fitted incidence per 100k by region/year/month |

**Expected runtime:** ~5 minutes per country.

---

#### Step 2 — Bootstrap statistics (`02_statistics/`)

Parallelised R²/AIC bootstrap across 6 models × 1000 iterations. Designed for SLURM (60 tasks), but individual tasks can be run locally.

```bash
# Run a single bootstrap task locally (task 1 of 60):
Rscript 02_statistics/02_bootstrap_parallel_full_sqcorr.R 1

# After all 60 tasks finish, aggregate results:
Rscript 02_statistics/02_aggreagte_full_sqcorr.R
```

Out-of-sample validation (5-year rolling windows):

```bash
# Single validation task (model 1, iterations 1–25):
Rscript 02_statistics/02_fit_outsample_validation.R 1 1 25 data/results/dengue_validation

# Aggregate after all tasks complete:
Rscript 02_statistics/02_aggregate_outsample_full.R
```

**Outputs:**
- `data/statistics/bootstrap_chunks_rsq/rsq_model<MM>_chunk<CCC>.csv` — per-task bootstrap results
- `data/statistics/all_models_bootstrap_rsquared_aic_presampled1000.csv` — combined summary
- `data/results/dengue_validation/all_models_validation_metrics_5year_windows.csv`

---

#### Step 3 — Sensitivity analysis (`03_sensitivity/`)

Refits models under actual (`act`) and counterfactual/natural (`nat`) climate to quantify the attribution of dengue burden to anthropogenic climate change. Requires the counterfactual ensemble files (`data/predict-all/ext-ens*.csv`, 525 files per variant).

```bash
# Run a single task locally — actual climate (task 1 of 1700):
Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_act.R 1

# Natural/counterfactual climate:
Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_nat.R 1

# Aggregate after all tasks complete:
Rscript 03_sensitivity/03_aggregate_results_act.R
Rscript 03_sensitivity/03_aggregate_results_nat.R
```

Variants holding temperature or precipitation fixed:

```bash
Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_tempfixed_act.R 1
Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_tempfixed_nat.R 1

Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_precfixed_act.R 1
Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_precfixed_nat.R 1
```

**Outputs (per model):**
- `data/sensitivity-act/predictions_regional_{model}.csv`
- `data/sensitivity-act/summary_stats_regional_{model}.csv`
- Equivalent files under `data/sensitivity-nat/`

---

#### Step 4 — Coefficients (`04_coefficients/`)

Parallelised coefficient bootstrap (16 models × 1000 iterations, 160 SLURM tasks) and Bayesian spatial model comparison with INLA.

```bash
# Single coefficient bootstrap task:
Rscript 04_coefficients/04_coefficients_full.R 1

# Aggregate:
Rscript 04_coefficients/04_aggregate_coefficients_full.R

# INLA Bayesian comparison (no task ID required):
Rscript 04_coefficients/04_inla_compare_statistics.R

# INLA for Mexico:
Rscript 04_coefficients/04_inla_model_mexico.R

# Year-by-year fits:
Rscript 04_coefficients/04_fitted_year_individually.R 1
```

**Outputs:**
- `data/{model_name}_coef_state_blockboot1000.csv` — coefficient bootstrap per model
- `data/all_models_original_coefficients.csv` — combined
- `04_coefficients/results/{model_tag}/` — INLA posterior summaries and samples

---

#### Step 5 — Figures

Open JupyterLab and run any `Figure*.ipynb` notebook. All notebooks read from the outputs of steps 1–4 and write figures to `img/`.

```bash
# Local (requires Python + Jupyter installed):
jupyter lab

# Docker (recommended):
docker compose up jupyter
# Open http://localhost:8888
```

**Main figures:**

| Notebook | Description |
|---|---|
| `Figure1.ipynb` | Temporal/spatial distribution and immunity proxies |
| `Figure2.ipynb` | Temperature response coefficients and model comparison |
| `Figure3.ipynb` | Attribution: actual vs. counterfactual climate |
| `Figure4.ipynb` | City-level sensitivity and immunity effects |

**Supplementary figures:** `FigureS1`, `FigureS3–S6`, `FigureS10–S17`, `FigureS21–S23`.

---

### Reproducing all manuscript results

Full reproduction requires an HPC cluster with SLURM. Steps 2–4 each consist of array jobs with hundreds to thousands of tasks; see the [Running on a cluster](#running-on-a-cluster-slurm) section below. Steps 1 and 5 can be run locally or via Docker. A complete description of the statistical methodology is provided in the Methods section of the accompanying manuscript.

---

## Running on a cluster (SLURM)

All compute-intensive steps ship with pre-written `.sbatch` files targeting the JASMIN HPC cluster. Adapt queue names and paths for other clusters.

```bash
# Step 2: Bootstrap R²/AIC (60 tasks)
sbatch 02_statistics/02_sqcorr_full.sbatch

# Step 2: Out-of-sample validation (80–83 tasks)
sbatch 02_statistics/02_outsample_validation.sbatch

# Step 3: Sensitivity — actual climate (1700 tasks)
sbatch 03_sensitivity/03_submit_bootstrap_act.sbatch

# Step 3: Sensitivity — natural/counterfactual climate (1700 tasks)
sbatch 03_sensitivity/03_submit_bootstrap_nat.sbatch

# Step 3: Sensitivity — temperature fixed (1700 tasks each)
sbatch 03_sensitivity/03_submit_bootstrap_tempfixed_act.sbatch
sbatch 03_sensitivity/03_submit_bootstrap_tempfixed_nat.sbatch

# Step 3: Sensitivity — precipitation fixed (1700 tasks each)
sbatch 03_sensitivity/03_submit_bootstrap_precfixed_act.sbatch
sbatch 03_sensitivity/03_submit_bootstrap_precfixed_nat.sbatch

# Step 4: Coefficient bootstrap (160 tasks)
sbatch 04_coefficients/04_coefficients_full.sbatch

# After each step's array completes, run the corresponding aggregate script:
Rscript 02_statistics/02_aggreagte_full_sqcorr.R
Rscript 03_sensitivity/03_aggregate_results_act.R
Rscript 03_sensitivity/03_aggregate_results_nat.R
Rscript 04_coefficients/04_aggregate_coefficients_full.R
```

Each SLURM task reads `$SLURM_ARRAY_TASK_ID` automatically, so no additional arguments are needed when submitting via `sbatch`.

---

## Running with Docker

The Docker image bundles R 4.5.1, all R packages (including INLA pre-built for x86_64), Python 3.9, and JupyterLab.

### Build

```bash
docker compose build
```

> First build takes ~30–60 minutes due to INLA compilation. Subsequent builds use the layer cache.

### Verify the environment

```bash
docker compose run --rm jupyter python3 tests/test_python_packages.py
docker compose run --rm jupyter Rscript tests/test_r_packages.R
```

### Run R scripts

```bash
# Step 1 – fit models
docker compose run --rm jupyter Rscript 01_model_fit/fit_model.R

# Step 3 – single sensitivity task
docker compose run --rm jupyter Rscript 03_sensitivity/03_fit_sensitivity_parallel_full_act.R 1
```

### Run Jupyter (figures)

```bash
docker compose up jupyter
```

Then open `http://localhost:8888` in your browser.

---

## Data

Input data are not versioned (excluded via `.gitignore`). After downloading from Mendeley Data (see [Demo](#demo) above), the `data/` directory should contain:

**Climate and Epidemiological Inputs:**
- `model_input_brazil_immunity_city_with_priorinf.csv` — Climate variables, epidemiological, and socioeconomic covariates by city

**Bootstrap and counterfactual ensemble:**
- `bootstrap_state_samples.csv` — pre-sampled bootstrap state indices (1000 iterations)
- `predict-all/ext-ens000.csv` … `ext-ens524.csv` — actual climate ensemble (525 files)
- `predict-all/ext-nat-ens000.csv` … `ext-nat-ens524.csv` — natural/counterfactual ensemble

---

## License

This code is released under the [MIT License](LICENSE). See the `LICENSE` file for details.

---

## Citation

If you use this code or data, please cite the accompanying preprint:

> [Author(s)], "The Role of Climate Change in the Expansion of Dengue", *medRxiv*, 2025. [https://doi.org/10.1101/2025.10.06.25337235](https://doi.org/10.1101/2025.10.06.25337235)

Data citation:

> [Author(s)], "Dengue Attribution Dataset", Mendeley Data, 2025. [https://doi.org/10.17632/dyx8r83ktx.1](https://doi.org/10.17632/dyx8r83ktx.1)
