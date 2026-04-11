#!/bin/bash
#SBATCH --job-name=dengue_coef_boot
#SBATCH --array=1-180                   # n_models(16) * n_chunks(10) = 160 tasks
#                                       # tasks 1-10    -> model 1,  chunks 1-10
#                                       # tasks 11-20   -> model 2,  chunks 1-10
#                                       # ...
#                                       # tasks 151-160 -> model 16, chunks 1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1               # single-threaded; parallelism from the array
#SBATCH --mem=16G                        # coefficient bootstrap is lighter than ensemble prediction
#SBATCH --time=06:00:00                 # 100 bootstrap fits per task; adjust to your runtime
#SBATCH --output=/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/coefficients/logs/coef_boot_%A_%a.out
#SBATCH --error=/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/coefficients/logs/coef_boot_%A_%a.err
#SBATCH --partition=short-serial        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long

# ---------------------------------------------------------------
# Environment
# ---------------------------------------------------------------
mkdir -p /gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/coefficients/logs
mkdir -p /gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/coefficients/bootstrap_chunks

module load jasr                        # Jasmin R module; check available version with: module avail
# If you use a personal R library, uncomment and set:
# export R_LIBS_USER=$HOME/R/library

# ---------------------------------------------------------------
# Run: pass the SLURM array task ID; the R script maps it to
# (model_idx, chunk_idx) internally.
# ---------------------------------------------------------------
echo "Starting task: $SLURM_ARRAY_TASK_ID"
Rscript coefficients_full.R $SLURM_ARRAY_TASK_ID
echo "Finished task: $SLURM_ARRAY_TASK_ID"