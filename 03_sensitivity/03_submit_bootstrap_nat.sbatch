#!/bin/bash
#SBATCH --job-name=dengue_bootstrap
#SBATCH --array=1-1700                  # n_bootstrap * n_models = 100 * 17 = 1600 tasks
#                                       # tasks 1-100   -> model 1, bootstraps 1-100
#                                       # tasks 101-200 -> model 2, bootstraps 1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1               # single-threaded per task; parallelism from the array
#SBATCH --mem=16G                       # adjust if ensemble files are large
#SBATCH --time=06:00:00                 # adjust to your expected runtime per task
#SBATCH --output=~/cpdn_nonnerc/aaim/dengue/sensitivity-nat/logs/bootstrap_%A_%a.out
#SBATCH --error=~/cpdn_nonnerc/aaim/dengue/sensitivity-nat/logs/bootstrap_%A_%a.err
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long

# ---------------------------------------------------------------
# Environment
# ---------------------------------------------------------------
mkdir -p ~/cpdn_nonnerc/aaim/dengue/sensitivity-nat/logs
mkdir -p ~/cpdn_nonnerc/aaim/dengue/sensitivity-nat/bootstrap_chunks

module load jasr                        # Jasmin R module; check available version with: module avail

# ---------------------------------------------------------------
# Run: pass the SLURM array task ID; the R script maps it to
# (model_idx, bootstrap_idx) internally.
# ---------------------------------------------------------------
echo "Starting task: $SLURM_ARRAY_TASK_ID"
Rscript 03_fit_sensitivity_parallel_full_nat.R $SLURM_ARRAY_TASK_ID
echo "Finished task: $SLURM_ARRAY_TASK_ID"
