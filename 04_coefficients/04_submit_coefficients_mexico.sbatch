#!/bin/bash
#SBATCH --job-name=dengue_mexico_coef_boot
#SBATCH --array=1-20                    # n_models(2) * n_chunks(10) = 20 tasks
#                                       # tasks 1-10  -> model 1 (Brazil + Mexico), chunks 1-10
#                                       # tasks 11-20 -> model 2 (Mexico only),     chunks 1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1               # single-threaded; parallelism from the array
#SBATCH --mem=16G
#SBATCH --time=06:00:00                 # 100 bootstrap fits per task; adjust to your runtime
#SBATCH --output=./statistics/logs/mexico_coef_boot_%A_%a.out
#SBATCH --error=./statistics/logs/mexico_coef_boot_%A_%a.err
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long

# ---------------------------------------------------------------
# Environment
# ---------------------------------------------------------------
mkdir -p ./statistics/logs
mkdir -p ./data/coefficients/bootstrap_chunks_mexico

module load jasr                        # Jasmin R module; check available version with: module avail
# If you use a personal R library, uncomment and set:
# export R_LIBS_USER=$HOME/R/library

# ---------------------------------------------------------------
# Run: pass the SLURM array task ID; the R script maps it to
# (model_idx, chunk_idx) internally.
# ---------------------------------------------------------------
echo "Starting task: $SLURM_ARRAY_TASK_ID"
Rscript 04_coefficients/04_coefficients_full_mexico.R $SLURM_ARRAY_TASK_ID
echo "Finished task: $SLURM_ARRAY_TASK_ID"
