#!/bin/bash
#SBATCH --job-name=dengue_rsq_aic_boot
#SBATCH --array=1-240                   # model 23 (newly added) x n_chunks(10) = tasks 221-230
#                                       # tasks 1-10   -> model 1, chunks 1-10
#                                       # tasks 11-20  -> model 2, chunks 1-10
#                                       # ...
#                                       # tasks 221-230 -> model 23, chunks 1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1               # single-threaded per task; parallelism from the array
#SBATCH --mem=16G                       # adjust if ensemble files are large
#SBATCH --time=10:00:00                 # adjust to your expected runtime per task
#SBATCH --output=./statistics/logs/bootstrap_%A_%a.out
#SBATCH --error=./statistics/logs/bootstrap_%A_%a.err
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long


# ---------------------------------------------------------------
# Environment
# ---------------------------------------------------------------
mkdir -p ./statistics/logs
mkdir -p ./statistics/bootstrap_chunks_rsq

module load jasr                        # Jasmin R module; check available version with: module avail
# If you use a personal R library, uncomment and set:
# export R_LIBS_USER=$HOME/R/library

# ---------------------------------------------------------------
# Run: pass the SLURM array task ID; the R script maps it to
# (model_idx, chunk_idx) internally.
# ---------------------------------------------------------------
echo "Starting task: $SLURM_ARRAY_TASK_ID"
Rscript 02_bootstrap_parallel_full_sqcorr.R $SLURM_ARRAY_TASK_ID
echo "Finished task: $SLURM_ARRAY_TASK_ID"