#!/bin/bash
#SBATCH --job-name=dengue_mexico_nat
#SBATCH --array=1-2                     # 2 models (Brazil+Mexico training, Mexico-only training)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00                 # single fit + 525 predict() calls per task, no bootstrap
#SBATCH --output=./statistics/logs/mexico_nat_%A_%a.out
#SBATCH --error=./statistics/logs/mexico_nat_%A_%a.err
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long

# ---------------------------------------------------------------
# Environment
# ---------------------------------------------------------------
mkdir -p ./statistics/logs
mkdir -p ./data/sensitivity-mexico-nat

module load jasr                        # Jasmin R module; check available version with: module avail

# ---------------------------------------------------------------
# Run: pass the SLURM array task ID directly as model_idx (1 or 2).
# ---------------------------------------------------------------
echo "Starting task: $SLURM_ARRAY_TASK_ID"
Rscript 03_sensitivity/03_fit_sensitivity_parallel_mexico_nat.R $SLURM_ARRAY_TASK_ID
echo "Finished task: $SLURM_ARRAY_TASK_ID"
