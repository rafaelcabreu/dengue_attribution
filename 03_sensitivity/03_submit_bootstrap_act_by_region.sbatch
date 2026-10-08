#!/bin/bash
#SBATCH --job-name=dengue_bootstrap_by_region
#SBATCH --array=1-1900                  # 19 models x 100 bootstraps (region loop happens inside each task)
#                                       # tasks 1-100    -> model 1, bootstraps 1-100
#                                       # tasks 101-200  -> model 2, bootstraps 1-100
#                                       # ...
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1               # single-threaded per task; parallelism from the array
#SBATCH --mem=16G                       # adjust if ensemble files are large
#SBATCH --time=12:00:00                 # each task now fits ~5 region models instead of 1 — longer than the national script
#SBATCH --output=~/cpdn_nonnerc/aaim/dengue/sensitivity-act/logs/bootstrap_by_region_%A_%a.out
#SBATCH --error=~/cpdn_nonnerc/aaim/dengue/sensitivity-act/logs/bootstrap_by_region_%A_%a.err
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long

# ---------------------------------------------------------------
# Environment
# ---------------------------------------------------------------
mkdir -p ~/cpdn_nonnerc/aaim/dengue/sensitivity-act/logs
mkdir -p ~/cpdn_nonnerc/aaim/dengue/sensitivity-act/bootstrap_chunks_by_region

module load jasr                        # Jasmin R module; check available version with: module avail

# ---------------------------------------------------------------
# Run: pass the SLURM array task ID; the R script maps it to
# (model_idx, bootstrap_idx) internally, then loops over regions.
# ---------------------------------------------------------------
echo "Starting task: $SLURM_ARRAY_TASK_ID"
Rscript 03_fit_sensitivity_parallel_full_act_by_region.R $SLURM_ARRAY_TASK_ID
echo "Finished task: $SLURM_ARRAY_TASK_ID"
