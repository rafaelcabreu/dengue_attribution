#!/bin/bash
#SBATCH --job-name=dengue_boot_year
#SBATCH --array=1-25                        # TODO: replace ?? with the number of years in your data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/boot_year_%a.out
#SBATCH --error=logs/boot_year_%a.err
#SBATCH --partition=short-serial        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --partition=standard        # Jasmin partition — change to lotus/high-mem if needed
#SBATCH --account=cpdn_nonnerc
#SBATCH --qos=long

# ---------------------------------------------------------------
# Create log directory if it doesn't exist
# ---------------------------------------------------------------
mkdir -p logs

# ---------------------------------------------------------------
# Load R module (adjust version/module name to your cluster)
# ---------------------------------------------------------------
module load jasr

# ---------------------------------------------------------------
# Run the bootstrap script for this array task
# SLURM_ARRAY_TASK_ID is passed as the task_id argument,
# which the R script maps to a calendar year automatically.
# ---------------------------------------------------------------
echo "Starting task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date)"

Rscript fitted_year_individually.R ${SLURM_ARRAY_TASK_ID}

EXIT_CODE=$?
echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date) with exit code ${EXIT_CODE}"
exit ${EXIT_CODE}