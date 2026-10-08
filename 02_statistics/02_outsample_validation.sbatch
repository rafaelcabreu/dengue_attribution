#!/bin/bash
#SBATCH --job-name=dengue_outsample_boot
#SBATCH --account=cpdn_nonnerc
#SBATCH --time=480:00
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --qos=long
#SBATCH --array=1-91   # model 22 (newly added), 4 chunks of 25 iterations (0-indexed)

# ── Directory configuration ───────────────────────────────────────────────────
LOG_DIR="/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/outsample-statistics/logs"
RESULTS_DIR="/gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/outsample-statistics/bootstrap_chunks"

# ── SLURM log paths ───────────────────────────────────────────────────────────
#SBATCH -o /gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/outsample-statistics/logs/%A_%a.out
#SBATCH -e /gws/ssde/j25a/cpdn_nonnerc/aaim/dengue/outsample-statistics/logs/%A_%a.err

# ── Array configuration ───────────────────────────────────────────────────────
# 19 models, each split into 4 chunks of 25 iterations
# Model indices are 1-based in R; array tasks are 0-based in SLURM
#
# Task → model index mapping (1-based):
#   tasks  0– 3  → model  1: Year + Month
#   tasks  4– 7  → model  2: Year|region + Month|region
#   tasks  8–11  → model  3: Climate(lag 1-5)
#   tasks 12–15  → model  4: Climate(lag 1-5) + Year|region
#   tasks 16–19  → model  5: Climate(lag 1-5) + Month|region
#   tasks 20–23  → model  6: Year + Month + P(Climate)
#   tasks 24–27  → model  7: Climate(lag 1) + Year|region + Month|region
#   tasks 28–31  → model  8: Climate(lag 1-3) + Year|region + Month|region
#   tasks 32–35  → model  9: Climate(lag 1-5) + Year|region + Month|region
#   tasks 36–39  → model 10: Climate(lag 1-5) + Year|region + Month|region + Socio
#   tasks 40–43  → model 11: Climate(lag 1-5) + Year|region + Month|region + PriorCases
#   tasks 44–47  → model 12: Climate(lag 1-5) + Year|region + Month|region + SeroRepla
#   tasks 48–51  → model 13: Climate(lag 1-5) + Year|region + Month|region + Immunity
#   tasks 52–55  → model 14: Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio + Immunity
#   tasks 56–59  → model 15: Climate(lag 1-5) + Socio
#   tasks 60–63  → model 16: Climate(lag 1-5) + PriorCases
#   tasks 64–67  → model 17: Climate(lag 1-5) + SeroRepla
#   tasks 68–71  → model 18: Climate(lag 1-5) + Immunity
#   tasks 72–75  → model 19: Climate(lag 1-5) + PriorCases + SeroRepla + Socio + Immunity
#   (models 20-21 added later, undocumented here — see R script for names)
#   tasks 84–87  → model 22: Climate(lag 1-5) + Year|region + Month|region + PriorCases + SeroRepla + Socio

CHUNKS_PER_MODEL=4
ITERS_PER_CHUNK=25

# Derive 1-based model index and chunk index from 0-based array task ID
MODEL_IDX=$(( SLURM_ARRAY_TASK_ID / CHUNKS_PER_MODEL + 1 ))
CHUNK_IDX=$(( SLURM_ARRAY_TASK_ID % CHUNKS_PER_MODEL ))

START_ITER=$(( CHUNK_IDX * ITERS_PER_CHUNK + 1 ))
END_ITER=$(( START_ITER + ITERS_PER_CHUNK - 1 ))

# Create directories if they don't exist
mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

echo "Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "Model index: ${MODEL_IDX}"
echo "Iterations:  ${START_ITER} to ${END_ITER}"
echo "Log dir:     ${LOG_DIR}"
echo "Results dir: ${RESULTS_DIR}"
echo "Started:     $(date)"

module load jasr

Rscript 02_fit_outsample_validation.R "${MODEL_IDX}" "${START_ITER}" "${END_ITER}" "${RESULTS_DIR}"

echo "Finished: $(date)"
