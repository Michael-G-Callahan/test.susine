#!/bin/bash
#SBATCH --job-name=susine_sim_grid_demo
#SBATCH --array=1-36
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=mgc5166@psu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output/slurm_prints/%x-%j-%a.out"
#SBATCH --error="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output/slurm_prints/%x-%j-%a.err"

set -euo pipefail

echo "[$(date -Is)] Starting task ${SLURM_ARRAY_TASK_ID} for job ${SLURM_JOB_NAME}"
JOB_ROOT="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output"
CONFIG_PATH="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output/run_history/susine_sim_grid_demo/job_config.json"
RUN_TASK_SCRIPT="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/inst/scripts/run_task.R"

module load r

export R_LIBS_USER="/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3"

Rscript "$RUN_TASK_SCRIPT" \
  --job-name "$SLURM_JOB_NAME" \
  --task-id "$SLURM_ARRAY_TASK_ID" \
  --job-root "$JOB_ROOT" \
  --config-path "$CONFIG_PATH"

# Reorganize stdout/stderr into nested folders
RAW_OUT_FILE="${SLURM_JOB_NAME}-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}.out"
RAW_ERR_FILE="${SLURM_JOB_NAME}-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}.err"
RAW_OUT_PATH="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output/slurm_prints/${RAW_OUT_FILE}"
RAW_ERR_PATH="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output/slurm_prints/${RAW_ERR_FILE}"
FINAL_DIR="/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output/slurm_prints/${SLURM_JOB_NAME}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}"
mkdir -p "$FINAL_DIR"
if [ -f "$RAW_OUT_PATH" ]; then mv "$RAW_OUT_PATH" "$FINAL_DIR/stdout.out"; fi
if [ -f "$RAW_ERR_PATH" ]; then mv "$RAW_ERR_PATH" "$FINAL_DIR/stderr.err"; fi
echo "[$(date -Is)] Completed task ${SLURM_ARRAY_TASK_ID}"
