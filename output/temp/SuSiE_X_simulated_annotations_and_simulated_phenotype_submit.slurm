#!/bin/bash
set -euo pipefail
#SBATCH --time=1:00:00
#SBATCH --job-name="SuSiE X, simulated annotations, and simulated phenotype"
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000
#SBATCH --mail-user=mgc5166@psu.edu
#SBATCH --mail-type=ALL
# Paths below are relative to the directory where you run `sbatch` (vignettes/).
#SBATCH --output=../output/slurm_output/%x-%A-%a-%j.out
#SBATCH --error=../output/slurm_output/%x-%A-%a-%j.err

cd /storage/work/mgc5166/test_susine/vignettes/
PROJ_ROOT=$(realpath ..)
PRINTS_DIR="$PROJ_ROOT/output/slurm_prints"
RUN_DIR="$PROJ_ROOT/output/run_history/SuSiE_X_simulated_annotations_and_simulated_phenotype_$SLURM_ARRAY_JOB_ID"
mkdir -p "$PRINTS_DIR" "$RUN_DIR"
echo "started $(date) on $HOSTNAME" > "$RUN_DIR/meta.txt"
echo "task=$SLURM_ARRAY_TASK_ID job=$SLURM_ARRAY_JOB_ID" >> "$RUN_DIR/meta.txt"

module load r

Rscript simulation_p_15_zoom_in.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID
