#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name="SuSiNE Simulations - informative hierarchical annotations"
#SBATCH --mem=1GB
#SBATCH --array=1-150
#SBATCH --output=/storage/work/mgc5166/simulation_batch_output/JOB_ID_%A/test-job-%j-task-%a.out

WORKDIR=/storage/work/mgc5166/simulation_batch_output/JOB_ID_$SLURM_ARRAY_JOB_ID
mkdir -p "$WORKDIR" || exit -1 

cd /storage/work/mgc5166/susine/vignettes/

# This may vary per HPC system. At USC's hpc system
# we use: source /usr/usc/R/default/setup.sh
module load r

Rscript simulation_p_10_hierarchical_sim_annotations.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID