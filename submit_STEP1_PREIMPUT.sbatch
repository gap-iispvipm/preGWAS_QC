#!/usr/bin/env bash
#SBATCH --job-name=STEP1_PREIMPUT
#SBATCH --output=STEP1_PREIMPUT_INSCRIPTION.out
#SBATCH --error=STEP1_PREIMPUT_INSCRIPTION.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G

set -euo pipefail

# Submit this job from the STEP1_PREIMPUT directory.
PIPELINE_DIR="${STEP1_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$PIPELINE_DIR"

#####################################################################
# JOB PARAMETERS
#####################################################################

# ALL: process every detected batch.
# SELECTED: process only the names written in BATCH_LIST.
export BATCH_MODE="SELECTED"

# Separate multiple batch names with spaces.
# Ignored when BATCH_MODE="ALL".
export BATCH_LIST="INSCRIPTION"

# STEP 1A QC always runs.
export RUN_STEP1B_DP="YES"
export RUN_STEP1C_PCA="YES"

echo "### SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "### Node: $(hostname)"
echo "### Working directory: $PIPELINE_DIR"
echo "### BATCH_MODE: $BATCH_MODE"
echo "### BATCH_LIST: $BATCH_LIST"
echo "### RUN_STEP1B_DP: $RUN_STEP1B_DP"
echo "### RUN_STEP1C_PCA: $RUN_STEP1C_PCA"

bash "${PIPELINE_DIR}/STEP1_PREIMPUT_QC_DP_JOB.sh"