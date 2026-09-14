#!/usr/bin/env bash
#SBATCH --job-name=STEP2_gepi
#SBATCH --output=STEP2_POSTIMPUT_gepi.out
#SBATCH --error=STEP2_POSTIMPUT_gepi.err
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=mr-06

set -euo pipefail

PIPELINE_DIR="${STEP2_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$PIPELINE_DIR"

export BATCH_MODE="SELECTED"
export BATCH_LIST="GEPI-BIOPSI"

echo "### SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "### Node: $(hostname)"
echo "### Working directory: $PIPELINE_DIR"
echo "### BATCH_MODE: $BATCH_MODE"
echo "### BATCH_LIST: $BATCH_LIST"

bash "${PIPELINE_DIR}/STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh"