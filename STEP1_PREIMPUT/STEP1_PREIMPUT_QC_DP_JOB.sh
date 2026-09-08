#!/usr/bin/env bash
#####################################################################
################ STEP1 PREIMPUT: QC + DP + PCA ######################
#####################################################################
# Coordinator script. The implementation is split into:
#   main_scripts/00_config_helpers.sh
#   main_scripts/01_step1A_QC.sh
#   main_scripts/02_step1B_DP.sh
#   main_scripts/03_step1C_PCA.sh
#####################################################################

set -euo pipefail # If something fails, the execution stops 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

source "main_scripts/00_config_helpers.sh"
source "main_scripts/01_step1A_QC.sh"
source "main_scripts/02_step1B_DP.sh"
source "main_scripts/03_step1C_PCA.sh"

#####################################################################
# BATCH DETECTION
#####################################################################
echo "Detecting batches inside: $INPUT_DIR"

ALL_BATCHES=()
for d in "${INPUT_DIR}"/*_geno_rawdata; do
    [[ -d "$d" ]] || continue
    batch=$(basename "$d" | sed 's/_geno_rawdata$//')
    ALL_BATCHES+=("$batch")
done

if [[ ${#ALL_BATCHES[@]} -eq 0 ]]; then
    echo "ERROR: No batches detected in ${INPUT_DIR}"
    echo "Expected: ${INPUT_DIR}/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_geno_rawdata.{ped/map or bed/bim/fam}"
    exit 1
fi

echo ""
echo "Detected batches:"
printf '  - %s\n' "${ALL_BATCHES[@]}"

#####################################################################
# NON-INTERACTIVE JOB PARAMETERS, set these variables in the sbatch launcher submit_STEP1_PREIMPUT.sbatch:
#--------------------------------------------------------------------
#   BATCH_MODE=ALL|SELECTED
#   BATCH_LIST="BATCH_NAME1 BATCH_NAME2"   # required when mode is SELECTED, ignored if ALL 

# Run DataProcessing (TOPMED cookbook)? 
#   RUN_STEP1B_DP=YES|NO

# Run Principal Component Analysis (with ancestry inference)? 
#   RUN_STEP1C_PCA=YES|NO
#####################################################################

BATCH_MODE="${BATCH_MODE:-SELECTED}"
BATCH_LIST="${BATCH_LIST:-}"
RUN_STEP1B_DP="${RUN_STEP1B_DP:-NO}"
RUN_STEP1C_PCA="${RUN_STEP1C_PCA:-NO}"

BATCH_MODE="${BATCH_MODE^^}"
RUN_STEP1B_DP="${RUN_STEP1B_DP^^}"
RUN_STEP1C_PCA="${RUN_STEP1C_PCA^^}"

if [[ "$RUN_STEP1B_DP" != "YES" && "$RUN_STEP1B_DP" != "NO" ]]; then
    echo "ERROR: RUN_STEP1B_DP must be YES or NO."
    exit 1
fi

if [[ "$RUN_STEP1C_PCA" != "YES" && "$RUN_STEP1C_PCA" != "NO" ]]; then
    echo "ERROR: RUN_STEP1C_PCA must be YES or NO."
    exit 1
fi

case "$BATCH_MODE" in
    ALL)
        SELECTED_BATCHES=("${ALL_BATCHES[@]}")
        ;;
    SELECTED)
        read -r -a USER_BATCHES <<< "$BATCH_LIST"
        if [[ ${#USER_BATCHES[@]} -eq 0 ]]; then
            echo "ERROR: BATCH_LIST is empty while BATCH_MODE=SELECTED."
            exit 1
        fi

        SELECTED_BATCHES=()
        for batch in "${USER_BATCHES[@]}"; do
            if [[ " ${ALL_BATCHES[*]} " =~ " ${batch} " ]]; then
                SELECTED_BATCHES+=("$batch")
            else
                echo "ERROR: Batch '$batch' not found."
                echo "Available batches are:"
                printf '  - %s\n' "${ALL_BATCHES[@]}"
                exit 1
            fi
        done
        ;;
    *)
        echo "ERROR: BATCH_MODE must be ALL or SELECTED."
        exit 1
        ;;
esac

echo ""
echo "Batches to process:"
printf '  -> %s\n' "${SELECTED_BATCHES[@]}"
echo ""

echo ""
echo "Selected downstream steps:"
echo "  STEP 1B_DP  -> $RUN_STEP1B_DP"
echo "  STEP 1C_PCA -> $RUN_STEP1C_PCA"
echo ""


#####################################################################
# MAIN LOOP
#####################################################################
for BATCH_NAME in "${SELECTED_BATCHES[@]}"; do
    echo "================================================="
    echo "Starting combined preimput pipeline for batch: $BATCH_NAME"
    echo "================================================="

    INPUT_BATCH_DIR="${INPUT_DIR}/${BATCH_NAME}_geno_rawdata"
    BATCH_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${BATCH_NAME}_output_preimput"
    QC_OUTPUT_DIR="${BATCH_OUTPUT_DIR}/1A_output_QC"
    DP_OUTPUT_DIR="${BATCH_OUTPUT_DIR}/1B_output_DP"
    PCA_OUTPUT_DIR="${BATCH_OUTPUT_DIR}/1C_output_PCA"
    BATCH_INFO_DIR="${INFO_BASE_DIR}/${BATCH_NAME}_info_preimput"
    INFO_QC_DIR="${BATCH_INFO_DIR}/_info_STEP1_QC"
    INFO_QC_LOG_DIR="${INFO_QC_DIR}/logs"
    INFO_DP_DIR="${BATCH_INFO_DIR}/_info_STEP2_DP"
    INFO_PCA_DIR="${BATCH_INFO_DIR}/_info_STEP3_PCA"
    DP_WORK_DIR="${BATCH_OUTPUT_DIR}/_work_STEP2_DP"
    PCA_WORK_DIR="${BATCH_OUTPUT_DIR}/_work_STEP3_PCA"

    mkdir -p "$BATCH_OUTPUT_DIR" "$QC_OUTPUT_DIR" "$DP_OUTPUT_DIR" "$PCA_OUTPUT_DIR" "$BATCH_INFO_DIR" "$INFO_QC_DIR" "$INFO_QC_LOG_DIR" "$INFO_DP_DIR" "$INFO_PCA_DIR" "$DP_WORK_DIR" "$PCA_WORK_DIR"

    run_step1a_qc
    run_step1b_dp
    run_step1c_pca
done


#####################################################################
# INFORMATION ABOUT EXECUTION FOR THE USER 
#####################################################################
echo "### PIPELINE SUMMARY ###"
echo "### STEP 1A_QC completed for selected batch(es)."

#--------------------------------------------------------------------
# Message for the user if 1B == YES 
#--------------------------------------------------------------------
if [[ "$RUN_STEP1B_DP" == "YES" ]]; then
    echo "### STEP 1B_DP completed: VCF files are ready for Michigan Imputation."
    echo "### MICHIGAN IMPUTATION ###"
    echo "Load QCtargetdata-updated-chr*.vcf.gz files to Michigan Imputation Server."
    echo "# Reference panel: HRC r1.1 2016 (GRCh37/hg19)"
    echo "# Array Build: GRCh37/hg19 (depending on the target base)"
    echo "# rsq filter: 0.3"
    echo "# Phasing: Eagle v2.4 (phased output)"
    echo "# Population: EUR"
    echo "# Mode: Quality Control & Imputation"
else
    echo "### STEP 1B_DP skipped: no imputation VCF files were created."
fi
#--------------------------------------------------------------------
# Message for the user if 1C == YES 
#--------------------------------------------------------------------
if [[ "$RUN_STEP1C_PCA" == "YES" ]]; then
    echo "### STEP 1C_PCA completed."
else
    echo "### STEP 1C_PCA skipped."
fi

if [[ "$RUN_STEP1B_DP" == "NO" && "$RUN_STEP1C_PCA" == "NO" ]]; then
    echo "### QC-only run finished."
else
    echo "### Selected preimputation steps finished."
fi
