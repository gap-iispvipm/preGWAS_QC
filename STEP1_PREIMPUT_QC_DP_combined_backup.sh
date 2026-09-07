#!/usr/bin/env bash
#####################################################################
#################### STEP1 PREIMPUT: QC + DP ########################
#####################################################################
# Combined pre-imputation pipeline:
#   1A_QC: original STEP1_QC_TARGET_DATA_PIPELINE order
#   1B_DP: original STEP2_PREIMPUTATION_DATA_PROCESSING order
#
# Folder layout:
#   STEP1_PREIMPUT/
#   |-- A_input_geno_rawdata/
#   |-- B_output_preimput/
#   |-- info_preimput/
#   |-- scripts_1A_QC/
#   |-- scripts_1B_DP/
#   `-- STEP1_PREIMPUT_QC_DP.sh
#
# The DP step no longer reads A_input_preimput_processing. It reads the
# PLINK dataset created by 1A_QC:
#   B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC
#####################################################################

set -euo pipefail

SCRIPT_VERSION="2026-06-14_QC_DP_PCA_v2_relatedness_summary"
ROOT_DIR="$(pwd)"
INPUT_DIR="A_input_geno_rawdata"
PCA_REF_DIR="1000G_ref_hg19"
OUTPUT_BASE_DIR="B_output_preimput"
INFO_BASE_DIR="info_preimput"
QC_SCRIPT_DIR="scripts/scripts_1A_QC"
DP_SCRIPT_DIR="scripts/scripts_1B_DP"
PCA_SCRIPT_DIR="scripts/scripts_1C_PCA"

echo "### Running STEP1_PREIMPUT_QC_DP.sh version: ${SCRIPT_VERSION}"

mkdir -p "$INPUT_DIR" "$PCA_REF_DIR" "$OUTPUT_BASE_DIR" "$INFO_BASE_DIR" "$QC_SCRIPT_DIR" "$DP_SCRIPT_DIR" "$PCA_SCRIPT_DIR"

#====================================================================
# MODULES NEEDED
#====================================================================
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
module load BioPerl/1.7.8-GCCcore-12.3.0
module load BCFtools/1.21-GCC-12.3.0

#====================================================================
# USER PARAMETERS
#====================================================================
runplink="plink"

# STEP 1A_QC parameters
step0fam="NO"
step1="YES"
step2="YES"
step3="YES"
step4="YES"
step5="YES"
step6="YES"

mymissigness1=0.2
mymissigness2=0.02
mymaf=0.01
myhwe1=1e-6
myhwe2=1e-10
mywindowsize=50
myshiftwindow=5
mypairwiser2=0.2
myrelatedness=0.125
mysplitx_build="b37"

# STEP 1B_DP parameters
s0liftoverTarget="NO"
s0liftoverReference="NO"
# runliftover="/path/to/liftOver"
# liftoverChain="/path/to/hg38ToHg19.over.chain"

#####################################################################
# SHARED HELPERS
#####################################################################
count_samples() {
    local base="$1"
    [[ -f "${base}.fam" ]] && wc -l < "${base}.fam" | awk '{print $1}' || echo "NA"
}

count_snps() {
    local base="$1"
    [[ -f "${base}.bim" ]] && wc -l < "${base}.bim" | awk '{print $1}' || echo "NA"
}

write_ids() {
    local base="$1"
    local out="$2"
    if [[ -f "${base}.fam" ]]; then
        awk '{print $1, $2}' "${base}.fam" | sort -u > "$out"
    else
        : > "$out"
    fi
}

write_snps() {
    local base="$1"
    local out="$2"
    if [[ -f "${base}.bim" ]]; then
        awk '{print $2}' "${base}.bim" | sort -u > "$out"
    else
        : > "$out"
    fi
}

diff_ids_removed() {
    local before="$1"
    local after="$2"
    local out="$3"
    comm -23 "$before" "$after" > "$out" || true
}

diff_snps_removed() {
    local before="$1"
    local after="$2"
    local out="$3"
    comm -23 "$before" "$after" > "$out" || true
}

count_related_flagged_individuals() {
    local related_file="$1"
    if [[ -f "$related_file" ]]; then
        awk 'NR>1 {print $1, $2}' "$related_file" | sort -u | wc -l | awk '{print $1}'
    else
        echo 0
    fi
}

count_related_flagged_pairs() {
    local related_file="$1"
    if [[ -f "$related_file" ]]; then
        awk 'NR>1 {n++} END {print int(n/2)}' "$related_file"
    else
        echo 0
    fi
}

cleanup_qc_tracking_files() {
    local info_dir="$1"
    find "$info_dir" -maxdepth 1 -type f \( \
        -name "*_before_ids.txt" -o \
        -name "*_after_ids.txt" -o \
        -name "*_before_snps.txt" -o \
        -name "*_after_snps.txt" -o \
        -name "initial_ids.txt" -o \
        -name "initial_snps.txt" -o \
        -name "final_ids.txt" -o \
        -name "final_snps.txt" \
    \) -delete
}

safe_pct_removed() {
    local start="$1"
    local end="$2"
    if [[ "$start" == "NA" || "$end" == "NA" || -z "$start" || -z "$end" || "$start" -eq 0 ]]; then
        echo "NA"
    else
        awk -v a="$start" -v b="$end" 'BEGIN{printf "%.2f", ((a-b)/a)*100}'
    fi
}

safe_removed() {
    local start="$1"
    local end="$2"
    if [[ "$start" == "NA" || "$end" == "NA" || -z "$start" || -z "$end" ]]; then
        echo "NA"
    else
        echo $((start - end))
    fi
}

append_state() {
    local table="$1"
    local step_name="$2"
    local base="$3"
    local reason_file="${4:-NA}"
    local n_ind n_snp
    n_ind=$(count_samples "$base")
    n_snp=$(count_snps "$base")
    printf "%s\t%s\t%s\t%s\n" "$step_name" "$n_ind" "$n_snp" "$reason_file" >> "$table"
}

make_transition_summary() {
    local states_file="$1"
    local out_file="$2"
    local title="$3"

    {
        echo "==============================================================="
        echo "$title"
        echo "==============================================================="
        printf "%-28s %-12s %-12s %-16s %-12s %-12s %-18s %-18s %s\n" \
            "STEP" "N_SAMPLES" "N_SNPS" "SAMPLES_REMOVED" "SNPS_REMOVED" \
            "PCT_IND" "PCT_SNP" "PCT_IND_START" "LIST_FILE"
        echo "--------------------------------------------------------------------------------------------------------------------------------"

        local start_ind="NA" start_snp="NA" prev_ind="NA" prev_snp="NA" first="YES"
        while IFS=$'\t' read -r step_name n_ind n_snp reason_file; do
            if [[ "$first" == "YES" ]]; then
                start_ind="$n_ind"
                start_snp="$n_snp"
                prev_ind="$n_ind"
                prev_snp="$n_snp"
                first="NO"
                printf "%-28s %-12s %-12s %-16s %-12s %-12s %-18s %-18s %s\n" \
                    "$step_name" "$n_ind" "$n_snp" "0" "0" "0" "0" "0" "$reason_file"
                continue
            fi

            local removed_ind removed_snp pct_ind pct_snp pct_ind_start
            removed_ind=$(safe_removed "$prev_ind" "$n_ind")
            removed_snp=$(safe_removed "$prev_snp" "$n_snp")
            pct_ind=$(safe_pct_removed "$prev_ind" "$n_ind")
            pct_snp=$(safe_pct_removed "$prev_snp" "$n_snp")
            pct_ind_start=$(safe_pct_removed "$start_ind" "$n_ind")

            printf "%-28s %-12s %-12s %-16s %-12s %-12s %-18s %-18s %s\n" \
                "$step_name" "$n_ind" "$n_snp" "$removed_ind" "$removed_snp" \
                "$pct_ind" "$pct_snp" "$pct_ind_start" "$reason_file"

            prev_ind="$n_ind"
            prev_snp="$n_snp"
        done < "$states_file"
        echo "==============================================================="
    } > "$out_file"
}

move_output_aux_to_info() {
    local output_dir="$1"
    local workdata="$2"
    local info_dir="$3"
    local prefix="$4"
    find "$output_dir" -name "${workdata}.*" \
        ! -name "*.bed" \
        ! -name "*.bim" \
        ! -name "*.fam" \
        -exec bash -c '
            f="$1"
            dest="$2"
            pfx="$3"
            base=$(basename "$f")
            mv "$f" "${dest}/${pfx}_${base}"
        ' _ {} "$info_dir" "$prefix" \;
}

require_r_script() {
    local script_path="$1"
    if [[ ! -f "$script_path" ]]; then
        echo "ERROR: Missing required R script: $script_path"
        exit 1
    fi
}

ask_yes_no() {
    local prompt="$1"
    local answer
    while true; do
        read -r -p "$prompt [YES/NO]: " answer
        case "${answer^^}" in
            YES|Y) echo "YES"; return ;;
            NO|N) echo "NO"; return ;;
            *) echo "Please answer YES or NO." ;;
        esac
    done
}

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
echo ""
echo "Select execution mode:"
echo "  1) Run ALL batches"
echo "  2) Run specific batch(es)"
read -r -p "Choose option [1/2]: " choice

if [[ "$choice" == "1" ]]; then
    SELECTED_BATCHES=("${ALL_BATCHES[@]}")
elif [[ "$choice" == "2" ]]; then
    echo "Enter batch names separated by spaces:"
    read -r -a USER_BATCHES
    if [[ ${#USER_BATCHES[@]} -eq 0 ]]; then
        echo "ERROR: No batch names provided."
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
else
    echo "Invalid option. Exiting."
    exit 1
fi

echo ""
echo "Batches to process:"
printf '  -> %s\n' "${SELECTED_BATCHES[@]}"
echo ""

RUN_STEP1B_DP=$(ask_yes_no "Do you want to run STEP 1B data processing for imputation?")
RUN_STEP1C_PCA=$(ask_yes_no "Do you want to run STEP 1C PCA?")

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

    #################################################################
    # STEP 1A_QC
    #################################################################
    QC_LOG="${INFO_QC_LOG_DIR}/${BATCH_NAME}_STEP1_QC.log"
    exec 3>&1 4>&2
    exec > >(tee -a "$QC_LOG") 2>&1

    echo "### STEP 1A_QC started for ${BATCH_NAME}"
    echo "Input directory : $INPUT_BATCH_DIR"
    echo "Output directory: $BATCH_OUTPUT_DIR"
    echo "Info directory  : $INFO_QC_DIR"
    echo "Log directory   : $INFO_QC_LOG_DIR"

    printf "### Script: STEP1_PREIMPUT_QC_DP
### Batch: %s
### Parameters:
step0fam=%s
step1_missingness=%s
mymissigness1=%s
mymissigness2=%s
step2_sex_discrepancy=%s
step3_maf=%s
mymaf=%s
step4_hwe=%s
myhwe1=%s
myhwe2=%s
step5_heterozygosity=%s
mywindowsize=%s
myshiftwindow=%s
mypairwiser2=%s
step6_relatedness=%s
myrelatedness=%s
mysplitx_build=%s
" "$BATCH_NAME" "$step0fam" "$step1" "$mymissigness1" "$mymissigness2" \
  "$step2" "$step3" "$mymaf" "$step4" "$myhwe1" "$myhwe2" "$step5" \
  "$mywindowsize" "$myshiftwindow" "$mypairwiser2" "$step6" "$myrelatedness" "$mysplitx_build" \
  > "${INFO_QC_DIR}/${BATCH_NAME}_STEP1_QC_parameters.txt"

    DATA="${INPUT_BATCH_DIR}/${BATCH_NAME}_geno_rawdata"
    WORKDATA="${BATCH_NAME}_working"
    WORKBASE="${BATCH_OUTPUT_DIR}/${WORKDATA}"
    QC_STATES="${INFO_QC_DIR}/${BATCH_NAME}_STEP1_QC_states.tsv"
    : > "$QC_STATES"

    ORIG_N_X="NA"
    ORIG_N_Y="NA"
    ORIG_N_XY_PAR="NA"
    ORIG_N_MT="NA"

    if [[ -f "${DATA}.bim" ]]; then
        ORIG_CHR_FILE="${DATA}.bim"
    elif [[ -f "${DATA}.map" ]]; then
        ORIG_CHR_FILE="${DATA}.map"
    else
        ORIG_CHR_FILE="NA"
    fi

    if [[ "$ORIG_CHR_FILE" != "NA" ]]; then
        ORIG_N_X=$(awk '$1 == "23" || $1 == "X" {n++} END {print n+0}' "$ORIG_CHR_FILE")
        ORIG_N_Y=$(awk '$1 == "24" || $1 == "Y" {n++} END {print n+0}' "$ORIG_CHR_FILE")
        ORIG_N_XY_PAR=$(awk '$1 == "25" || $1 == "XY" {n++} END {print n+0}' "$ORIG_CHR_FILE")
        ORIG_N_MT=$(awk '$1 == "26" || $1 == "MT" || $1 == "M" {n++} END {print n+0}' "$ORIG_CHR_FILE")
    fi

    SPLITX_INITIAL_STATUS="NA"
    SPLITX_INITIAL_BUILD="NA"
    SPLITX_FINAL_STATUS="NA"
    SPLITX_FINAL_BUILD="NA"

    if [[ -f "${DATA}.bed" && -f "${DATA}.bim" && -f "${DATA}.fam" ]]; then
        echo "Detected binary PLINK files (.bed/.bim/.fam) -> copying to working dataset"
        "$runplink" --bfile "$DATA" --make-bed --out "$WORKBASE"
        SKIP_CONVERSION="YES"
    elif [[ -f "${DATA}.ped" && -f "${DATA}.map" ]]; then
        echo "Detected PED/MAP files -> conversion needed"
        "$runplink" --file "$DATA" --make-bed --out "$WORKBASE"
        SKIP_CONVERSION="NO"
    else
        echo "ERROR: No valid input files found for ${BATCH_NAME}"
        echo "Expected ${DATA}.ped/.map or ${DATA}.bed/.bim/.fam"
        exit 1
    fi

    append_state "$QC_STATES" "0A_input_converted" "$WORKBASE" "NA"

    echo "### Starting Step 0: raw files script adaptation ###"

    write_snps "$WORKBASE" "${INFO_QC_DIR}/0C_before_snps.txt"
    "$runplink" --bfile "$WORKBASE" --snps-only 'just-acgt' --make-bed --out "$WORKBASE"
    move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "0C"
    write_snps "$WORKBASE" "${INFO_QC_DIR}/0C_after_snps.txt"
    diff_snps_removed "${INFO_QC_DIR}/0C_before_snps.txt" "${INFO_QC_DIR}/0C_after_snps.txt" "${INFO_QC_DIR}/0C_removed_non_acgt_multiallelic_snps.txt"
    append_state "$QC_STATES" "0C_snps_only_acgt" "$WORKBASE" "0C_removed_non_acgt_multiallelic_snps.txt"

    if [[ "$step0fam" == "YES" ]]; then
        echo "# Running Step 0 phenotypes"
        require_r_script "${QC_SCRIPT_DIR}/QC0_Target_data_phenotype.R"
        Rscript --no-save "${QC_SCRIPT_DIR}/QC0_Target_data_phenotype.R"

        "$runplink" \
            --bfile "$WORKBASE" \
            --pheno 0A_phenotypes.txt \
            --allow-no-sex \
            --make-bed \
            --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "0A"

        awk '$3==-9 {print $1,$2}' 0A_phenotypes.txt > "${INFO_QC_DIR}/0B_removed_missing_phenotype_iid.txt"
        write_ids "$WORKBASE" "${INFO_QC_DIR}/0B_before_ids.txt"
        "$runplink" \
            --bfile "$WORKBASE" \
            --remove "${INFO_QC_DIR}/0B_removed_missing_phenotype_iid.txt" \
            --allow-no-sex \
            --make-bed \
            --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "0B"
        write_ids "$WORKBASE" "${INFO_QC_DIR}/0B_after_ids.txt"
        append_state "$QC_STATES" "0B_missing_phenotype" "$WORKBASE" "0B_removed_missing_phenotype_iid.txt"
    else
        echo "# NO performed Step 0 phenotypes"
    fi

    echo "### End of Step 0: raw files script adaptation ###"

    write_ids "$WORKBASE" "${INFO_QC_DIR}/initial_ids.txt"
    write_snps "$WORKBASE" "${INFO_QC_DIR}/initial_snps.txt"

    #===============================================================
    # Step 1: MISSINGNESS
    #===============================================================
    if [[ "$step1" == "NO" ]]; then
        echo "### NO performed Step 1: Missingness ###"
    else
        echo "### Starting Step 1: Missingness ###"
        require_r_script "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"

        "$runplink" --bfile "$WORKBASE" --missing
        Rscript --no-save "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"
        for f in plink.* hist*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/1A_$f"; done
        append_state "$QC_STATES" "1A_missingness_report" "$WORKBASE" "NA"

        write_snps "$WORKBASE" "${INFO_QC_DIR}/1B_before_snps.txt"
        "$runplink" --bfile "$WORKBASE" --geno "$mymissigness1" --make-bed --out "$WORKBASE"
        "$runplink" --bfile "$WORKBASE" --missing
        Rscript --no-save "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"
        for f in plink.* hist*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/1B_$f"; done
        write_snps "$WORKBASE" "${INFO_QC_DIR}/1B_after_snps.txt"
        diff_snps_removed "${INFO_QC_DIR}/1B_before_snps.txt" "${INFO_QC_DIR}/1B_after_snps.txt" "${INFO_QC_DIR}/1B_removed_snp_missingness_gt_${mymissigness1}.txt"
        append_state "$QC_STATES" "1B_geno_${mymissigness1}" "$WORKBASE" "1B_removed_snp_missingness_gt_${mymissigness1}.txt"

        write_ids "$WORKBASE" "${INFO_QC_DIR}/1C_before_ids.txt"
        "$runplink" --bfile "$WORKBASE" --mind "$mymissigness1" --make-bed --out "$WORKBASE"
        "$runplink" --bfile "$WORKBASE" --missing
        Rscript --no-save "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"
        Rscript --no-save "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"
        write_ids "$WORKBASE" "${INFO_QC_DIR}/1C_after_ids.txt"
        diff_ids_removed "${INFO_QC_DIR}/1C_before_ids.txt" "${INFO_QC_DIR}/1C_after_ids.txt" "${INFO_QC_DIR}/1C_removed_ind_missingness_gt_${mymissigness1}.txt"
        append_state "$QC_STATES" "1C_mind_${mymissigness1}" "$WORKBASE" "1C_removed_ind_missingness_gt_${mymissigness1}.txt"

        write_snps "$WORKBASE" "${INFO_QC_DIR}/1D_before_snps.txt"
        "$runplink" --bfile "$WORKBASE" --geno "$mymissigness2" --make-bed --out "$WORKBASE"
        "$runplink" --bfile "$WORKBASE" --missing
        Rscript --no-save "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"
        for f in plink.* hist*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/1D_$f"; done
        write_snps "$WORKBASE" "${INFO_QC_DIR}/1D_after_snps.txt"
        diff_snps_removed "${INFO_QC_DIR}/1D_before_snps.txt" "${INFO_QC_DIR}/1D_after_snps.txt" "${INFO_QC_DIR}/1D_removed_snp_missingness_gt_${mymissigness2}.txt"
        append_state "$QC_STATES" "1D_geno_${mymissigness2}" "$WORKBASE" "1D_removed_snp_missingness_gt_${mymissigness2}.txt"

        write_ids "$WORKBASE" "${INFO_QC_DIR}/1E_before_ids.txt"
        "$runplink" --bfile "$WORKBASE" --mind "$mymissigness2" --make-bed --out "$WORKBASE"
        "$runplink" --bfile "$WORKBASE" --missing
        Rscript --no-save "${QC_SCRIPT_DIR}/QC1_Hist_miss.R"
        for f in plink.* hist*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/1E_$f"; done
        write_ids "$WORKBASE" "${INFO_QC_DIR}/1E_after_ids.txt"
        diff_ids_removed "${INFO_QC_DIR}/1E_before_ids.txt" "${INFO_QC_DIR}/1E_after_ids.txt" "${INFO_QC_DIR}/1E_removed_ind_missingness_gt_${mymissigness2}.txt"
        append_state "$QC_STATES" "1E_mind_${mymissigness2}" "$WORKBASE" "1E_removed_ind_missingness_gt_${mymissigness2}.txt"

        echo "### End of Step 1: Missingness ###"
    fi

    #===============================================================
    # Step 2: SEX DISCREPANCY
    #===============================================================
    if [[ "$step2" == "NO" ]]; then
        echo "### NO performed Step 2: Sex discrepancy ###"
    else
        echo "### Starting Step 2: Sex discrepancy ###"
        require_r_script "${QC_SCRIPT_DIR}/QC2_Gender_check.R"
        SEX_FILE="${INPUT_BATCH_DIR}/${BATCH_NAME}_pheno_sex.txt"

        if [[ -f "$SEX_FILE" ]]; then
            echo "Sex file detected -> updating sex from file"
            "$runplink" --bfile "$WORKBASE" --update-sex "$SEX_FILE" --make-bed --out "$WORKBASE"
        else
            echo "No sex file found -> keeping sex from .fam"
            "$runplink" --bfile "$WORKBASE" --make-bed --out "$WORKBASE"
        fi

        SEXCHECK_DATA="${WORKBASE}_splitx_for_sexcheck"
        if awk '$1 == "25" || $1 == "XY" {found=1} END {exit !found}' "${WORKBASE}.bim"; then
            echo "XY/PAR region already detected -> running --check-sex directly"
            SEXCHECK_BFILE="$WORKBASE"
            SPLITX_INITIAL_STATUS="NO_ALREADY_XY_PAR"
            SPLITX_INITIAL_BUILD="not_applied_configured_${mysplitx_build}"
        else
            echo "No XY/PAR region detected -> creating temporary split-X dataset"
            "$runplink" --bfile "$WORKBASE" --split-x "$mysplitx_build" no-fail --make-bed --out "$SEXCHECK_DATA"
            SEXCHECK_BFILE="$SEXCHECK_DATA"
            SPLITX_INITIAL_STATUS="YES_SPLITX_APPLIED"
            SPLITX_INITIAL_BUILD="$mysplitx_build"
        fi

        "$runplink" --bfile "$SEXCHECK_BFILE" --check-sex
        Rscript --no-save "${QC_SCRIPT_DIR}/QC2_Gender_check.R"
        for f in plink.* Gender_* Men_* Women*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/2_$f"; done
        for f in "${SEXCHECK_DATA}".*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/2_splitx_$(basename "$f")"; done

        grep "PROBLEM" "${INFO_QC_DIR}/2_plink.sexcheck" | awk '{print $1,$2}' > "${INFO_QC_DIR}/2_removed_sex_discrepancy.txt" || true

        write_ids "$WORKBASE" "${INFO_QC_DIR}/2_before_ids.txt"
        "$runplink" --bfile "$WORKBASE" --remove "${INFO_QC_DIR}/2_removed_sex_discrepancy.txt" --make-bed --out "$WORKBASE"
        [[ -f "${WORKBASE}.log" ]] && cp "${WORKBASE}.log" "${INFO_QC_DIR}/2_removal_${WORKDATA}.log"
        write_ids "$WORKBASE" "${INFO_QC_DIR}/2_after_ids.txt"
        append_state "$QC_STATES" "2_sex_discrepancy" "$WORKBASE" "2_removed_sex_discrepancy.txt"

        SEXCHECK_DATA="${WORKBASE}_splitx_for_sexcheck_after_removal"
        if awk '$1 == "25" || $1 == "XY" {found=1} END {exit !found}' "${WORKBASE}.bim"; then
            echo "XY/PAR region already detected -> running final --check-sex directly"
            SEXCHECK_BFILE="$WORKBASE"
            SPLITX_FINAL_STATUS="NO_ALREADY_XY_PAR"
            SPLITX_FINAL_BUILD="not_applied_configured_${mysplitx_build}"
        else
            echo "No XY/PAR region detected -> creating temporary split-X dataset for final check"
            "$runplink" --bfile "$WORKBASE" --split-x "$mysplitx_build" no-fail --make-bed --out "$SEXCHECK_DATA"
            SEXCHECK_BFILE="$SEXCHECK_DATA"
            SPLITX_FINAL_STATUS="YES_SPLITX_APPLIED"
            SPLITX_FINAL_BUILD="$mysplitx_build"
        fi

        "$runplink" --bfile "$SEXCHECK_BFILE" --check-sex
        Rscript --no-save "${QC_SCRIPT_DIR}/QC2_Gender_check.R"
        for f in plink.* Gender_* Men_* Women*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/2_final_$f"; done
        for f in "${SEXCHECK_DATA}".*; do [[ -e "$f" ]] && mv "$f" "${INFO_QC_DIR}/2_final_splitx_$(basename "$f")"; done
        echo "### End of Step 2: Sex discrepancy ###"
    fi

    #===============================================================
    # Step 3: MAF
    #===============================================================
    if [[ "$step3" == "NO" ]]; then
        echo "### NO performed Step 3: minor allele frequency ###"
    else
        echo "### Starting Step 3: minor allele frequency ###"
        require_r_script "${QC_SCRIPT_DIR}/QC3_MAF_check.R"

        write_snps "$WORKBASE" "${INFO_QC_DIR}/3A_before_snps.txt"
        awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' "${WORKBASE}.bim" > "${INFO_QC_DIR}/3_snp_1_22.txt"
        "$runplink" --bfile "$WORKBASE" --extract "${INFO_QC_DIR}/3_snp_1_22.txt" --make-bed --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "3A"
        write_snps "$WORKBASE" "${INFO_QC_DIR}/3A_after_snps.txt"
        diff_snps_removed "${INFO_QC_DIR}/3A_before_snps.txt" "${INFO_QC_DIR}/3A_after_snps.txt" "${INFO_QC_DIR}/3A_removed_non_autosomal_snps.txt"
        "$runplink" --bfile "$WORKBASE" --freq --out "${INFO_QC_DIR}/3A_MAF_check"
        Rscript --no-save "${QC_SCRIPT_DIR}/QC3_MAF_check.R" "$INFO_QC_DIR" "3A"
        append_state "$QC_STATES" "3A_autosomal_only" "$WORKBASE" "3A_removed_non_autosomal_snps.txt"

        write_snps "$WORKBASE" "${INFO_QC_DIR}/3B_before_snps.txt"
        "$runplink" --bfile "$WORKBASE" --maf "$mymaf" --make-bed --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "3B"
        write_snps "$WORKBASE" "${INFO_QC_DIR}/3B_after_snps.txt"
        diff_snps_removed "${INFO_QC_DIR}/3B_before_snps.txt" "${INFO_QC_DIR}/3B_after_snps.txt" "${INFO_QC_DIR}/3B_removed_maf_lt_${mymaf}.txt"
        "$runplink" --bfile "$WORKBASE" --freq --out "${INFO_QC_DIR}/3B_MAF_check"
        Rscript --no-save "${QC_SCRIPT_DIR}/QC3_MAF_check.R" "$INFO_QC_DIR" "3B"
        append_state "$QC_STATES" "3B_maf_${mymaf}" "$WORKBASE" "3B_removed_maf_lt_${mymaf}.txt"

        echo "### End of Step 3: minor allele frequency ###"
    fi

    #===============================================================
    # Step 4: HWE
    #===============================================================
    if [[ "$step4" == "NO" ]]; then
        echo "### NO performed Step 4: hardy-weinberg equilibrium ###"
    else
        echo "### Starting Step 4: hardy-weinberg equilibrium ###"
        require_r_script "${QC_SCRIPT_DIR}/QC4_HWE.R"

        "$runplink" --bfile "$WORKBASE" --hardy --out "${INFO_QC_DIR}/4A_hwe"
        Rscript --no-save "${QC_SCRIPT_DIR}/QC4_HWE.R" "$INFO_QC_DIR" "4A"

        write_snps "$WORKBASE" "${INFO_QC_DIR}/4A_before_snps.txt"
        "$runplink" --bfile "$WORKBASE" --hwe "$myhwe1" --make-bed --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "4A"
        write_snps "$WORKBASE" "${INFO_QC_DIR}/4A_after_snps.txt"
        diff_snps_removed "${INFO_QC_DIR}/4A_before_snps.txt" "${INFO_QC_DIR}/4A_after_snps.txt" "${INFO_QC_DIR}/4A_removed_hwe_controls_lt_${myhwe1}.txt"
        append_state "$QC_STATES" "4A_hwe_controls_${myhwe1}" "$WORKBASE" "4A_removed_hwe_controls_lt_${myhwe1}.txt"

        write_snps "$WORKBASE" "${INFO_QC_DIR}/4B_before_snps.txt"
        "$runplink" --bfile "$WORKBASE" --hwe "$myhwe2" --hwe-all --make-bed --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "4B"
        write_snps "$WORKBASE" "${INFO_QC_DIR}/4B_after_snps.txt"
        diff_snps_removed "${INFO_QC_DIR}/4B_before_snps.txt" "${INFO_QC_DIR}/4B_after_snps.txt" "${INFO_QC_DIR}/4B_removed_hwe_all_lt_${myhwe2}.txt"
        append_state "$QC_STATES" "4B_hwe_all_${myhwe2}" "$WORKBASE" "4B_removed_hwe_all_lt_${myhwe2}.txt"

        echo "### End of Step 4: hardy-weinberg equilibrium ###"
    fi

    #===============================================================
    # Step 5: HETEROZYGOSITY
    #===============================================================
    if [[ "$step5" == "NO" ]]; then
        echo "### NO performed Step 5: heterozygosity ###"
    else
        echo "### Starting Step 5: heterozygosity ###"
        require_r_script "${QC_SCRIPT_DIR}/QC5_Check_heterozygosity_rate.R"

        "$runplink" --bfile "$WORKBASE" \
            --exclude "${QC_SCRIPT_DIR}/QC5_Inversion.txt" --range \
            --indep-pairwise "$mywindowsize" "$myshiftwindow" "$mypairwiser2" \
            --out "${INFO_QC_DIR}/5A_indepSNP" \
            > "${INFO_QC_DIR}/5A_indepSNP.log" 2>&1

        "$runplink" --bfile "$WORKBASE" \
            --extract "${INFO_QC_DIR}/5A_indepSNP.prune.in" \
            --het --out "${INFO_QC_DIR}/5_pruned.data"

        Rscript --no-save "${QC_SCRIPT_DIR}/QC5_Check_heterozygosity_rate.R" "${INFO_QC_DIR}/5_pruned.data.het" "$INFO_QC_DIR"

        sed 's/"//g' "${INFO_QC_DIR}/5_fail-het-qc.txt" | awk '{print $1, $2}' > "${INFO_QC_DIR}/5A_removed_heterozygosity_outliers.txt"
        write_ids "$WORKBASE" "${INFO_QC_DIR}/5_before_ids.txt"
        "$runplink" --bfile "$WORKBASE" \
            --remove "${INFO_QC_DIR}/5A_removed_heterozygosity_outliers.txt" \
            --make-bed \
            --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "5"
        write_ids "$WORKBASE" "${INFO_QC_DIR}/5_after_ids.txt"
        append_state "$QC_STATES" "5_heterozygosity" "$WORKBASE" "5A_removed_heterozygosity_outliers.txt"

        echo "### End of Step 5: heterozygosity ###"
    fi

    #===============================================================
    # Step 6: RELATEDNESS
    # Related individuals are identified but NOT removed.
    #===============================================================
    RELATED_FLAGGED="${INFO_QC_DIR}/${BATCH_NAME}_related_px.txt"
    if [[ "$step6" == "NO" ]]; then
        echo "### NO performed Step 6: relatedness ###"
        printf "FID\tIID\tRELATED_FID\tRELATED_IID\tRT\tPI_HAT\n" > "$RELATED_FLAGGED"
        append_state "$QC_STATES" "6_relatedness_not_run" "$WORKBASE" "${BATCH_NAME}_related_px.txt"
    else
        echo "### Starting Step 6: relatedness ###"
        require_r_script "${QC_SCRIPT_DIR}/QC6_Relatedness.R"

        "$runplink" --bfile "$WORKBASE" \
            --extract "${INFO_QC_DIR}/5A_indepSNP.prune.in" \
            --genome --min "$myrelatedness" \
            --out "${INFO_QC_DIR}/6A_pihat_min${myrelatedness}" \
            > "${INFO_QC_DIR}/6A_pihat.log" 2>&1

        Rscript --no-save "${QC_SCRIPT_DIR}/QC6_Relatedness.R" "${INFO_QC_DIR}/6A_pihat_min${myrelatedness}.genome" "$INFO_QC_DIR"

        "$runplink" --bfile "$WORKBASE" --filter-founders --make-bed --out "$WORKBASE"
        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "6A"

        "$runplink" --bfile "$WORKBASE" \
            --extract "${INFO_QC_DIR}/5A_indepSNP.prune.in" \
            --genome --min "$myrelatedness" \
            --out "${INFO_QC_DIR}/6B_pihat_min${myrelatedness}" \
            > "${INFO_QC_DIR}/6B_pihat.log" 2>&1

        {
            printf "FID\tIID\tRELATED_FID\tRELATED_IID\tRT\tPI_HAT\n"
            awk 'BEGIN {OFS="\t"}
                 NR>1 {
                     print $1,$2,$3,$4,$5,$10
                     print $3,$4,$1,$2,$5,$10
                 }' "${INFO_QC_DIR}/6B_pihat_min${myrelatedness}.genome" | sort -u
        } > "${INFO_QC_DIR}/6B_related_individuals_not_removed.txt"
        cp "${INFO_QC_DIR}/6B_related_individuals_not_removed.txt" "$RELATED_FLAGGED"
        N_RELATED_PAIRS=$(count_related_flagged_pairs "$RELATED_FLAGGED")
        N_RELATED=$(count_related_flagged_individuals "$RELATED_FLAGGED")
        echo "Related pairs flagged (NOT removed): $N_RELATED_PAIRS"
        echo "Related individuals flagged (both members listed, NOT removed): $N_RELATED -> $RELATED_FLAGGED"

        move_output_aux_to_info "$BATCH_OUTPUT_DIR" "$WORKDATA" "$INFO_QC_DIR" "6C"
        append_state "$QC_STATES" "6_relatedness_flagged_only" "$WORKBASE" "${BATCH_NAME}_related_px.txt"

        echo "### End of Step 6: relatedness ###"
    fi

    FINAL_QC_BASE="${QC_OUTPUT_DIR}/${BATCH_NAME}_preimput_QC"
    mv "${WORKBASE}.bed" "${FINAL_QC_BASE}.bed"
    mv "${WORKBASE}.bim" "${FINAL_QC_BASE}.bim"
    mv "${WORKBASE}.fam" "${FINAL_QC_BASE}.fam"

    write_ids "$FINAL_QC_BASE" "${INFO_QC_DIR}/final_ids.txt"
    write_snps "$FINAL_QC_BASE" "${INFO_QC_DIR}/final_snps.txt"
    diff_ids_removed "${INFO_QC_DIR}/initial_ids.txt" "${INFO_QC_DIR}/final_ids.txt" "${INFO_QC_DIR}/${BATCH_NAME}_excluded_individuals_TOTAL.txt"
    diff_snps_removed "${INFO_QC_DIR}/initial_snps.txt" "${INFO_QC_DIR}/final_snps.txt" "${INFO_QC_DIR}/${BATCH_NAME}_excluded_snps_TOTAL.txt"

    make_transition_summary "$QC_STATES" "${INFO_QC_DIR}/${BATCH_NAME}_STEP1_QC_summary.txt" "STEP 1A QC SUMMARY"
    {
        echo ""
        echo "--- RELATEDNESS (STEP 6) ---"
        echo "Related individuals are IDENTIFIED but NOT removed at pre-imputation QC."
        echo "Flagged individuals PI_HAT > ${myrelatedness}: $(count_related_flagged_individuals "$RELATED_FLAGGED")"
        echo "List: $(basename "$RELATED_FLAGGED")"
        echo ""
        echo "--- SEX CHECK / PSEUDOAUTOSOMAL REGION HANDLING ---"
        echo "Original 23/X variants: ${ORIG_N_X}"
        echo "Original 24/Y variants: ${ORIG_N_Y}"
        echo "Original 25/XY/PAR variants: ${ORIG_N_XY_PAR}"
        echo "Original 26/MT/M variants: ${ORIG_N_MT}"
        echo "Configured split-X build: ${mysplitx_build}"
        echo "Initial sex check split-X status: ${SPLITX_INITIAL_STATUS}"
        echo "Initial sex check split-X build used: ${SPLITX_INITIAL_BUILD}"
        echo "Final sex check split-X status: ${SPLITX_FINAL_STATUS}"
        echo "Final sex check split-X build used: ${SPLITX_FINAL_BUILD}"
    } >> "${INFO_QC_DIR}/${BATCH_NAME}_STEP1_QC_summary.txt"

    cleanup_qc_tracking_files "$INFO_QC_DIR"

    rm -f "${BATCH_OUTPUT_DIR}"/*~
    rm -f "${INFO_QC_DIR}"/*~
    mv "${INFO_QC_DIR}"/*.log "$INFO_QC_LOG_DIR/" 2>/dev/null || true

    echo "### STEP 1A_QC completed for ${BATCH_NAME}"
    echo "QC PLINK dataset: ${FINAL_QC_BASE}.bed/.bim/.fam"
    exec 1>&3 2>&4
    exec 3>&- 4>&-

    #################################################################
    # STEP 1B_DP
    #################################################################
    TOTAL_FINAL_SNPS="NA"
    if [[ "$RUN_STEP1B_DP" == "YES" ]]; then
    DP_LOG="${INFO_DP_DIR}/${BATCH_NAME}_STEP2_DP.log"
    exec 3>&1 4>&2
    exec > >(tee -a "$DP_LOG") 2>&1

    echo "### STEP 1B_DP started for ${BATCH_NAME}"
    echo "Input PLINK dataset from STEP 1A_QC: $FINAL_QC_BASE"
    echo "Output directory: $BATCH_OUTPUT_DIR"
    echo "Info directory  : $INFO_DP_DIR"
    echo "Working directory: $DP_WORK_DIR"

    printf "### Script: STEP1_PREIMPUT_QC_DP
### Batch: %s
### STEP 1B_DP parameters:
s0liftoverTarget=%s
s0liftoverReference=%s
mywindowsize=%s
myshiftwindow=%s
mypairwiser2=%s
" "$BATCH_NAME" "$s0liftoverTarget" "$s0liftoverReference" "$mywindowsize" "$myshiftwindow" "$mypairwiser2" \
  > "${INFO_DP_DIR}/${BATCH_NAME}_STEP2_DP_parameters.txt"

    if [[ ! -f "${FINAL_QC_BASE}.bed" || ! -f "${FINAL_QC_BASE}.bim" || ! -f "${FINAL_QC_BASE}.fam" ]]; then
        echo "ERROR: Missing STEP 1A_QC PLINK files for ${BATCH_NAME}: $FINAL_QC_BASE"
        exit 1
    fi

    rm -rf "$DP_WORK_DIR"
    mkdir -p "$DP_WORK_DIR"
    cd "$DP_WORK_DIR"

    DP_STATES="${ROOT_DIR}/${INFO_DP_DIR}/${BATCH_NAME}_STEP2_DP_states.tsv"
    : > "$DP_STATES"

    cp "${ROOT_DIR}/${FINAL_QC_BASE}.bed" QCtargetdata.bed
    cp "${ROOT_DIR}/${FINAL_QC_BASE}.bim" QCtargetdata.bim
    cp "${ROOT_DIR}/${FINAL_QC_BASE}.fam" QCtargetdata.fam
    append_state "$DP_STATES" "0_input_from_STEP1_QC" "QCtargetdata" "NA"

    #---------------------------------------------------------
    # 0A) Remove palindromic SNPs
    #---------------------------------------------------------
    awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || \
         ($5=="C" && $6=="G") || ($5=="G" && $6=="C")' QCtargetdata.bim \
         | cut -f2 > "${ROOT_DIR}/${INFO_DP_DIR}/0A_removed_palindromic_snps.txt"

    "$runplink" --bfile QCtargetdata --exclude "${ROOT_DIR}/${INFO_DP_DIR}/0A_removed_palindromic_snps.txt" --make-bed --out QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0A_$i"; done
    append_state "$DP_STATES" "0A_no_palindromic" "QCtargetdata" "0A_removed_palindromic_snps.txt"

    #---------------------------------------------------------
    # 0B) Assign missing rsIDs
    #---------------------------------------------------------
    "$runplink" --bfile QCtargetdata --set-missing-var-ids @:#:\$1:\$2 --make-bed --out QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0B_$i"; done
    append_state "$DP_STATES" "0B_set_missing_var_ids" "QCtargetdata" "NA"

    #---------------------------------------------------------
    # 0C) Remove duplicated rsIDs
    #---------------------------------------------------------
    "$runplink" --bfile QCtargetdata --write-snplist --out 0B_targetID.snp
    awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' \
        0B_targetID.snp.snplist 0B_targetID.snp.snplist | sort | uniq \
        > "${ROOT_DIR}/${INFO_DP_DIR}/0C_removed_duplicated_rsids.txt"

    "$runplink" --bfile QCtargetdata --exclude "${ROOT_DIR}/${INFO_DP_DIR}/0C_removed_duplicated_rsids.txt" --make-bed --out QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0C_$i"; done
    append_state "$DP_STATES" "0C_no_duplicated_rsids" "QCtargetdata" "0C_removed_duplicated_rsids.txt"

    #---------------------------------------------------------
    # 0D) Remove same position (CH, POS, A1, A2)
    #---------------------------------------------------------
    require_r_script "${ROOT_DIR}/${DP_SCRIPT_DIR}/PRE0D_Duplicates.R"
    "$runplink" --bfile QCtargetdata --list-duplicate-vars --freq --out 0C_targetID.duplicated
    Rscript --no-save "${ROOT_DIR}/${DP_SCRIPT_DIR}/PRE0D_Duplicates.R"
    cp 0C_targetID.duplicated.exclude.txt "${ROOT_DIR}/${INFO_DP_DIR}/0D_removed_duplicated_position_low_frequency.txt" 2>/dev/null || true
    "$runplink" --bfile QCtargetdata --exclude 0C_targetID.duplicated.exclude.txt --make-bed --out QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0D_$i"; done
    append_state "$DP_STATES" "0D_no_same_position" "QCtargetdata" "0D_removed_duplicated_position_low_frequency.txt"

    #---------------------------------------------------------
    # 0E) Set new ID as Chr:BP
    #---------------------------------------------------------
    awk '{print $2, $1":"$4}' QCtargetdata.bim > 0D_targetIDcoordinates.txt
    "$runplink" --bfile QCtargetdata --update-name 0D_targetIDcoordinates.txt --make-bed --out QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0E_$i"; done
    append_state "$DP_STATES" "0E_ids_chr_bp" "QCtargetdata" "NA"

    #---------------------------------------------------------
    # 0F) Remove same rsID after renaming
    #---------------------------------------------------------
    "$runplink" --bfile QCtargetdata --write-snplist --out 0E_targetID.snp
    awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' 0E_targetID.snp.snplist 0E_targetID.snp.snplist \
        | sort | uniq > "${ROOT_DIR}/${INFO_DP_DIR}/0F_removed_duplicated_chr_bp_ids.txt"

    "$runplink" --bfile QCtargetdata --exclude "${ROOT_DIR}/${INFO_DP_DIR}/0F_removed_duplicated_chr_bp_ids.txt" --make-bed --out 0F_QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0F_$i"; done
    append_state "$DP_STATES" "0F_no_duplicated_chr_bp" "0F_QCtargetdata" "0F_removed_duplicated_chr_bp_ids.txt"
    cp 0F_QCtargetdata.bed QCtargetdata.bed
    cp 0F_QCtargetdata.bim QCtargetdata.bim
    cp 0F_QCtargetdata.fam QCtargetdata.fam

    #---------------------------------------------------------
    # 0.3) Optional LiftOver blocks, kept in original position
    #---------------------------------------------------------
    if [[ "$s0liftoverTarget" == "YES" ]]; then
        echo "### Starting Step 0: LiftOver TARGET ###"
        require_r_script "${ROOT_DIR}/${DP_SCRIPT_DIR}/PRE0_FID_IID.R"
        "$runplink" --bfile QCtargetdata --recode --tab --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0A_$i"; done
        awk '{print "chr"$1,$4,$4+1,$2}' 0A_QCtargetdata.map > 0B_QCtargetdata.bed
        chmod +x "$runliftover"
        chmod +x "$liftoverChain"
        "$runliftover" 0B_QCtargetdata.bed "$liftoverChain" 0C_QCtargetdata.bed 0C_unlifted.bed
        awk 'BEGIN {OFS="\t"} ; {print substr($1,4),$4,0,$2}' 0C_QCtargetdata.bed > 0D_QCtargetdata.map
        awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5 || $1==6 || $1==7 || $1==8 || $1==9 || $1==10 || $1==11 || $1==12 || $1==13 || $1==14 || $1==15 || $1==16 || $1==17 || $1==18 || $1==19 || $1==20 || $1==21 || $1==22' 0D_QCtargetdata.map > 0E_QCtargetdata.map
        awk '{print $2,$4}' 0E_QCtargetdata.map > 0E_ID-POS.txt
        awk '{print $2}' 0E_QCtargetdata.map > 0E_ID.txt
        "$runplink" --bfile QCtargetdata --extract 0E_ID.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0F_$i"; done
        "$runplink" --bfile QCtargetdata --update-map 0E_ID-POS.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0G_$i"; done
        awk '{print $2, $1":"$4}' QCtargetdata.bim > 0G_IDcoordinates.txt
        "$runplink" --bfile QCtargetdata --update-name 0G_IDcoordinates.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0H_$i"; done
        Rscript --no-save "${ROOT_DIR}/${DP_SCRIPT_DIR}/PRE0_FID_IID.R"
        "$runplink" --bfile QCtargetdata --update-ids 0H_FID.IID.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0I_$i"; done
        echo "### End of Step 0: LiftOver TARGET ###"
    elif [[ "$s0liftoverReference" == "YES" ]]; then
        echo "### Starting Step 0: LiftOver REFERENCE ###"
        echo "WARNING: Reference LiftOver block is kept from original script but QCreferencedata is not part of this combined target-data pipeline."
    fi

    #================================================================
    # STEP 5: Allele update according TopMed / HRC
    #================================================================
    "$runplink" --bfile 0F_QCtargetdata --freq --out QCtargetdata
    for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "5A_$i"; done

    wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
    unzip HRC-1000G-check-bim-v4.3.0.zip
    rm -f HRC-1000G-check-bim-v4.3.0.zip

    wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
    gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

    perl HRC-1000G-check-bim.pl -b QCtargetdata.bim -f 5A_QCtargetdata.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
    sh Run-plink.sh
    echo "### End of Step 5: allele update according Topmed ###"

    #================================================================
    # STEP 6: Split by chromosome & VCF creation
    #================================================================
    TOTAL_FINAL_SNPS=0
    for chr in {1..22}; do
        bcftools sort "QCtargetdata-updated-chr${chr}.vcf" -Oz -o "${ROOT_DIR}/${DP_OUTPUT_DIR}/QCtargetdata-updated-chr${chr}.vcf.gz"
        bcftools index --tbi "${ROOT_DIR}/${DP_OUTPUT_DIR}/QCtargetdata-updated-chr${chr}.vcf.gz"
        NSNPS=$(bcftools index -n "${ROOT_DIR}/${DP_OUTPUT_DIR}/QCtargetdata-updated-chr${chr}.vcf.gz")
        TOTAL_FINAL_SNPS=$((TOTAL_FINAL_SNPS + NSNPS))
    done

    #================================================================
    # STEP 7/8: DP summaries
    #================================================================
    SUMMARY_FILE="${ROOT_DIR}/${INFO_DP_DIR}/${BATCH_NAME}_STEP2_DP_summary.txt"
    make_transition_summary "$DP_STATES" "$SUMMARY_FILE" "STEP 1B DP SUMMARY"

    {
        echo ""
        echo "--- FINAL VCF COUNTS ---"
        printf "%-10s %-15s\n" "CHR" "N_SNPS"
        for chr in {1..22}; do
            VCF="${ROOT_DIR}/${DP_OUTPUT_DIR}/QCtargetdata-updated-chr${chr}.vcf.gz"
            if [[ -f "$VCF" ]]; then
                NSNPS=$(bcftools index -n "$VCF")
            else
                NSNPS=0
            fi
            printf "%-10s %-15s\n" "chr${chr}" "$NSNPS"
        done
        printf "%-10s %-15s\n" "TOTAL" "$TOTAL_FINAL_SNPS"
    } >> "$SUMMARY_FILE"

    echo "### Organising STEP 1B_DP info files ###"
    mkdir -p "${ROOT_DIR}/${INFO_DP_DIR}/logs" "${ROOT_DIR}/${INFO_DP_DIR}/QC_files" "${ROOT_DIR}/${INFO_DP_DIR}/HRC_helpers"
    mv *QCtargetdata-HRC.txt "${ROOT_DIR}/${INFO_DP_DIR}/QC_files/" 2>/dev/null || true
    mv LOG-QCtargetdata-HRC.txt "${ROOT_DIR}/${INFO_DP_DIR}/logs/" 2>/dev/null || true
    mv Strand-Flip-QCtargetdata-HRC.txt "${ROOT_DIR}/${INFO_DP_DIR}/QC_files/" 2>/dev/null || true
    mv 0A_* 0B_* 0C_* 0D_* 0E_* 0F_* 5A_QCtargetdata.frq "${ROOT_DIR}/${INFO_DP_DIR}/QC_files/" 2>/dev/null || true
    mv LICENSE.txt Run-plink.sh "${ROOT_DIR}/${INFO_DP_DIR}/HRC_helpers/" 2>/dev/null || true
    mv *.log "${ROOT_DIR}/${INFO_DP_DIR}/logs/" 2>/dev/null || true
    rm -f *.bed *.bim *.fam QCtargetdata.* QCtargetdata-updated-chr*.vcf HRC.r1-1.GRCh37.wgs.mac5.sites.tab HRC-1000G-check-bim.pl *~

    cd "$ROOT_DIR"
    rm -rf "$DP_WORK_DIR"

    COMBINED_SUMMARY="${BATCH_INFO_DIR}/${BATCH_NAME}_preimput_combined_summary.txt"
    {
        echo "==============================================================="
        echo "COMBINED PREIMPUT SUMMARY"
        echo "==============================================================="
        echo "Batch: $BATCH_NAME"
        echo "Date : $(date)"
        echo ""
        echo "STEP 1A QC summary : _info_STEP1_QC/${BATCH_NAME}_STEP1_QC_summary.txt"
        echo "STEP 1B DP summary : _info_STEP2_DP/${BATCH_NAME}_STEP2_DP_summary.txt"
        echo "Relatedness list   : _info_STEP1_QC/${BATCH_NAME}_related_px.txt"
        echo "QC PLINK dataset   : ../../${QC_OUTPUT_DIR}/${BATCH_NAME}_preimput_QC.bed/.bim/.fam"
        echo "Final VCF outputs  : ../../${DP_OUTPUT_DIR}/QCtargetdata-updated-chr1-22.vcf.gz(.tbi)"
        echo ""
        echo "STEP 1A final samples: $(count_samples "$FINAL_QC_BASE")"
        echo "STEP 1A final SNPs   : $(count_snps "$FINAL_QC_BASE")"
        echo "STEP 1B final SNPs   : $TOTAL_FINAL_SNPS"
        echo "==============================================================="
    } > "$COMBINED_SUMMARY"

    echo "### STEP 1B_DP completed for ${BATCH_NAME}"
    echo "Final VCF outputs saved in: ${DP_OUTPUT_DIR}"
    echo "Combined info saved in: ${BATCH_INFO_DIR}"

    exec 1>&3 2>&4
    exec 3>&- 4>&-
    else
        echo "### STEP 1B_DP skipped for ${BATCH_NAME}"
        {
            echo "==============================================================="
            echo "STEP 1B DP SUMMARY"
            echo "==============================================================="
            echo "Batch: $BATCH_NAME"
            echo "Status: SKIPPED by user"
            echo "Input available for later steps: ${FINAL_QC_BASE}.bed/.bim/.fam"
            echo "==============================================================="
        } > "${INFO_DP_DIR}/${BATCH_NAME}_STEP2_DP_summary.txt"
    fi

    COMBINED_SUMMARY="${BATCH_INFO_DIR}/${BATCH_NAME}_preimput_combined_summary.txt"
    if [[ ! -f "$COMBINED_SUMMARY" ]]; then
        {
            echo "==============================================================="
            echo "COMBINED PREIMPUT SUMMARY"
            echo "==============================================================="
            echo "Batch: $BATCH_NAME"
            echo "Date : $(date)"
            echo ""
            echo "STEP 1A QC summary : _info_STEP1_QC/${BATCH_NAME}_STEP1_QC_summary.txt"
            echo "STEP 1B DP summary : _info_STEP2_DP/${BATCH_NAME}_STEP2_DP_summary.txt"
            echo "Relatedness list   : _info_STEP1_QC/${BATCH_NAME}_related_px.txt"
            echo "QC PLINK dataset   : ../../${QC_OUTPUT_DIR}/${BATCH_NAME}_preimput_QC.bed/.bim/.fam"
            echo ""
            echo "STEP 1A final samples: $(count_samples "$FINAL_QC_BASE")"
            echo "STEP 1A final SNPs   : $(count_snps "$FINAL_QC_BASE")"
            echo "STEP 1B status       : SKIPPED"
            echo "==============================================================="
        } > "$COMBINED_SUMMARY"
    fi

    #################################################################
    # STEP 1C_PCA
    #################################################################
    if [[ "$RUN_STEP1C_PCA" == "YES" ]]; then
        PCA_LOG_DIR="${INFO_PCA_DIR}/logs"
        PCA_QC_DIR="${INFO_PCA_DIR}/QC_files"
        PCA_RESULTS_DIR="${INFO_PCA_DIR}/PCA_results"
        mkdir -p "$PCA_LOG_DIR" "$PCA_QC_DIR" "$PCA_RESULTS_DIR"

        PCA_LOG="${PCA_LOG_DIR}/${BATCH_NAME}_STEP3_PCA.log"
        exec 3>&1 4>&2
        exec > >(tee -a "$PCA_LOG") 2>&1

        echo "### STEP 1C_PCA started for ${BATCH_NAME}"
        echo "Target input from STEP 1A_QC: $FINAL_QC_BASE"
        echo "Reference input directory: $PCA_REF_DIR"
        echo "Info directory: $INFO_PCA_DIR"
        echo "Working directory: $PCA_WORK_DIR"

        require_r_script "${PCA_SCRIPT_DIR}/PRE0D_Duplicates.R"
        require_r_script "${PCA_SCRIPT_DIR}/PRE3_PCA.R"

        if [[ ! -f "${FINAL_QC_BASE}.bed" || ! -f "${FINAL_QC_BASE}.bim" || ! -f "${FINAL_QC_BASE}.fam" ]]; then
            echo "ERROR: Missing STEP 1A_QC PLINK files for PCA: $FINAL_QC_BASE"
            exit 1
        fi

        if [[ -f "${PCA_REF_DIR}/1000G_QC_PRC_v4.bed" && -f "${PCA_REF_DIR}/1000G_QC_PRC_v4.bim" && -f "${PCA_REF_DIR}/1000G_QC_PRC_v4.fam" ]]; then
            REF_BASE="${PCA_REF_DIR}/1000G_QC_PRC_v4"
        elif [[ -f "${PCA_REF_DIR}/QCreferencedata.bed" && -f "${PCA_REF_DIR}/QCreferencedata.bim" && -f "${PCA_REF_DIR}/QCreferencedata.fam" ]]; then
            REF_BASE="${PCA_REF_DIR}/QCreferencedata"
        else
            mapfile -t REF_BEDS < <(find "$PCA_REF_DIR" -maxdepth 1 -type f \( -name '1000G_QC_PRC*.bed' -o -name 'QCreferencedata*.bed' \) | sort)
            if [[ ${#REF_BEDS[@]} -eq 0 ]]; then
                echo "ERROR: No reference PLINK dataset detected in ${PCA_REF_DIR}"
                echo "Expected 1000G_QC_PRC_v4.bed/.bim/.fam or QCreferencedata.bed/.bim/.fam"
                exit 1
            fi
            REF_BASE="${REF_BEDS[0]%.bed}"
        fi

        if [[ ! -f "${REF_BASE}.bed" || ! -f "${REF_BASE}.bim" || ! -f "${REF_BASE}.fam" ]]; then
            echo "ERROR: Incomplete reference PLINK dataset: $REF_BASE"
            exit 1
        fi

        if [[ -f "${PCA_REF_DIR}/integrated_call_samples_v3.20200731.ALL.ped" ]]; then
            PED_REF="${PCA_REF_DIR}/integrated_call_samples_v3.20200731.ALL.ped"
        else
            PED_REF="${PCA_REF_DIR}/PEDreferencedata_1000G.ped"
        fi
        if [[ ! -f "$PED_REF" ]]; then
            echo "ERROR: Missing PCA pedigree/population helper: $PED_REF"
            echo "Expected integrated_call_samples_v3.20200731.ALL.ped or PEDreferencedata_1000G.ped in ${PCA_REF_DIR}"
            exit 1
        fi

        rm -rf "$PCA_WORK_DIR"
        mkdir -p "$PCA_WORK_DIR"
        cd "$PCA_WORK_DIR"

        cp "${ROOT_DIR}/${FINAL_QC_BASE}.bed" QCtargetdata.bed
        cp "${ROOT_DIR}/${FINAL_QC_BASE}.bim" QCtargetdata.bim
        cp "${ROOT_DIR}/${FINAL_QC_BASE}.fam" QCtargetdata.fam
        cp "${ROOT_DIR}/${REF_BASE}.bed" QCreferencedata.bed
        cp "${ROOT_DIR}/${REF_BASE}.bim" QCreferencedata.bim
        cp "${ROOT_DIR}/${REF_BASE}.fam" QCreferencedata.fam
        cp "${ROOT_DIR}/${PED_REF}" PEDreferencedata_1000G.ped

        echo "### Starting STEP 1C.0: PCA bfiles adaptation ###"

        awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C")' QCtargetdata.bim | cut -f2 > 0A_palindromic_snps.txt
        "$runplink" --bfile QCtargetdata --exclude 0A_palindromic_snps.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0A_$i"; done

        "$runplink" --bfile QCtargetdata --set-missing-var-ids @:#:\$1:\$2 --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0B_$i"; done
        "$runplink" --bfile QCreferencedata --set-missing-var-ids @:#:\$1:\$2 --make-bed --out QCreferencedata
        for i in $(find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0B_$i"; done

        "$runplink" --bfile QCtargetdata --write-snplist --out 0B_targetID.snp
        awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' 0B_targetID.snp.snplist 0B_targetID.snp.snplist | sort | uniq > 0B_targetID.duplicated.txt
        "$runplink" --bfile QCtargetdata --exclude 0B_targetID.duplicated.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0C_$i"; done

        "$runplink" --bfile QCreferencedata --write-snplist --out 0B_referenceID.snp
        awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' 0B_referenceID.snp.snplist 0B_referenceID.snp.snplist | sort | uniq > 0B_referenceID.duplicated.txt
        "$runplink" --bfile QCreferencedata --exclude 0B_referenceID.duplicated.txt --make-bed --out QCreferencedata
        for i in $(find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0C_$i"; done

        "$runplink" --bfile QCtargetdata --list-duplicate-vars --freq --out 0C_targetID.duplicated
        "$runplink" --bfile QCreferencedata --list-duplicate-vars --freq --out 0C_referenceID.duplicated
        Rscript --no-save "${ROOT_DIR}/${PCA_SCRIPT_DIR}/PRE0D_Duplicates.R"

        "$runplink" --bfile QCtargetdata --exclude 0C_targetID.duplicated.exclude.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0D_$i"; done
        "$runplink" --bfile QCreferencedata --exclude 0C_referenceID.duplicated.exclude.txt --make-bed --out QCreferencedata
        for i in $(find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0D_$i"; done

        awk '{print $2, $1":"$4}' QCtargetdata.bim > 0D_targetIDcoordinates.txt
        "$runplink" --bfile QCtargetdata --update-name 0D_targetIDcoordinates.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0E_$i"; done
        awk '{print $2, $1":"$4}' QCreferencedata.bim > 0D_referenceIDcoordinates.txt
        "$runplink" --bfile QCreferencedata --update-name 0D_referenceIDcoordinates.txt --make-bed --out QCreferencedata
        for i in $(find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0E_$i"; done

        "$runplink" --bfile QCtargetdata --write-snplist --out 0E_targetID.snp
        awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' 0E_targetID.snp.snplist 0E_targetID.snp.snplist | sort | uniq > 0E_targetID.duplicated.txt
        "$runplink" --bfile QCtargetdata --exclude 0E_targetID.duplicated.txt --make-bed --out 0F_QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0F_$i"; done
        cp 0F_QCtargetdata.bed QCtargetdata.bed
        cp 0F_QCtargetdata.bim QCtargetdata.bim
        cp 0F_QCtargetdata.fam QCtargetdata.fam

        "$runplink" --bfile QCreferencedata --write-snplist --out 0E_referenceID.snp
        awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' 0E_referenceID.snp.snplist 0E_referenceID.snp.snplist | sort | uniq > 0E_referenceID.duplicated.txt
        "$runplink" --bfile QCreferencedata --exclude 0E_referenceID.duplicated.txt --make-bed --out QCreferencedata
        for i in $(find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "0F_$i"; done

        echo "### End STEP 1C.0: PCA bfiles adaptation ###"

        echo "### Starting STEP 1C.1: prepare target and reference for merge ###"
        awk '{print $2,$5}' QCreferencedata.bim > 1A_QCreferencedata-list.txt
        "$runplink" --bfile 0F_QCtargetdata --reference-allele 1A_QCreferencedata-list.txt --make-bed --out QCtargetdata \
            > 1A_reference_allele.log 2>&1
        grep "Impossible A1 allele assignment" 1A_reference_allele.log > 1A_reference_allele_impossible_A1.txt || true
        N_A1_WARN_1A=$(wc -l < 1A_reference_allele_impossible_A1.txt | awk '{print $1}')
        echo "Reference allele assignment warnings before flip: ${N_A1_WARN_1A} (details saved in PCA QC files)"
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "1A_$i"; done

        awk '{print $2,$5,$6}' QCreferencedata.bim > 1A_QCreferencedata_strand.txt
        awk '{print $2,$5,$6}' QCtargetdata.bim > 1A_QCtargetdata_strand.txt
        sort 1A_QCreferencedata_strand.txt 1A_QCtargetdata_strand.txt | uniq -u > 1A_strand_differences.txt
        awk '{print $1}' 1A_strand_differences.txt | sort -u > 1A_flip_list.txt
        "$runplink" --bfile QCtargetdata --flip 1A_flip_list.txt --reference-allele 1A_QCreferencedata-list.txt --make-bed --out QCtargetdata \
            > 1B_flip_reference_allele.log 2>&1
        grep "Impossible A1 allele assignment" 1B_flip_reference_allele.log > 1B_flip_reference_allele_impossible_A1.txt || true
        N_A1_WARN_1B=$(wc -l < 1B_flip_reference_allele_impossible_A1.txt | awk '{print $1}')
        echo "Reference allele assignment warnings after flip: ${N_A1_WARN_1B} (details saved in PCA QC files)"
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "1B_$i"; done

        awk '{print $2,$5,$6}' QCtargetdata.bim > 1B_QCtargetdata_problematics.txt
        sort 1A_QCreferencedata_strand.txt 1B_QCtargetdata_problematics.txt | uniq -u > 1B_QCtargetdata_uncorresponded.txt
        awk '{print $1}' 1B_QCtargetdata_uncorresponded.txt | sort -u > 1B_QCtargetdata_excluded.txt

        "$runplink" --bfile QCtargetdata --exclude 1B_QCtargetdata_excluded.txt --make-bed --out QCtargetdata
        for i in $(find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "1C_$i"; done
        "$runplink" --bfile QCreferencedata --exclude 1B_QCtargetdata_excluded.txt --make-bed --out QCreferencedata
        for i in $(find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "1C_$i"; done
        echo "### End STEP 1C.1: prepare target and reference for merge ###"

        echo "### Starting STEP 1C.2: merge and PCA ###"
        awk '{print $1,$2}' QCtargetdata.fam > 1C_QCtargetdata_individualsID.txt
        awk '{print $1,$2}' QCreferencedata.fam > 1C_QCreferencedata_individualsID.txt
        "$runplink" --bfile QCtargetdata --bmerge QCreferencedata.bed QCreferencedata.bim QCreferencedata.fam --allow-no-sex --make-bed --out PCA_datasets
        for i in $(find . -name 'PCA_datasets.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/'); do mv "$i" "2_$i"; done

        "$runplink" --bfile PCA_datasets --indep-pairwise "$mywindowsize" "$myshiftwindow" "$mypairwiser2" --out 3A_PCA_datasets
        "$runplink" --bfile PCA_datasets --extract 3A_PCA_datasets.prune.in --make-bed --out 3B_PCA_datasets
        "$runplink" --bfile 3B_PCA_datasets --pca --out 3_PCA_results
        Rscript --no-save "${ROOT_DIR}/${PCA_SCRIPT_DIR}/PRE3_PCA.R"

        "$runplink" --bfile 0F_QCtargetdata --indep-pairwise "$mywindowsize" "$myshiftwindow" "$mypairwiser2" --out 4_QCtargetdata
        "$runplink" --bfile 0F_QCtargetdata --extract 4_QCtargetdata.prune.in --make-bed --out 4_QCtargetdata_PCA
        "$runplink" --bfile 4_QCtargetdata_PCA --pca 10 --out 4_PCA_results
        cp 4_PCA_results.eigenvec 4_PCA_results.covariate
        cp 4_PCA_results.covariate "${ROOT_DIR}/${PCA_OUTPUT_DIR}/${BATCH_NAME}_PCA_results.covariate"
        echo "### End STEP 1C.2: merge and PCA ###"

        SUMMARY_FILE="${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_STEP3_PCA_summary.txt"
        {
            echo "==============================================================="
            echo "STEP 1C PCA SUMMARY"
            echo "==============================================================="
            echo "Batch: $BATCH_NAME"
            echo "Date : $(date)"
            echo "Target input: ${FINAL_QC_BASE}.bed/.bim/.fam"
            echo "Reference input: ${REF_BASE}.bed/.bim/.fam"
            echo "Reference PED: ${PED_REF}"
            echo ""
            echo "Target samples after PCA prep: $(wc -l < QCtargetdata.fam)"
            echo "Target SNPs after PCA prep: $(wc -l < QCtargetdata.bim)"
            echo "Reference samples after PCA prep: $(wc -l < QCreferencedata.fam)"
            echo "Reference SNPs after PCA prep: $(wc -l < QCreferencedata.bim)"
            echo "Merged PCA samples: $(wc -l < PCA_datasets.fam)"
            echo "Merged PCA SNPs: $(wc -l < PCA_datasets.bim)"
            echo "Pruned merged PCA SNPs: $(wc -l < 3B_PCA_datasets.bim)"
            echo "Reference allele warnings before flip: ${N_A1_WARN_1A:-NA}"
            echo "Reference allele warnings after flip: ${N_A1_WARN_1B:-NA}"
            echo ""
            echo "Target nonEUR by EUR +/- 3SD: $(wc -l < 3_target_nonEUR_3SD.txt 2>/dev/null || echo 0)"
            echo "Target nonEUR by centroid 15pct: $(wc -l < 3_target_nonEUR_centroid_15pct.txt 2>/dev/null || echo 0)"
            echo "Target uncertain by centroid 15pct: $(wc -l < 3_target_uncertain_centroid_15pct.txt 2>/dev/null || echo 0)"
            echo "Target covariate output: ../../${PCA_OUTPUT_DIR}/${BATCH_NAME}_PCA_results.covariate"
            echo "==============================================================="
        } > "$SUMMARY_FILE"

        cp 3_target_ancestry_assignments.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_ancestry_assignments.txt" 2>/dev/null || true
        cp 3_target_nonEUR_3SD.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_nonEUR_3SD.txt" 2>/dev/null || true
        cp 3_target_EUR_3SD.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_EUR_3SD.txt" 2>/dev/null || true
        cp 3_target_nonEUR_centroid_15pct.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_nonEUR_centroid_15pct.txt" 2>/dev/null || true
        cp 3_target_EUR_centroid_15pct.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_EUR_centroid_15pct.txt" 2>/dev/null || true
        cp 3_target_uncertain_centroid_15pct.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_uncertain_centroid_15pct.txt" 2>/dev/null || true
        cp 1B_QCtargetdata_excluded.txt "${ROOT_DIR}/${INFO_PCA_DIR}/${BATCH_NAME}_PCA_target_reference_uncorresponded_snps.txt" 2>/dev/null || true

        mv 3_PCA_results.* 4_PCA_results.* 3_*.pdf "${ROOT_DIR}/${PCA_RESULTS_DIR}/" 2>/dev/null || true
        mv *.log "${ROOT_DIR}/${PCA_LOG_DIR}/" 2>/dev/null || true
        mv *.txt *.frq *.hh *.nosex *.prune.in *.prune.out "${ROOT_DIR}/${PCA_QC_DIR}/" 2>/dev/null || true

        cd "$ROOT_DIR"
        rm -rf "$PCA_WORK_DIR"

        {
            echo ""
            echo "STEP 1C PCA summary : _info_STEP3_PCA/${BATCH_NAME}_STEP3_PCA_summary.txt"
            echo "STEP 1C PCA covariate: ../../${PCA_OUTPUT_DIR}/${BATCH_NAME}_PCA_results.covariate"
        } >> "${BATCH_INFO_DIR}/${BATCH_NAME}_preimput_combined_summary.txt"

        echo "### STEP 1C_PCA completed for ${BATCH_NAME}"
        exec 1>&3 2>&4
        exec 3>&- 4>&-
    else
        echo "### STEP 1C_PCA skipped for ${BATCH_NAME}"
        {
            echo "==============================================================="
            echo "STEP 1C PCA SUMMARY"
            echo "==============================================================="
            echo "Batch: $BATCH_NAME"
            echo "Status: SKIPPED by user"
            echo "==============================================================="
        } > "${INFO_PCA_DIR}/${BATCH_NAME}_STEP3_PCA_summary.txt"
    fi
done

echo "### PIPELINE SUMMARY ###"
echo "### STEP 1A_QC completed for selected batch(es)."

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
