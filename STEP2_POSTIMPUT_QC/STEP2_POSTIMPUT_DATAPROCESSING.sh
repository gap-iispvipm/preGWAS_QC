#!/bin/bash

###################################################################
######################### START ANALYSIS ##########################
###################################################################
################### Andrea Perez 2026-06-15 #######################
###################################################################
# ===================================
# MODULES NEEDED
# ===================================

module load PLINK/1.9b_6.21-x86_64
module load BCFtools/1.18
module load R/4.3.2-gfbf-2023a
module load p7zip/17.04-GCCcore-11.3.0

    # In case something fails the script stops 
set -euo pipefail

# ===================================
# USER CONFIG
# ===================================

runplink="plink"

    #Thresholds
s4mymaf=0.001
s5mymissingness1=0.2
s5mymissingness2=0.02

    #Declare directories 
INPUT_DIR="A_input_postimput_processing"
INFO_BASE_DIR="info_postimput_processing"
OUTPUT_BASE_DIR="output_postimput_processing"
TRASH_BASE_DIR="trash"
PASSWORD_FILE="A_input_postimput_processing/batch_passwords.txt"

# ===================================
# LOAD PASSWORDS
# FORMAT(WARNING - spaces or tabs):
# BATCH=password
# ===================================

declare -A BATCH_PASSWORDS

while IFS='=' read -r batch pass || [ -n "$batch" ]; do

    batch=$(echo "$batch" | tr -d '[:space:]')
    pass=$(echo "$pass" | tr -d '\r')

    [ -z "$batch" ] && continue

    BATCH_PASSWORDS["$batch"]="$pass"

done < "$PASSWORD_FILE"

# ===================================
# DETECT AVAILABLE BATCHES
# ===================================

echo "Detecting batches inside: ${INPUT_DIR}"

AVAILABLE_BATCHES=()

while IFS= read -r d; do
    base=$(basename "$d")

    batch="${base%_input_postimput_processing}"

    AVAILABLE_BATCHES+=("$batch")

done < <(
    find "${INPUT_DIR}" -maxdepth 1 -type d \
    -name "*_input_postimput_processing"
)

echo "Detected batches:"
printf ' - %s\n' "${AVAILABLE_BATCHES[@]}"

if [ ${#AVAILABLE_BATCHES[@]} -eq 0 ]; then
    echo "ERROR: no batches detected."
    exit 1
fi

echo ""
echo "======================================================="
echo "AVAILABLE BATCHES"
echo "======================================================="

printf " - %s\n" "${AVAILABLE_BATCHES[@]}"

echo ""
echo "======================================================="
echo "PROCESSING MODE"
echo "======================================================="
echo "1) All batches"
echo "2) Selected batches"

read -p "Select option: " PROCESS_MODE

# ===================================
# SELECT BATCHES
# ===================================

SELECTED_BATCHES=()

if [ "$PROCESS_MODE" == "1" ]; then
    SELECTED_BATCHES=("${AVAILABLE_BATCHES[@]}")

    

elif [ "$PROCESS_MODE" == "2" ]; then
    read -p "Enter batch name(s) separated by spaces: " -a INPUT_BATCHES

    for b in "${INPUT_BATCHES[@]}"; do
        if [[ " ${AVAILABLE_BATCHES[*]} " == *" ${b} "* ]]; then
            SELECTED_BATCHES+=("$b")
        else
            echo "WARNING: batch '${b}' not found among detected batches, skipping."
        fi
    done

    if [ ${#SELECTED_BATCHES[@]} -eq 0 ]; then
        echo "ERROR: no valid batches selected."
        exit 1
    fi

    
else
    echo "Invalid option"
    exit 1
fi

###################################################################
################### MAIN  FUNCTION ################################
###################################################################

run_pipeline() {

BATCH_NAME="$1"

echo ""
echo "#######################################################"
echo "### PROCESSING BATCH: ${BATCH_NAME}"
echo "#######################################################"

# ===================================
# PASSWORD
# ===================================

PASSWORD="${BATCH_PASSWORDS[$BATCH_NAME]:-}"

if [ -z "$PASSWORD" ]; then
    echo "ERROR: no password found for ${BATCH_NAME}"
    return
fi

# ===================================
# DIRECTORIES
# ===================================

D="${INPUT_DIR}/${BATCH_NAME}_input_postimput_processing"

O="${OUTPUT_BASE_DIR}/${BATCH_NAME}_output_postimput_processing"

I="${INFO_BASE_DIR}/${BATCH_NAME}_info_postimput_processing"

L="${I}/${BATCH_NAME}_step2_logs"

Q="${I}/${BATCH_NAME}_step2_QC"

T="${TRASH_BASE_DIR}/${BATCH_NAME}_trash"

mkdir -p "$O"
mkdir -p "$I"
mkdir -p "$L"
mkdir -p "$Q"
mkdir -p "$T"


# ===================================
# VALIDATION
# ===================================

if [ ! -d "$D" ]; then
    echo "ERROR: input directory not found:"
    echo "$D"
    return
fi

# ===================================
# STEP 0: UNZIP
#===================================

echo ""
echo "======================================================="
echo "STEP 0: UNZIPPING CHROMOSOMES"
echo "======================================================="
for chr in {1..22}; do

    ZIPFILE="${D}/chr_${chr}.zip"

    if [ ! -f "$ZIPFILE" ]; then
        echo "WARNING: missing ${ZIPFILE}"
        continue
    fi
    echo "Unzipping chr_${chr}.zip"

    7z x "$ZIPFILE" -o"$D" -p"$PASSWORD" -y

done


echo "STEP 0 completed"

# ===================================
# STEP 0B: EXTRACT VCF HEADER INFO
# ===================================

echo ""
echo "======================================================="
echo "STEP 0B: EXTRACTING VCF HEADER CHARACTERISTICS"
echo "======================================================="

# Use chr1 as reference; fallback to first available chromosome
HEADER_SOURCE=""
for chr in {1..22}; do
    if [ -f "${D}/chr${chr}.dose.vcf.gz" ]; then
        HEADER_SOURCE="${D}/chr${chr}.dose.vcf.gz"
        echo "Using chr${chr}.dose.vcf.gz for header extraction"
        break
    fi
done

if [ -z "$HEADER_SOURCE" ]; then
    echo "WARNING: no VCF found for header extraction"
else
    bcftools view -h "$HEADER_SOURCE" | grep '^##' \
        > "${Q}/${BATCH_NAME}_TOPMED_characteristics.txt"

    echo "Header characteristics saved to: ${Q}/${BATCH_NAME}_TOPMED_characteristics.txt"
fi

# ===================================
# STEP 0C: CHROMOSOME DIMENSIONS
# ===================================

echo ""
echo "======================================================="
echo "STEP 0C: CHROMOSOME DIMENSIONS (pre-merge)"
echo "======================================================="

CHR_DIM_FILE="${Q}/${BATCH_NAME}_chr_dim_postimputed.txt"

printf "%-6s %15s %15s\n" "CHR" "Length_bp" "N_variants" > "$CHR_DIM_FILE"
echo "--------------------------------------" >> "$CHR_DIM_FILE"

TOTAL_VARIANTS=0

for chr in {1..22}; do

    VCFFILE="${D}/chr${chr}.dose.vcf.gz"

    if [ ! -f "$VCFFILE" ]; then
        printf "%-6s %15s %15s\n" "$chr" "NA" "NA" >> "$CHR_DIM_FILE"
        continue
    fi

    # Index if not already indexed
    if [ ! -f "${VCFFILE}.tbi" ] && [ ! -f "${VCFFILE}.csi" ]; then
        bcftools index "$VCFFILE"
    fi

    # Length in bp: from ##contig header line
    CHR_LEN=$(bcftools view -h "$VCFFILE" \
        | grep "##contig=<ID=${chr}," \
        | grep -oP 'length=\K[0-9]+' \
        || echo "NA")

    # Number of variants: from index (fast, no decompression)
    N_VAR=$(bcftools index -n "$VCFFILE")

    printf "%-6s %15s %15s\n" "$chr" "$CHR_LEN" "$N_VAR" >> "$CHR_DIM_FILE"

   TOTAL_VARIANTS=$(( TOTAL_VARIANTS + N_VAR ))

done

echo "--------------------------------------" >> "$CHR_DIM_FILE"
printf "%-6s %15s %15s\n" "TOTAL" "-" "$TOTAL_VARIANTS" >> "$CHR_DIM_FILE"

echo "Chromosome dimensions saved to: $CHR_DIM_FILE"
cat "$CHR_DIM_FILE"



# ===================================
# STEP 1: MERGE CHROMOSOMES
# ===================================

echo ""
echo "======================================================="
echo "STEP 1: MERGING CHROMOSOMES"
echo "======================================================="

bcftools concat -O z -o "${D}/1_AllChromosomes.vcf.gz" \
    "${D}/chr1.dose.vcf.gz" \
    "${D}/chr2.dose.vcf.gz" \
    "${D}/chr3.dose.vcf.gz" \
    "${D}/chr4.dose.vcf.gz" \
    "${D}/chr5.dose.vcf.gz" \
    "${D}/chr6.dose.vcf.gz" \
    "${D}/chr7.dose.vcf.gz" \
    "${D}/chr8.dose.vcf.gz" \
    "${D}/chr9.dose.vcf.gz" \
    "${D}/chr10.dose.vcf.gz" \
    "${D}/chr11.dose.vcf.gz" \
    "${D}/chr12.dose.vcf.gz" \
    "${D}/chr13.dose.vcf.gz" \
    "${D}/chr14.dose.vcf.gz" \
    "${D}/chr15.dose.vcf.gz" \
    "${D}/chr16.dose.vcf.gz" \
    "${D}/chr17.dose.vcf.gz" \
    "${D}/chr18.dose.vcf.gz" \
    "${D}/chr19.dose.vcf.gz" \
    "${D}/chr20.dose.vcf.gz" \
    "${D}/chr21.dose.vcf.gz" \
    "${D}/chr22.dose.vcf.gz"
    
    #Count nvariants after merge
    bcftools index "${D}/1_AllChromosomes.vcf.gz"
    SNPs_1_AllChromosomes=$(bcftools index -n "${D}/1_AllChromosomes.vcf.gz")

echo "STEP 1 completed"
bcftools index "${D}/1_AllChromosomes.vcf.gz" 2>/dev/null || true
 SNPs_1_AllChromosomes=$(bcftools index -n "${D}/1_AllChromosomes.vcf.gz")
 
# ===================================
# STEP 2: INFO SCORE FILTER
# ===================================

echo ""
echo "======================================================="
echo "STEP 2: INFO SCORE FILTER"
echo "======================================================="

bcftools filter -Oz -i 'INFO/R2>0.9' \
    "${D}/1_AllChromosomes.vcf.gz" \
    > "${D}/2_Target.imputated.vcf.gz"
    
$runplink \
   --vcf "${D}/2_Target.imputated.vcf.gz" \
   --make-bed \
   --const-fid \
   --out "${D}/Target.imputated"

mv "${D}/Target.imputated.log" \
   "${L}/2_Target.imputated.log" \
   2>/dev/null || true



# Count variants after R2 filter
SNPs_2_AfterR2=$(wc -l < "${D}/Target.imputated.bim")
SNPs_2_RemovedR2=$((SNPs_1_AllChromosomes-SNPs_2_AfterR2))

echo "STEP 2 completed"

#===================================
# STEP 3A: PALINDROMIC SNP CHECK ONLY
# ===================================

echo ""
echo "======================================================="
echo "STEP 3A: PALINDROMIC SNP CHECK (NOT REMOVED)"
echo "======================================================="

awk '($5=="A" && $6=="T") || \
     ($5=="T" && $6=="A") || \
     ($5=="C" && $6=="G") || \
     ($5=="G" && $6=="C")' \
     "${D}/Target.imputated.bim" \
     | cut -f2 \
     > "${D}/3A_palindromic_snps_identification.txt"

SNPs_3A_Palindromic=$(wc -l < "${D}/3A_palindromic_snps_identification.txt")
SNPs_3A_After=$(wc -l < "${D}/Target.imputated.bim")

echo "Palindromic SNPs detected but NOT removed: ${SNPs_3A_Palindromic}"
echo "STEP 3A completed"

#===================================
# STEP 3B: Assign IDs only to variants with missing IDsset-missing-var-ids — First of any step, (only rename the ones with "." ID)
#===================================
$runplink --bfile "${D}/Target.imputated" \
    --set-missing-var-ids @:#:\$1:\$2 \
    --make-bed --out "${D}/Target.imputated"
    
#===================================
# STEP 3C: Remove duplicated rsIDs (1st round: before renaming)
#===================================
# Save information about duplicated rsIDs that will be removed
{
    echo "# Variants removed because their rsID was duplicated"
    echo -e "CHR\tSNP_ID\tCM\tPOS\tA1\tA2"
    awk 'NR==FNR {dup[$1]; next} ($2 in dup)' \
        "${D}/3C_IDduplicated.txt" \
        "${D}/Target.imputated.bim"
} > "${Q}/3C_removed_duplicated_rsIDs.txt"

$runplink --bfile "${D}/Target.imputated" \
    --write-snplist --out "${D}/3C_ID"
awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' \
    "${D}/3C_ID.snplist" "${D}/3C_ID.snplist" \
    | sort | uniq > "${D}/3C_IDduplicated.txt"

$runplink --bfile "${D}/Target.imputated" --exclude "${D}/3C_IDduplicated.txt" --make-bed --out "${D}/Target.imputated"    
    
SNPs_3C_Duplicated=$(wc -l < "${D}/3C_IDduplicated.txt")
SNPs_3C_After=$(wc -l < "${D}/Target.imputated.bim")

#===================================
# 3D: Remove duplicated positions — we keep those with higher frequency
#===================================
$runplink --bfile "${D}/Target.imputated" \
    --list-duplicate-vars \
    --freq \
    --out "${D}/3D_IDduplicated"

Rscript scripts/PRE0D_Duplicates.R "${D}"

$runplink --bfile "${D}/Target.imputated" \
    --exclude "${D}/3D_IDduplicated.exclude.txt" \
    --make-bed \
    --out "${D}/Target.imputated"
SNPs_3D_Duplicated=$(wc -l < "${D}/3D_IDduplicated.exclude.txt")
SNPs_3D_After=$(wc -l < "${D}/Target.imputated.bim")

# IMPORTANT.: We comment all these lines because we've already check that there's no SNP dropping if we rename, so we dont't have to do the mess of process of renaming and mapping 

##===================================
## Save original rsIDs (or ID) into a file to do mapping later
##===================================
#awk 'BEGIN{OFS="\t"} {print $2, $1":"$4":"$5":"$6}' \
#    "${D}/Target.imputated.bim" \
#    > "${D}/ID_map_original_to_new.txt"
#    
## Invertir el mapa 
#awk 'BEGIN{OFS="\t"} {print $2, $1}' \
#    "${D}/ID_map_original_to_new.txt" \
#    > "${D}/ID_map_new_to_rsID.txt"
    
##===================================
## 3E) Renombrar TODOS a chr:pos:ref:alt (PLINK 1.9 compatible)
##===================================
#$runplink --bfile "${D}/Target.imputated" \
#    --update-name "${D}/ID_map_original_to_new.txt" \
#    --make-bed \
#    --out "${D}/Target.imputated"
#SNPs_3E_After=$SNPs_3D_After

##===================================
## 3F) Remove duplicates once renamed (2n round dropping duplicates) 
##===================================
#$runplink --bfile "${D}/Target.imputated" \
#    --write-snplist \
#    --out "${D}/3F_ID"
#
#awk 'NR==FNR{a[$1]++;next}{if(a[$1]>1)print}' \
#    "${D}/3F_ID.snplist" "${D}/3F_ID.snplist" \
#    | sort | uniq \
#    > "${D}/3F_IDduplicated.txt"
#
#$runplink --bfile "${D}/Target.imputated" \
#    --exclude "${D}/3F_IDduplicated.txt" \
#    --make-bed \
#    --out "${D}/Target.imputated"
#SNPs_3F_Duplicated=$(wc -l < "${D}/3F_IDduplicated.txt")
#SNPs_3F_After=$(wc -l < "${D}/Target.imputated.bim")
echo "STEP 3 completed"

# ===================================
# STEP 4: MAF FILTER
# ===================================

echo ""
echo "======================================================="
echo "STEP 4: MAF FILTER"
echo "======================================================="

$runplink \
    --bfile "${D}/Target.imputated" \
    --freq \
    --out "${D}/3E_MAF"

mv "${D}/3E_MAF.log" "${L}/3E_MAF.log" 2>/dev/null || true

Rscript --no-save scripts/POS4_MAF_check.R "${D}/3E_MAF.frq"

mv MAF_distribution.pdf \
   "${Q}/3E_MAF_distribution.pdf" \
   2>/dev/null || true

$runplink \
    --bfile "${D}/Target.imputated" \
    --maf $s4mymaf \
    --make-bed \
    --out "${D}/Target.imputated_maf"

mv "${D}/Target.imputated_maf.log" \
   "${L}/4_Target.imputated_maf.log" \
   2>/dev/null || true

$runplink \
    --bfile "${D}/Target.imputated_maf" \
    --freq \
    --out "${D}/4_MAF_check"

mv "${D}/4_MAF_check.log" \
   "${L}/4_MAF_check.log" \
   2>/dev/null || true

Rscript --no-save scripts/POS4_MAF_check.R "${D}/4_MAF_check.frq"

mv MAF_distribution.pdf \
   "${Q}/4_MAF_distribution.pdf" \
   2>/dev/null || true

mv "${D}/3E_MAF.frq" "$Q/" 2>/dev/null || true
mv "${D}/4_MAF_check.frq" "$Q/" 2>/dev/null || true
SNPs_4_After=$(wc -l < "${D}/Target.imputated_maf.bim")
SNPs_4_RemovedMAF=$(( SNPs_3D_After - SNPs_4_After ))
echo "STEP 4 completed"

# ===================================
# STEP 5: MISSINGNESS REPORT ONLY
# ===================================

echo ""
echo "======================================================="
echo "STEP 5: MISSINGNESS REPORT ONLY - no SNPs removed"
echo "======================================================="

$runplink \
    --bfile "${D}/Target.imputated_maf" \
    --missing \
    --out "${D}/5_missing_report"

mv "${D}/5_missing_report.log" \
   "${L}/5_missing_report.log" \
   2>/dev/null || true

#Rscript --no-save scripts/POS5_Hist_miss.R \
#    "${D}/5_missing_report.imiss" \
#    "${D}/5_missing_report.lmiss"

mv "${D}/5_missing_report.imiss" "$Q/" 2>/dev/null || true
mv "${D}/5_missing_report.lmiss" "$Q/" 2>/dev/null || true
mv histimiss.pdf "${Q}/5_histimiss.pdf" 2>/dev/null || true
mv histlmiss.pdf "${Q}/5_histlmiss.pdf" 2>/dev/null || true

SNPs_5A_After=$SNPs_4_After
SNPs_5A_Removed=0
SNPs_5B_After=$SNPs_4_After
SNPs_5B_Removed=0

$runplink \
    --bfile "${D}/Target.imputated_maf" \
    --missing \
    --out "${D}/5_missing_report"

awk -v thr="$s5mymissingness1" 'NR==1 || $6 > thr' \
    "${D}/5_missing_report.imiss" \
    > "${Q}/5_individuals_fail_mind_${s5mymissingness1}.txt"

awk -v thr="$s5mymissingness2" 'NR==1 || $5 > thr' \
    "${D}/5_missing_report.lmiss" \
    > "${Q}/5_snps_fail_geno_${s5mymissingness2}.txt"

Individuals_5_Fail=$(($(wc -l < "${Q}/5_individuals_fail_mind_${s5mymissingness1}.txt") - 1))
SNPs_5_Fail=$(($(wc -l < "${Q}/5_snps_fail_geno_${s5mymissingness2}.txt") - 1))

mv "${D}/5_missing_report.log" "${L}/5_missing_report.log" 2>/dev/null || true
mv "${D}/5_missing_report.imiss" "$Q/" 2>/dev/null || true
mv "${D}/5_missing_report.lmiss" "$Q/" 2>/dev/null || true

SNPs_5A_After=$SNPs_4_After
SNPs_5A_Removed=0
SNPs_5B_After=$SNPs_4_After
SNPs_5B_Removed=0

echo "Individuals that would fail missingness > ${s5mymissingness1}: ${Individuals_5_Fail}"
echo "SNPs that would fail missingness > ${s5mymissingness2}: ${SNPs_5_Fail}"

{
echo "Individuals that would fail missingness > ${s5mymissingness1}: ${Individuals_5_Fail}"
echo "SNPs that would fail missingness > ${s5mymissingness2}: ${SNPs_5_Fail}"
echo "STEP 5 completed: missingness was reported but no SNPs or individuals were removed"
} | tee "${Q}/5_missingness_report.txt"

echo "STEP 5 completed: missingness was reported but no SNPs or individuals were removed"




echo "STEP 5 completed: missingness was reported but no SNPs were removed"


## =========================================
## Verify the rename to rsID  
## =========================================
#TOTAL_BIM=$(wc -l < "${D}/Target.imputated_maf.bim")
#TOTAL_MAP=$(wc -l < "${D}/ID_map_new_to_rsID.txt")

## Cuantos IDs del .bim actual tienen entrada en el mapa
#MATCHED=$(awk 'NR==FNR{map[$1]=1;next} $2 in map {c++} END{print c+0}' \
#    "${D}/ID_map_new_to_rsID.txt" \
#    "${D}/Target.imputated_maf.bim")
#
## Cuantos de esos tienen rsID real (no chr:pos format)
#WITH_RSID=$(awk 'NR==FNR{map[$2]=$1;next} \
#    ($2 in map) && (map[$2] !~ /^[0-9]+:/) {c++} END{print c+0}' \
#    "${D}/ID_map_new_to_rsID.txt" \
#    "${D}/Target.imputated_maf.bim")

#WITHOUT_RSID=$(( MATCHED - WITH_RSID ))
#UNMATCHED=$(( TOTAL_BIM - MATCHED ))

#echo "-----------------------------------------------"
#echo "  VERIFICACION: mapa rsID antes de restaurar"
#echo "-----------------------------------------------"
#echo "  SNPs en .bim actual         : $TOTAL_BIM"
#echo "  Entradas en mapa (ID_map)   : $TOTAL_MAP"
#echo "  SNPs con match en mapa      : $MATCHED"
#echo "  >> Con rsID real            : $WITH_RSID"
#echo "  >> Solo chr:pos (sin rsID)  : $WITHOUT_RSID"
#echo "  SNPs SIN match en mapa      : $UNMATCHED"
#if [ "$UNMATCHED" -gt 0 ]; then
#    echo "  WARNING: $UNMATCHED SNPs no encontrados en el mapa."
#    echo "  Posible reordenacion de alelos en pasos intermedios."
#    echo "  Ejemplos de IDs sin match:"
#    awk 'NR==FNR{map[$1]=1;next} !($2 in map) {print "   "$2}' \
#        "${D}/ID_map_new_to_rsID.txt" \
#        "${D}/Target.imputated_maf.bim" | head -5
#else
#    echo "  OK: todos los SNPs tienen match en el mapa."
#fi
#echo "-----------------------------------------------"
## Tornar al ID original rsID 
#$runplink \
#    --bfile "${D}/Target.imputated_maf" \
#    --update-name "${D}/ID_map_new_to_rsID.txt" \
#    --make-bed \
##    --out "${D}/Target.imputated_maf_rsID"


# ===================================
# FINAL OUTPUT
# ===================================

echo ""
echo "======================================================="
echo "FINALISING OUTPUT"
echo "======================================================="

mv "${D}/Target.imputated_maf.bed" \
   "${O}/${BATCH_NAME}_TOPMED_POSTimputed.bed"

mv "${D}/Target.imputated_maf.bim" \
   "${O}/${BATCH_NAME}_TOPMED_POSTimputed.bim"

mv "${D}/Target.imputated_maf.fam" \
   "${O}/${BATCH_NAME}_TOPMED_POSTimputed.fam"

rm -f "${D}"/*~ 2>/dev/null || true


# ===================================
# SUPER SUMMARY
# ===================================

SNPs_3A_Removed=0
SNPs_3C_Removed=$(( SNPs_3A_After - SNPs_3C_After ))
SNPs_3D_Removed=$(( SNPs_3C_After - SNPs_3D_After ))
SNPs_5_Removed=0
SNPs_Final=$SNPs_4_After
SNPs_Total_Removed=$(( SNPs_1_AllChromosomes - SNPs_Final ))

{
echo ""
echo "======================================================="
echo "   SUPER SUMMARY: ${BATCH_NAME}"
echo "======================================================="
printf "%-50s %12s %12s\n" "STEP" "SNPs kept" "SNPs removed"
echo "--------------------------------------------------------------------------"

printf "%-50s %12s %12s\n" "1.  AllChromosomes merged VCF" "$SNPs_1_AllChromosomes" "-"
printf "%-50s %12s %12s\n" "2.  After R2 > 0.9 filter" "$SNPs_2_AfterR2" "$SNPs_2_RemovedR2"
printf "%-50s %12s %12s\n" "3A. Palindromic SNPs checked only" "$SNPs_3A_After" "$SNPs_3A_Removed"
printf "%-50s %12s %12s\n" "3C. After duplicated rsID removal" "$SNPs_3C_After" "$SNPs_3C_Removed"
printf "%-50s %12s %12s\n" "3D. After duplicated position removal" "$SNPs_3D_After" "$SNPs_3D_Removed"
printf "%-50s %12s %12s\n" "4.  After MAF > ${s4mymaf} filter" "$SNPs_4_After" "$SNPs_4_RemovedMAF"
printf "%-50s %12s %12s\n" "5.  Missingness report only" "$SNPs_Final" "$SNPs_5_Removed"

echo "--------------------------------------------------------------------------"
printf "%-50s %12s %12s\n" "FINAL POST-imputed bed/bim/fam" "$SNPs_Final" "$SNPs_Total_Removed"
echo "======================================================="


} | tee "${Q}/${BATCH_NAME}_SNP_summary.txt"

echo "Summary saved to: ${Q}/${BATCH_NAME}_SNP_summary.txt"

# ===================================
# CLEANUP: remove temporary bed/bim/fam from QC folder
# ===================================

echo "Cleaning up temporary bed/bim/fam from QC folder..."

find "$Q" -maxdepth 1 \( -name "*.bed" -o -name "*.bim" -o -name "*.fam" \) -delete
# ===================================
# CLEANUP INPUT DIRECTORY
# ===================================

echo ""
echo "======================================================="
echo "CLEANUP: INPUT DIRECTORY"
echo "======================================================="

# Move statistics/ folder to info
if [ -d "${D}/statistics" ]; then
    mv "${D}/statistics" "${I}/${BATCH_NAME}_statistics"
    echo "Moved statistics/ to info"
fi

# If the output final files are correctly generated, then remove intermediate plink files 
if [[ -f "${O}/${BATCH_NAME}_TOPMED_POSTimputed.bed" && \
      -f "${O}/${BATCH_NAME}_TOPMED_POSTimputed.bim" && \
      -f "${O}/${BATCH_NAME}_TOPMED_POSTimputed.fam" ]]; then

    # Remove heavy intermediate PLINK files that should not remain in the input folder
    rm -f "${D}"/Target.imputated.*
    rm -f "${D}"/Target.imputated_maf.*
    rm -f "${D}"/*.nosex

    echo "Removed intermediate PLINK files from input directory"
else
    echo "WARNING: final output not found, intermediate PLINK files were kept"
fi


# Move .txt and .log files to info
find "$D" -maxdepth 1 -name "*.txt" -exec mv {} "$I/" \; 2>/dev/null || true
find "$D" -maxdepth 1 -name "*.log" -exec mv {} "$L/" \; 2>/dev/null || true

echo "Moved .txt and .log files to info"

# Move .vcf.gz and index files to trash
find "$D" -maxdepth 1 -name "*.vcf.gz"     -exec mv {} "$T/" \; 2>/dev/null || true
find "$D" -maxdepth 1 -name "*.vcf.gz.tbi" -exec mv {} "$T/" \; 2>/dev/null || true
find "$D" -maxdepth 1 -name "*.vcf.gz.csi" -exec mv {} "$T/" \; 2>/dev/null || true

echo "Moved VCF files to trash: ${T}"

# Delete trash contents
rm -rf "${T:?}"/*
echo "Trash cleared"

echo "Cleanup completed"

echo ""
echo "#######################################################"
echo "### COMPLETED: ${BATCH_NAME}"
echo "#######################################################"

}

###################################################################
######################## RUN ######################################
###################################################################

for BATCH in "${SELECTED_BATCHES[@]}"; do

    run_pipeline "$BATCH"

done

echo ""
echo "======================================================="
echo "ALL PIPELINES FINISHED"
echo "======================================================="