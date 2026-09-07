# STEP1_PREIMPUT User Guide

This guide explains how to prepare inputs, run the pipeline, and find the results.

Use this document when you want to process one or more batches. Use `README.md` for the internal workflow and script-level details.

## 1. Before You Start

Work from the `STEP1_PREIMPUT` directory:

```bash
cd /path/to/STEP1_PREIMPUT
```

Expected top-level structure:

```text
STEP1_PREIMPUT/
|-- STEP1_PREIMPUT_QC_DP.sh
|-- STEP1_PREIMPUT_QC_DP_JOB.sh
|-- submit_STEP1_PREIMPUT.sbatch
|-- main_scripts/
|-- scripts_1A_QC/
|-- scripts_1B_DP/
|-- scripts_1C_PCA/
|-- 1000G_ref_hg19/
|-- A_input_geno_rawdata/
|-- B_output_preimput/
`-- info_preimput/
```

The user normally prepares:

```text
A_input_geno_rawdata/
1000G_ref_hg19/      # only needed if PCA will be run
```

`B_output_preimput/` and `info_preimput/` are created or filled by the pipeline.

## 2. Prepare One Batch

Create a folder using this exact pattern:

```text
A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/
```

Inside it, the PLINK prefix must also use the batch name:

```text
A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_geno_rawdata.bed
A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_geno_rawdata.bim
A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_geno_rawdata.fam
```

Example:

```text
A_input_geno_rawdata/GOBIS_GSA25_geno_rawdata/GOBIS_GSA25_geno_rawdata.bed
A_input_geno_rawdata/GOBIS_GSA25_geno_rawdata/GOBIS_GSA25_geno_rawdata.bim
A_input_geno_rawdata/GOBIS_GSA25_geno_rawdata/GOBIS_GSA25_geno_rawdata.fam
```

PED/MAP is also accepted:

```text
{BATCH_NAME}_geno_rawdata.ped
{BATCH_NAME}_geno_rawdata.map
```

## 3. Optional: Keep Only Selected IIDs

If you only want to process some individuals from a batch, create:

```text
A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_IIDs.txt
```

Format:

```text
IID_1
IID_2
IID_3
```

One IID per line. These values are matched against column 2 of the `.fam` file.

If this file does not exist, the pipeline processes the full batch.

## 4. Check Sex Information

Before running QC, check whether sex is coded in the `.fam` file. Column 5 is sex:

```bash
awk '{print $5}' A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_geno_rawdata.fam | sort | uniq -c
```

PLINK sex coding:

```text
1 = male
2 = female
0 = unknown
```

If sex is missing or needs to be corrected, add:

```text
A_input_geno_rawdata/{BATCH_NAME}_geno_rawdata/{BATCH_NAME}_pheno_sex.txt
```

Format:

```text
FID IID SEX
```

Example:

```text
F001 ID001 2
F002 ID002 1
```

The pipeline uses this file with PLINK `--update-sex` before the sex discrepancy check.

## 5. Prepare 1000G Reference For PCA

This is only required if you will run `STEP 1C_PCA`.

Expected directory:

```text
1000G_ref_hg19/
```

Preferred files:

```text
1000G_ref_hg19/1000G_QC_PRC_v4.bed
1000G_ref_hg19/1000G_QC_PRC_v4.bim
1000G_ref_hg19/1000G_QC_PRC_v4.fam
1000G_ref_hg19/integrated_call_samples_v3.20200731.ALL.ped
```

The population metadata file must contain sample IDs and 1000G population labels.

## 6. Choose What To Run

`STEP 1A_QC` always runs.

Then you choose:

```text
STEP 1B_DP  -> create VCFs for Michigan Imputation Server
STEP 1C_PCA -> run PCA against 1000G reference
```

Common choices:

```text
Only QC:
STEP 1B_DP = NO
STEP 1C_PCA = NO

QC + PCA:
STEP 1B_DP = NO
STEP 1C_PCA = YES

QC + imputation VCF preparation:
STEP 1B_DP = YES
STEP 1C_PCA = NO

Full pre-imputation run:
STEP 1B_DP = YES
STEP 1C_PCA = YES
```

## 7. Run Interactively

From the `STEP1_PREIMPUT` directory:

```bash
bash STEP1_PREIMPUT_QC_DP.sh
```

The script will:

1. Detect available batches.
2. Ask whether to run all batches or selected batches.
3. Ask whether to run `STEP 1B_DP`.
4. Ask whether to run `STEP 1C_PCA`.

## 8. Run As A SLURM Job

Edit:

```text
submit_STEP1_PREIMPUT.sbatch
```

Set:

```bash
export BATCH_MODE="SELECTED"
export BATCH_LIST="GOBIS_GSA25"
export RUN_STEP1B_DP="NO"
export RUN_STEP1C_PCA="YES"
```

For several selected batches:

```bash
export BATCH_LIST="GOBIS_GSA25 cegen2021 cegen2022"
```

For all detected batches:

```bash
export BATCH_MODE="ALL"
```

Submit:

```bash
sbatch submit_STEP1_PREIMPUT.sbatch
```

Monitor:

```bash
squeue -u $USER
tail -f STEP1_PREIMPUT_<JOBID>.out
tail -f STEP1_PREIMPUT_<JOBID>.err
```

## 9. Important Outputs

For each batch, outputs are created under:

```text
B_output_preimput/{BATCH}_output_preimput/
info_preimput/{BATCH}_info_preimput/
```

### STEP 1A_QC Final PLINK Dataset

```text
B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.bed
B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.bim
B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.fam
```

This is the main QC-cleaned pre-imputation target dataset.

### STEP 1A_QC Summary

```text
info_preimput/{BATCH}_info_preimput/_info_STEP1_QC/{BATCH}_STEP1_QC_summary.txt
```

### Relatedness Outputs

```text
info_preimput/{BATCH}_info_preimput/_info_STEP1_QC/6B_pihat_min0.125.genome
info_preimput/{BATCH}_info_preimput/_info_STEP1_QC/6B_related_individuals_not_removed.txt
```

Related individuals are identified but not removed from the final `1A_output_QC` dataset.

### STEP 1B_DP Final VCFs

```text
B_output_preimput/{BATCH}_output_preimput/1B_output_DP/QCtargetdata-updated-chr1.vcf.gz
B_output_preimput/{BATCH}_output_preimput/1B_output_DP/QCtargetdata-updated-chr1.vcf.gz.tbi
...
B_output_preimput/{BATCH}_output_preimput/1B_output_DP/QCtargetdata-updated-chr22.vcf.gz
B_output_preimput/{BATCH}_output_preimput/1B_output_DP/QCtargetdata-updated-chr22.vcf.gz.tbi
```

These files are the ones to upload to Michigan Imputation Server.

### STEP 1C_PCA Outputs

```text
B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/{BATCH}_STEP3_PCA_summary.txt
B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/{BATCH}_PCA_results.covariate
B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/{BATCH}_PCA_3_target_ancestry_ALL_METHODS_details.txt
B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/{BATCH}_PCA_3_ancestry_method_summary.txt
B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/PCA_results/3_PCA_INTEGRATED_ANCESTRY_REPORT.pdf
```

The integrated PCA report includes:

```text
EUR +/- 3SD
Nearest 1000G centroid
EUR PC1-PC2 ellipse
Mahalanobis
Consensus across methods
```

## 10. Check That A Run Finished Correctly

Check the final terminal or job output:

```text
### STEP 1A_QC completed
### STEP 1B_DP completed/skipped
### STEP 1C_PCA completed/skipped
```

Check final files:

```bash
ls B_output_preimput/{BATCH}_output_preimput/1A_output_QC/
ls info_preimput/{BATCH}_info_preimput/_info_STEP1_QC/
```

If `STEP 1B_DP=YES`:

```bash
ls B_output_preimput/{BATCH}_output_preimput/1B_output_DP/*.vcf.gz
```

If `STEP 1C_PCA=YES`:

```bash
ls B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/
ls B_output_preimput/{BATCH}_output_preimput/1C_output_PCA/PCA_results/
```

## 11. Cleanup After Checking Results

After confirming that QC completed correctly, the final `STEP 1A_QC` dataset can be removed only if it will not be needed again for `STEP 1B_DP` or `STEP 1C_PCA`.

Files:

```text
B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.bed
B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.bim
B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.fam
```

Only delete these after:

1. The QC summary has been reviewed.
2. Any desired `STEP 1B_DP` or `STEP 1C_PCA` run has already completed.
3. The final downstream outputs have been checked.

Example:

```bash
rm -f B_output_preimput/{BATCH}_output_preimput/1A_output_QC/{BATCH}_preimput_QC.{bed,bim,fam}
```

## 12. Common Problems

### Windows line endings

If you see errors like `$'\r': command not found` or `set: pipefail: invalid option`, run:

```bash
sed -i 's/\r$//' STEP1_PREIMPUT_QC_DP.sh STEP1_PREIMPUT_QC_DP_JOB.sh submit_STEP1_PREIMPUT.sbatch main_scripts/*.sh scripts_1C_PCA/*.R
```

### Missing R script

If you see:

```text
ERROR: Missing required R script
```

Check that the expected script exists in the correct folder:

```bash
ls scripts_1A_QC/
ls scripts_1B_DP/
ls scripts_1C_PCA/
```

### Rscript or bcftools not found

Load the required modules and check paths:

```bash
module load modulepath/UPF/apps
module load modulepath/noarch
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
module load BCFtools/1.21-GCC-12.3.0

which plink
which Rscript
which bcftools
```

### No batches detected

Check the directory names:

```bash
find A_input_geno_rawdata -maxdepth 2 -type f | head
```

The folder must end with:

```text
_geno_rawdata
```

and the PLINK prefix must also be:

```text
{BATCH_NAME}_geno_rawdata
```

