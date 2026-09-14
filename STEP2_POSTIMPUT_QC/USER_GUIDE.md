# STEP2 Post-Imputation QC User Guide

## 1. Prepare the project directory

Copy the two recommended scripts into the same directory:

```text
submit_STEP2.sbatch
STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh
```

Create the required folders (if necessary):

```bash
mkdir -p scripts
mkdir -p A_input_postimput_processing
mkdir -p info_postimput_processing
mkdir -p output_postimput_processing
mkdir -p trash
```

Place these helper scripts in `scripts/`:

```text
PRE0D_Duplicates.R
POS4_MAF_check.R
```

`POS5_Hist_miss.R` is not currently required because its call is commented out. Missingness tables are still produced.

## 2. Prepare a batch

Create one input directory per batch. For YOUR BATCH :

```bash
mkdir -p A_input_postimput_processing/{BATCH_NAME}_input_postimput_processing
```

Copy all 22 Michigan/TOPMed result archives into it:

```text
chr_1.zip
chr_2.zip
...
chr_22.zip
```

The names must match exactly. After extraction, each archive must produce the corresponding file:

```text
chr1.dose.vcf.gz
chr2.dose.vcf.gz
...
chr22.dose.vcf.gz
```

## 3. Configure archive passwords

Create:

```text
A_input_postimput_processing/batch_passwords.txt
```

Add one entry per batch:

```text
{BATCH_NAME}=<password>
```

The text before `=` must exactly match the batch name used in the directory. 

## 4. Select the batches

Edit `submit_STEP2.sbatch`.

To process only your new/selected batch:

```bash
export BATCH_MODE="SELECTED"
export BATCH_LIST="{BATCH_NAME}"
```

To process several batches sequentially in the same SLURM job:

```bash
export BATCH_MODE="SELECTED"
export BATCH_LIST="BATCH_NAME1 BATCH_NAME2"
```

To process every detected batch:

```bash
export BATCH_MODE="ALL"
export BATCH_LIST=""
```

The pipeline skips names in `BATCH_LIST` that do not have a matching input directory and stops if no valid batch remains.

## 5. Review SLURM resources

The supplied submission file requests:

```text
Time:       15 hours
CPUs:       4
Memory:     24 GB
Partition:  mr-06
```

Adjust these directives to the cluster and dataset size:

```bash
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=mr-06
```

Although four CPUs are requested, most commands are not currently passed an explicit thread count.

## 6. Perform preflight checks

From the project directory, verify the files:

```bash
test -f submit_STEP2.sbatch
test -f STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh
test -f scripts/PRE0D_Duplicates.R
test -f scripts/POS4_MAF_check.R
test -f A_input_postimput_processing/batch_passwords.txt
```

Check that all chromosome archives exist for a batch:

```bash
for chr in {1..22}; do
    test -f "A_input_postimput_processing/GEPI-BIOPSI_input_postimput_processing/chr_${chr}.zip" \
        || echo "Missing chr_${chr}.zip"
done
```

The unzip loop only warns about missing archives, but chromosome concatenation requires all 22 extracted VCFs and will fail if any are absent.

## 7. Submit and monitor

Submit:

```bash
sbatch submit_STEP2.sbatch
```

Record the job ID printed by `sbatch`. Monitor it with the cluster's normal SLURM commands:

```bash
squeue -j <job_id>
sacct -j <job_id> --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

You can also change teh names of the output reports in the submit file in order to monitor the progress.: 
```text
STEP2_POSTIMPUT_{BATCH_NAME}.out
STEP2_POSTIMPUT_{BATCH_NAME}.err
```

The job uses `set -euo pipefail`, so an unhandled command error, unset variable, or failed pipeline normally stops execution.

## 8. Validate the results

A successful batch must contain all three final files:

```text
output_postimput_processing/<BATCH>_output_postimput_processing/
├── <BATCH>_TOPMED_POSTimputed.bed
├── <BATCH>_TOPMED_POSTimputed.bim
└── <BATCH>_TOPMED_POSTimputed.fam
```

Check that they are non-empty:

```bash
ls -lh output_postimput_processing/{BATCH_NAME}_output_postimput_processing/
wc -l output_postimput_processing/{BATCH_NAME}_output_postimput_processing/*.bim
wc -l output_postimput_processing/{BATCH_NAME}_output_postimput_processing/*.fam
```

Review the summary:

```text
info_postimput_processing/<BATCH>_info_postimput_processing/
└── <BATCH>_step2_QC/<BATCH>_SNP_summary.txt
```

Also inspect:

- The number of variants before and after `R2 > 0.9`.
- Duplicate-rsID and duplicate-position exclusions.
- MAF plots before and after filtering.
- Individuals with `F_MISS > 0.2`.
- SNPs with `F_MISS > 0.02`.
- PLINK and SLURM error logs.

Missingness failures are reports only; the corresponding samples and variants remain in the final files.

## 9. Understand cleanup

After successful final output generation, the pipeline:

- Removes intermediate `Target.imputated.*` and `Target.imputated_maf.*` files.
- Moves text files and logs into the information directories.
- Moves merged and chromosome VCF files and indexes to `trash/<BATCH>_trash/`.
- Permanently deletes the contents of that trash directory.

Copy any VCFs that must be preserved to an archive location before submission. The ZIP files are retained and can be used to recreate them.

## 10. Troubleshooting

### No batches detected

Confirm that input directories are directly inside `A_input_postimput_processing/` and end with:

```text
_input_postimput_processing
```

### No password found

Ensure that the password-file key exactly equals the derived batch name. Avoid spaces around the batch name. The parser preserves spaces in the password but removes carriage returns.

### `unzip` fails

Check the password, ZIP integrity, and availability of the `unzip` command. Test one archive manually in a temporary directory.

### Concatenation fails

Verify that every `chrN.dose.vcf.gz` file exists and that chromosome VCF headers and sample columns are compatible. Missing ZIPs are warnings during extraction but become fatal at concatenation.

### `INFO/R2` is undefined

Inspect the VCF header saved as `<BATCH>_TOPMED_characteristics.txt`. Confirm that the imputation-quality field is named `R2`; a different imputation-server output may use another field.

### R duplicate processing fails

Confirm that `scripts/PRE0D_Duplicates.R` exists, its packages are installed, and it writes:

```text
3D_IDduplicated.exclude.txt
```

inside the batch input directory.

### MAF plot is missing

Confirm that `scripts/POS4_MAF_check.R` accepts the `.frq` filepath as its first argument and writes `MAF_distribution.pdf` in the job's working directory.

### Final output exists but missingness failures are present

This is expected behavior. Step 5 only creates reports; it intentionally sets the numbers removed by missingness to zero.

## Recommended script

Use `STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh` for SLURM execution. Do not substitute `STEP2_POSTIMPUT_DATAPROCESSING.sh` without reviewing its interactive selection logic, its `7z` dependency, and the ordering of duplicate-rsID commands.
