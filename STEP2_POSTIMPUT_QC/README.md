# STEP2 Post-Imputation QC

`STEP2_POSTIMPUT_QC` processes chromosome-level TOPMed/Michigan Imputation Server results and produces a quality-controlled PLINK dataset for each selected batch.

The recommended execution pair is:

- `submit_STEP2.sbatch`: SLURM submission configuration and batch selection.
- `STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh`: non-interactive post-imputation pipeline executed by SLURM.

`STEP2_POSTIMPUT_DATAPROCESSING.sh` is an interactive variant.

## Workflow

For every selected batch, the pipeline:

1. Extracts the password-protected chromosome archives.
2. Records the VCF header and chromosome dimensions.
3. Concatenates chromosomes 1–22.
4. Retains variants with `INFO/R2 > 0.9`.
5. Converts the filtered VCF to PLINK BED/BIM/FAM format.
6. Reports palindromic SNPs without removing them.
7. assigns IDs to variants with missing IDs.
8. Removes duplicated rsIDs.
9. Removes duplicated positions, retaining the variant selected by frequency.
10. Applies `MAF >= 0.001`.
11. Reports sample and SNP missingness without removing either.
12. Writes final PLINK files, QC reports, logs, and an SNP summary.

## Requirements

The scripts target a SLURM-based Linux/HPC environment and load:

```text
PLINK/1.9b_6.21-x86_64
BCFtools/1.21-GCC-12.3.0
R/4.3.2-gfbf-2023a
```

They also require Bash, `unzip`, GNU command-line utilities, and the R packages used by the helper scripts.

## Expected project structure

```text
STEP2_project/
├── submit_STEP2.sbatch
├── STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh
├── scripts/
│   ├── PRE0D_Duplicates.R
│   └── POS4_MAF_check.R
├── A_input_postimput_processing/
│   ├── batch_passwords.txt
│   └── GEPI-BIOPSI_input_postimput_processing/
│       ├── chr_1.zip
│       ├── chr_2.zip
│       ├── ...
│       └── chr_22.zip
├── info_postimput_processing/
├── output_postimput_processing/
└── trash/
```

Each input directory must follow this exact pattern:

```text
<BATCH_NAME>_input_postimput_processing
```

## Batch selection

The job version reads two environment variables:

```bash
BATCH_MODE="SELECTED"   # SELECTED or ALL
BATCH_LIST="BATCH_NAME1"
```

`submit_STEP2.sbatch` currently exports these values and therefore processes them. Multiple selected batches must be separated by spaces:

```bash
export BATCH_LIST="BATCH_NAME1 BATCH_NAME2 BATCH_NAMEZ"
```

To process all detected batch directories:

```bash
export BATCH_MODE="ALL"
export BATCH_LIST=""
```

## Password file

`A_input_postimput_processing/batch_passwords.txt` uses one `BATCH=password` entry per line:

```text
BATCH_NAME1=<password>
BATCH_NAME2=<password>
```


## Running the pipeline

Submit the job from the project directory:

```bash
sbatch submit_STEP2.sbatch
```

The submission script uses its own directory as the pipeline directory. It can also be overridden explicitly:

```bash
STEP2_DIR=/absolute/path/to/STEP2_project sbatch submit_STEP2.sbatch
```

See [USERGUIDE.md](USERGUIDE.md) for setup, validation, monitoring, and troubleshooting instructions.

## Outputs

For a batch named `GEPI-BIOPSI`, the final genotype files are:

```text
output_postimput_processing/
└── GEPI-BIOPSI_output_postimput_processing/
    ├── GEPI-BIOPSI_TOPMED_POSTimputed.bed
    ├── GEPI-BIOPSI_TOPMED_POSTimputed.bim
    └── GEPI-BIOPSI_TOPMED_POSTimputed.fam
```

QC information is stored under:

```text
info_postimput_processing/
└── GEPI-BIOPSI_info_postimput_processing/
    ├── GEPI-BIOPSI_step2_logs/
    ├── GEPI-BIOPSI_step2_QC/
    └── GEPI-BIOPSI_statistics/    # when supplied by the server
```

Important reports include:

- `<BATCH>_TOPMED_characteristics.txt`
- `<BATCH>_chr_dim_postimputed.txt`
- `3C_removed_duplicated_rsIDs.txt`
- `3E_MAF_distribution.pdf`
- `4_MAF_distribution.pdf`
- `5_missingness_report.txt`
- `<BATCH>_SNP_summary.txt`

SLURM writes scheduler output and errors to the filenames configured in `submit_STEP2.sbatch`.

## Filtering policy

| Check | Threshold | Action |
|---|---:|---|
| Imputation quality | `INFO/R2 > 0.9` | Variants failing the threshold are removed. |
| Palindromic variants | A/T or C/G | Reported but retained. |
| Duplicate rsIDs | Repeated ID | Removed. |
| Duplicate positions | Same PLINK duplicate position | Lower-priority entry selected by `PRE0D_Duplicates.R` is removed. |
| Minor allele frequency | `MAF >= 0.001` | Variants below the threshold are removed. |
| Individual missingness | `F_MISS > 0.2` | Reported only. |
| SNP missingness | `F_MISS > 0.02` | Reported only. |

## Cleanup behavior

Once all three final BED/BIM/FAM files exist, intermediate PLINK files are removed from the batch input directory. VCFs and their indexes are moved to the batch trash directory and then permanently deleted. Original ZIP archives remain in the input directory.

Because the VCF cleanup is destructive, archive any VCFs that must be retained before running the pipeline.

## Script versions

Use `STEP2_POSTIMPUT_DATAPROCESSING_JOB.sh` with SLURM. Compared with `STEP2_POSTIMPUT_DATAPROCESSING.sh`, it provides non-interactive batch selection and uses `unzip`, consistent with the current job submission file.
