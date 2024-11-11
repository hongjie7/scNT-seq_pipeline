# scNT-seq_pipeline
## Dynast Pipeline Tutorial for Single-Cell Metabolically Labeled New RNA Tagging Sequencing (scNT-seq) Analysis

This tutorial outlines the process for using the **Dynast pipeline** to analyze scNT-seq data. This includes building a STAR index, aligning reads, generating consensus reads, and performing background correction.

## 1. Build the STAR Index

Start by building the STAR index using the reference genome and annotation files.

```bash
# Build the STAR index using Dynast
dynast ref -i STAR GRCm39.genome.fa.gz gencode.vM29.annotation.gtf.gz
```

### Input Files:
- `GRCm39.genome.fa.gz`: Reference genome (FASTA file).
- `gencode.vM29.annotation.gtf.gz`: Gene annotation file (GTF format).

## 2. Define Paths and Variables

Set the paths for the STAR index, GTF annotation, and fastq file directory.

```bash
STAR_index=/path/to/your/STAR_index
gtf=/path/to/your/gtf/gencode/gencode_mouse/release_M29/gencode.vM29.annotation.gtf
path=/path/to/your/fastq
```

## 3. Run Alignment, Consensus, and Counting Steps

Loop through the samples to align reads, generate consensus reads, and count with Dynast.

```bash
for sample in sample1 sample2 sample3 control-sample;
do
    mkdir -p consensus/${sample}
    mkdir -p align/${sample}
    mkdir -p count/${sample}

    # Align reads
    dynast align -t 16 --strand forward -i ${STAR_index} -o align/${sample} \
    --STAR-overrides '--soloCellFilter TopCells 5000' -x dropseq \
    ${path}/${sample}_R2_001.fastq.gz ${path}/${sample}_R1_001.fastq.gz

    # Generate consensus reads
    dynast consensus -t 16 --strand forward -g ${gtf} --barcode-tag CB --umi-tag UB \
    align/${sample}/Aligned.sortedByCoord.out.bam -o consensus/${sample} \
    --barcodes align/${sample}/Solo.out/Gene/filtered/barcodes.tsv

    # Count reads with conversion
    dynast count -t 16 --strand forward -g ${gtf} --barcode-tag CB --umi-tag UB \
    consensus/${sample}/consensus.bam -o count/${sample} --conversion TC \
    --barcodes align/${sample}/Solo.out/Gene/filtered/barcodes.tsv

done
```

### Step Details:
- **`dynast align`**: Aligns reads using STAR.
- **`dynast consensus`**: Generates consensus reads.
- **`dynast count`**: Counts reads with specified conversion (e.g., TC).

## 4. Background Correction use Control Sample (without 4sU labeling)

Create directories for control counting and estimation.

```bash
mkdir -p control_count
mkdir -p control_estimate

ctrlBam=/path/to/your/consensus/control-sample/consensus.bam

# Perform control counting
dynast count -t 16 --control --snp-threshold 0.5 -o control_count \
--barcode-tag CB --umi-tag UB --conversion TC -g ${gtf} ${ctrlBam} \
--barcodes align/control-sample/Solo.out/Gene/filtered/barcodes.tsv

# Estimate the background conversion rate
dynast estimate -t 16 --control -o control_estimate control_count
```

### Step Details:
- **`--control`**: Specifies control counting to generate background SNP data.
- **`dynast estimate`**: Calculates background conversion rates and outputs to `p_e.csv`.

## 5. Alpha Correction for Experimental Samples

Run alpha correction to account for background conversion rates.

```bash
for sample in sample1 sample2 sample3;
do
    mkdir -p estimate_count/${sample}
    mkdir -p bgcorrect_count/${sample}

    # Perform count correction using control SNPs
    dynast count -t 16 --snp-csv control_count/snps.csv -o bgcorrect_count/${sample} \
    consensus/${sample}/consensus.bam -g ${gtf} --barcode-tag CB --umi-tag UB --conversion TC \
    --barcodes align/${sample}/Solo.out/Gene/filtered/barcodes.tsv

    # Alpha correction step
    dynast estimate -t 16 --p-e control_estimate/p_e.csv -o estimate_count/${sample} \
    bgcorrect_count/${sample} --method alpha --reads total

done
```

### Explanation:
- **`--snp-csv`**: Uses the control SNP CSV for correction.
- **`--method alpha`**: Applies the alpha correction method to estimate conversion rates.

## 6. Workflow and Output
<img width="467" alt="image" src="https://github.com/user-attachments/assets/33b9d7ed-eee3-403b-98ff-4d56aa1b4f7d">
<img width="363" alt="image" src="https://github.com/user-attachments/assets/cdab1e67-7971-433c-a951-64adee8b71fb">




---

This tutorial provides a complete guide for running the Dynast pipeline from raw data to count matrix output with T to C substitution, enabling RNA velocity and RNA kinetics analysis and insights into single-cell gene expression regulatory dynamics.

#### Reference:
https://dynast-release.readthedocs.io/en/latest/pipeline_usage.html \
https://anndata.readthedocs.io/en/latest/
