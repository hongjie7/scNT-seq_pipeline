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
STAR_index=/mnt/data1/hongjie/reference/dropseq/gencode/gencode_mouse/release_M29/STAR
gtf=/mnt/data1/hongjie/reference/dropseq/gencode/gencode_mouse/release_M29/gencode.vM29.annotation.gtf
path=/mnt/data1/hongjie/project/scNTseq/E16/fastq
```

## 3. Run Alignment, Consensus, and Counting Steps

Loop through the samples to align reads, generate consensus reads, and count with Dynast.

```bash
for sample in E16Ctx1-1_S1 E16Ctx1-2_S2 E16Ctx2-1_S3 E16Ctx2-2_S4 E16Ctx3-saline_S5; do
    mkdir -p consensus/${sample}
    mkdir -p align/${sample}
    mkdir -p unCorrect_count/${sample}

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
    consensus/${sample}/consensus.bam -o unCorrect_count/${sample} --conversion TC \
    --barcodes align/${sample}/Solo.out/Gene/filtered/barcodes.tsv

done
```

### Step Details:
- **`dynast align`**: Aligns reads using STAR.
- **`dynast consensus`**: Generates consensus reads.
- **`dynast count`**: Counts reads with specified conversion (e.g., TC).

## 4. Background Correction and Control Analysis

Create directories for control counting and estimation.

```bash
mkdir -p control_count
mkdir -p control_estimate

ctrlBam=/mnt/data2/hongjie/project/invivoscNTseq/data/dynast_output/E16Ctx/consensus/E16Ctx3-saline_S5/consensus.bam

# Perform control counting
dynast count -t 16 --control --snp-threshold 0.5 -o control_count \
--barcode-tag CB --umi-tag UB --conversion TC -g ${gtf} ${ctrlBam} \
--barcodes align/E16Ctx3-saline_S5/Solo.out/Gene/filtered/barcodes.tsv

# Estimate the background conversion rate
dynast estimate -t 16 --control -o control_estimate control_count
```

### Step Details:
- **`--control`**: Specifies control counting to generate background SNP data.
- **`dynast estimate`**: Calculates background conversion rates and outputs to `p_e.csv`.

## 5. Alpha Correction for Experimental Samples

Run alpha correction to account for background conversion rates.

```bash
for sample in E16Ctx1-1_S1 E16Ctx1-2_S2 E16Ctx2-1_S3 E16Ctx2-2_S4; do
    mkdir -p correct_estimate_alpha/${sample}

    # Perform count correction using control SNPs
    dynast count -t 16 --snp-csv control_count/snps.csv -o correct_count/${sample} \
    consensus/${sample}/consensus.bam -g ${gtf} --barcode-tag CB --umi-tag UB --conversion TC \
    --barcodes align/${sample}/Solo.out/Gene/filtered/barcodes.tsv

    # Alpha correction step
    dynast estimate -t 16 --p-e control_estimate/p_e.csv -o correct_estimate_alpha/${sample} \
    correct_count/${sample} --method alpha --reads total

done
```

### Explanation:
- **`--snp-csv`**: Uses the control SNP CSV for correction.
- **`--method alpha`**: Applies the alpha correction method to estimate conversion rates.

---

This tutorial provides a complete guide for running the Dynast pipeline from raw data to corrected output, enabling RNA velocity analysis and insights into single-cell gene expression dynamics.

