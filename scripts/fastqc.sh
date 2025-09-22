#!/bin/bash
set -euo pipefail

# PARAMETERS
RAW_DIR="data/raw"
TRIM_DIR="data/trimmed"
QC_DIR="output/qc"
THREADS=6

mkdir -p "$TRIM_DIR" "$QC_DIR"/raw "$QC_DIR"/trimmed

# Run FastQC on raw files
for fq in "$RAW_DIR"/*.fastq.gz; do
    fastqc -t "$THREADS" -o "$QC_DIR/raw" "$fq"
done

# Run MultiQC to summarize all FastQC
multiqc "$QC_DIR"/raw -o "$QC_DIR"/raw

# Run fastp for trimming
for fq1 in "$RAW_DIR"/*_1.fastq.gz; do
    fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
    base=$(basename "$fq1" _1.fastq.gz)
    fastp -i "$fq1" -I "$fq2" \
          -o "$TRIM_DIR/${base}_1.trimmed.fastq.gz" \
          -O "$TRIM_DIR/${base}_2.trimmed.fastq.gz" \
          -w "$THREADS" \
          --html "$QC_DIR/trimmed/${base}_fastp.html" 
done

# Run FastQC on trimmed files
for fq in "$TRIM_DIR"/*.trimmed.fastq.gz; do
    fastqc -t "$THREADS" -o "$QC_DIR/trimmed" "$fq"
done

# Run MultiQC to summarize all FastQC
multiqc "$QC_DIR"/trimmed -o "$QC_DIR"/trimmed

echo "[INFO] QC and trimming complete."