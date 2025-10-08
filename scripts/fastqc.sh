#!/bin/bash
set -euo pipefail

# PARAMETERS
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

RAW_DIR="${PROJECT_ROOT}/data/raw"
TRIM_DIR="${PROJECT_ROOT}/data/trimmed"
QC_DIR="${PROJECT_ROOT}/output/qc"
THREADS=18

shopt -s nullglob

for exe in fastqc fastp multiqc; do
    if ! command -v "$exe" >/dev/null 2>&1; then
        echo "ERROR: Required executable '$exe' not found in PATH" >&2
        exit 1
    fi
done

mkdir -p "$TRIM_DIR" "$QC_DIR"/raw "$QC_DIR"/trimmed

# Run FastQC on raw files
raw_fastqs=("$RAW_DIR"/*.fastq.gz)
if [ ${#raw_fastqs[@]} -eq 0 ]; then
    echo "ERROR: No FASTQ files matching *.fastq.gz found in $RAW_DIR" >&2
    exit 1
fi

for fq in "${raw_fastqs[@]}"; do
    fastqc -t "$THREADS" -o "$QC_DIR/raw" "$fq"
done

# Run MultiQC to summarize all FastQC
if ! multiqc "$QC_DIR"/raw -o "$QC_DIR"/raw; then
    echo "[WARN] MultiQC failed for raw reads; continuing without summary (check MultiQC installation)." >&2
fi

# Remove raw FastQC artifacts to keep only the MultiQC summary
find "$QC_DIR/raw" -maxdepth 1 -name '*_fastqc.html' -delete
find "$QC_DIR/raw" -maxdepth 1 -name '*_fastqc.zip' -delete

# Run fastp for trimming
for fq1 in "$RAW_DIR"/*_1.fastq.gz; do
    fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
    if [ ! -f "$fq2" ]; then
        echo "[WARN] Mate file for $fq1 not found (expected $fq2); skipping pair." >&2
        continue
    fi

    base=$(basename "$fq1" _1.fastq.gz)
    fastp -i "$fq1" -I "$fq2" \
          -o "$TRIM_DIR/${base}_1.trimmed.fastq.gz" \
          -O "$TRIM_DIR/${base}_2.trimmed.fastq.gz" \
          -w "$THREADS" \
          --html "$QC_DIR/trimmed/${base}_fastp.html"
done

# Run FastQC on trimmed files
trimmed_fastqs=("$TRIM_DIR"/*.trimmed.fastq.gz)
if [ ${#trimmed_fastqs[@]} -eq 0 ]; then
    echo "[WARN] No trimmed FASTQ files found to run FastQC on." >&2
else
    for fq in "${trimmed_fastqs[@]}"; do
        fastqc -t "$THREADS" -o "$QC_DIR/trimmed" "$fq"
    done
fi

# Run MultiQC to summarize all FastQC
if ! multiqc "$QC_DIR"/trimmed -o "$QC_DIR"/trimmed; then
    echo "[WARN] MultiQC failed for trimmed reads; continuing without summary (check MultiQC installation)." >&2
fi

# Remove trimmed FastQC artifacts to keep only the MultiQC summary
find "$QC_DIR/trimmed" -maxdepth 1 -name '*_fastqc.html' -delete
find "$QC_DIR/trimmed" -maxdepth 1 -name '*_fastqc.zip' -delete

echo "[INFO] QC and trimming complete."
