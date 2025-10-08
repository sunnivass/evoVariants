#!/bin/bash

set -euo pipefail

# Ensure required tools are available
if ! command -v bcftools >/dev/null 2>&1; then
    echo "[ERROR] bcftools is not available on PATH." >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "[ERROR] samtools is not available on PATH." >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 is not available on PATH." >&2
    exit 1
fi

# Parameters
REF_FASTA="ref/fasta/sacCer3_SFS01decoy_PfRpn11.masked2.fa"
GFF3="ref/gff/genomic_PfRPN11.gff"
SAMPLESHEET="data/samplesheet_PfRPN11.csv"
VCF_DIR="output/vcf"
OUT_DIR="output/annotated_vcf"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Track per-sample outputs for downstream aggregation
declare -a SAMPLE_NAMES=()
declare -A SAMPLE_TSV_PATHS=()

# Index reference if not already
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Indexing reference genome..."
    samtools faidx "$REF_FASTA"
fi

# Annotate with bcftools csq
# Read the CSV file, skip the header - modified to avoid subshell issues
while IFS=, read -r SAMPLE GROUP FASTQ1 FASTQ2; do
    IN_VCF="$VCF_DIR/${SAMPLE}.filtered.vcf"
    OUT_VCF="$OUT_DIR/${SAMPLE}.annotated.vcf.gz"
    OUT_TSV="$OUT_DIR/${SAMPLE}.annotated.tsv"

    echo "[INFO] Annotating variants for ${SAMPLE}..."
    #bcftools reheader -s <(echo "$SAMPLE") "$IN_VCF" | \
    bcftools csq \
        -f "$REF_FASTA" \
        -g "$GFF3" \
        -v \
        -Oz -o "$OUT_VCF" \
        "$IN_VCF"

    # Index the compressed VCF (required for merging)
    bcftools index "$OUT_VCF"

    SAMPLE_NAMES+=("$SAMPLE")
    SAMPLE_TSV_PATHS["$SAMPLE"]="$OUT_TSV"

    # Convert to TSV for individual sample (include header)
    printf 'CHROM\tPOS\tREF\tALT\tBCSQ\tDP\tAF\n' > "$OUT_TSV"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\t%INFO/DP\t%INFO/AF\n' "$OUT_VCF" >> "$OUT_TSV"

    echo "[INFO] Annotation completed: ${SAMPLE}"
done < <(tail -n +2 "$SAMPLESHEET")

if [ ${#SAMPLE_NAMES[@]} -eq 0 ]; then
    echo "[ERROR] No samples found in $SAMPLESHEET." >&2
    exit 1
fi

# Merge all annotated VCFs
shopt -s nullglob
ANNOTATED_VCFS=("$OUT_DIR"/*.annotated.vcf.gz)
shopt -u nullglob

if [ ${#ANNOTATED_VCFS[@]} -eq 0 ]; then
    echo "[ERROR] No annotated VCFs found in $OUT_DIR; aborting merge." >&2
    exit 1
fi

echo "[INFO] Merging annotated VCFs..."
bcftools merge \
    -O z \
    -o "$OUT_DIR/merged.vcf.gz" \
    "${ANNOTATED_VCFS[@]}"

# Index merged VCF
bcftools index "$OUT_DIR/merged.vcf.gz"

# Convert merged per-sample TSV files into a combined table with DP/AF per sample
MERGED_VARIANTS_TSV="$OUT_DIR/merged.variants.tsv"
MERGE_SAMPLE_INFO=""
for SAMPLE in "${SAMPLE_NAMES[@]}"; do
    TSV_PATH="${SAMPLE_TSV_PATHS[$SAMPLE]}"
    if [ -z "$MERGE_SAMPLE_INFO" ]; then
        MERGE_SAMPLE_INFO="${SAMPLE}=${TSV_PATH}"
    else
        MERGE_SAMPLE_INFO+=";${SAMPLE}=${TSV_PATH}"
    fi
done

export MERGE_SAMPLE_INFO
export MERGED_VARIANTS_TSV
