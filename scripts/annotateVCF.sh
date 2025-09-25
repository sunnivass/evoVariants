#!/bin/bash

set -euo pipefail

# Parameters
REF_FASTA="ref/fasta/sacCer3_SFS01decoy_PfRpn11.masked.fa"
GFF3="ref/gff/genomic_PfRPN11.gff"
SAMPLESHEET="data/samplesheet_PfRPN11.csv"
VCF_DIR="output/vcf"
OUT_DIR="output/annotated_vcf"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Index reference if not already
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Indexing reference genome..."
    samtools faidx "$REF_FASTA"
fi

# Annotate with bcftools csq
# Read the CSV file, skip the header
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r SAMPLE GROUP FASTQ1 FASTQ2
do
    OUT_VCF="$OUT_DIR/${SAMPLE}.annotated.vcf.gz"
    OUT_TSV="$OUT_DIR/${SAMPLE}.annotated.tsv"

    echo "[INFO] Annotating variants for ${SAMPLE}..."
    bcftools csq \
        -f "$REF_FASTA" \
        -g "$GFF3" \
        -Oz -o "$OUT_VCF" \
        "$VCF_DIR/${SAMPLE}.filtered.vcf"

    # Index the compressed VCF (required for merging)
    bcftools index "$OUT_VCF"

    # Convert to TSV for individual sample
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\n' "$OUT_VCF" > "$OUT_TSV"

    echo "[INFO] Annotation completed: ${SAMPLE}"
done

# Merge all annotated VCFs
echo "[INFO] Merging annotated VCFs..."
bcftools merge \
    -O z \
    -o "$OUT_DIR/merged.vcf.gz" \
    "$OUT_DIR"/*.annotated.vcf.gz

# Index merged VCF
bcftools index "$OUT_DIR/merged.vcf.gz"

# Convert merged VCF to TSV with sample information
echo "[INFO] Converting merged VCF to TSV..."
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ[\t%SAMPLE=%GT]\n' \
    "$OUT_DIR/merged.vcf.gz" > "$OUT_DIR/merged.variants.tsv"

echo "[INFO] All done! Merged variants saved to: $OUT_DIR/merged.variants.tsv"