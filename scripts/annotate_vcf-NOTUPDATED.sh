#!/bin/bash

set -euo pipefail

# Parameters
REF_FASTA="ref/sacCer3"
GFF3="data/reference/PfRPN11-ref.gff3"
SAMPLESHEET="data/samplesheet_PfRPN11.csv"
VCF_DIR="results/vcf"
OUT_DIR="output/annotated_vcf"
# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# Index reference if not already
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Indexing reference genome..."
    samtools faidx $REF_FASTA
fi

# Annotate with bcftools csq
# Read the CSV file, skip the header
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r SAMPLE GROUP FASTQ1 FASTQ2
do
    bcftools csq \
    -f $REF_FASTA \
    -g $GFF3 \
    -Oz -o $OUT_DIR/${SAMPLE}.annotated.vcf \
    $VCF_DIR/${SAMPLE}.filtered.vcf

    bcftools view $OUT_DIR/${SAMPLE}.annotated.vcf |
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\n' > $OUT_DIR/${SAMPLE}.annotated.tsv

    echo "Annotation completed: ${SAMPLE}"
done