#!/bin/bash

set -euo pipefail

# PARAMETERS
REF="ref/SFS01decoy_PfRPN11.fa"
SAMPLESHEET="data/samplesheet_PfRPN11.csv"
RESULTS_DIR="output"
THREADS=8
AF_CUTOFF=0.01

mkdir -p $RESULTS_DIR/bam
mkdir -p $RESULTS_DIR/vcf

# Index reference genome if not already done
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index $REF
fi

# Read the CSV file, skip the header
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r SAMPLE GROUP FASTQ_1 FASTQ_2
do
    echo "Processing sample: $SAMPLE"

    # Check that files exist
    if [[ ! -f "$FASTQ_1" || ! -f "$FASTQ_2" ]]; then
        echo "Error: FASTQ files for $SAMPLE not found!"
        exit 1
    fi
    
    # Alignment
    bwa mem -t $THREADS $REF "$FASTQ_1" "$FASTQ_2" \
        | samtools view -Sb - \
        | samtools sort -@ $THREADS -o ${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam
    
    samtools index ${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam

    # Add indel qualities (required by LoFreq)
    lofreq indelqual --dindel -f $REF \
        -o ${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam \
        ${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam

    samtools index ${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam

    # Variant calling
    lofreq call-parallel\
        --pp-threads $THREADS \
        -f $REF \
        -o ${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf \
        ${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam

    # Filtering by allele frequency
    bcftools filter -i "AF>${AF_CUTOFF}" \
        ${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf \
        -o ${RESULTS_DIR}/vcf/${SAMPLE}.filtered.vcf

    echo "Finished $SAMPLE."
done

echo "All samples processed successfully."
