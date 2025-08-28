#!/bin/bash
set -euo pipefail

# PARAMETERS
REF="ref/SFS01decoy_PfRPN11.fa"
SAMPLESHEET="data/test-samplesheet_PfRPN11.csv"
RESULTS_DIR="output"
THREADS=8
AF_CUTOFF=0.01

# ---- sanity checks ----
for exe in bwa samtools lofreq bcftools; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: $exe not found in PATH"; exit 1; }
done

mkdir -p "$RESULTS_DIR"/{bam,vcf,logs}

# Index reference for BWA if not already done
if [ ! -f "${REF}.bwt" ]; then
    echo "[INFO] Indexing reference for BWA..."
    bwa index "$REF" 2> "${RESULTS_DIR}/logs/bwa_index.err"
fi

# Index reference for SAMTOOLS/lofreq (.fai) -- REQUIRED by lofreq
if [ ! -f "${REF}.fai" ]; then
    echo "[INFO] Indexing reference for samtools (FAI)..."
    samtools faidx "$REF" 2> "${RESULTS_DIR}/logs/samtools_faidx.err"
fi

# Read the CSV file, skip the header
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r SAMPLE GROUP FASTQ_1 FASTQ_2
do
    echo "[INFO] Processing sample: $SAMPLE"

    # Check that files exist
    if [[ ! -f "$FASTQ_1" || ! -f "$FASTQ_2" ]]; then
        echo "ERROR: FASTQ files for $SAMPLE not found!"
        exit 1
    fi
    
    # Alignment -> sorted BAM
    bwa mem -t "$THREADS" "$REF" "$FASTQ_1" "$FASTQ_2" 2> "${RESULTS_DIR}/logs/${SAMPLE}.bwa.err" \
      | samtools view -Sb - 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_view.err" \
      | samtools sort -@ "$THREADS" -o "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_sort.err"
    
    samtools index "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_index_sorted.err"

    # Add indel qualities (required by LoFreq)
    lofreq indelqual --dindel -f "$REF" \
        -o "${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam" \
        "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_indelqual.err"

    samtools index "${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_index_indelq.err"

    # Delete sorted BAM to save space (keep indelq only)
    rm -f "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" \
          "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam.bai"

    # Variant calling
    # Try parallel first; if it fails, fall back to single-process to surface errors better.
    set +e
    lofreq call-parallel \
        --pp-threads "$THREADS" \
        -f "$REF" \
        -o "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
        "${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_call_parallel.err"
    rc=$?
    set -e

    if [ $rc -ne 0 ]; then
        echo "[WARN] lofreq call-parallel failed for ${SAMPLE} (see logs). Retrying with single-threaded lofreq call..."
        lofreq call \
            -f "$REF" \
            -o "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
            "${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_call_single.err"
    fi

    # Filtering by allele frequency
    bcftools filter -i "AF>${AF_CUTOFF}" \
        "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
        -o "${RESULTS_DIR}/vcf/${SAMPLE}.filtered.vcf" 2> "${RESULTS_DIR}/logs/${SAMPLE}.bcftools_filter.err"
    
    # Delete raw VCF after filtering
    rm -f "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf"

    echo "[INFO] Finished $SAMPLE."
done

echo "[INFO] All samples processed successfully."
