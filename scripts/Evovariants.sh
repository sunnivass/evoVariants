#!/bin/bash
set -euo pipefail

# PARAMETERS
REF="ref/fasta/sacCer3_SFS01decoy_PfRpn11.fa"
SAMPLESHEET="data/samplesheet_PfRPN11.csv"
RESULTS_DIR="output"
TRIM_DIR="data/trimmed"
READ_SOURCE="${READ_SOURCE:-trimmed}" # raw|trimmed (defaults to trimmed output from fastp)
THREADS=18
AF_CUTOFF=0.01
MAX_JOBS=1 # Number of samples to process in parallel

# ---- sanity checks ----
for exe in bwa-mem2 samtools lofreq bcftools; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: $exe not found in PATH"; exit 1; }
done

mkdir -p "$RESULTS_DIR"/{bam,vcf,logs}

# Index reference for BWA if not already done
if [ ! -f "${REF}.bwt" ]; then
    echo "[INFO] Indexing reference for BWA..."
    bwa-mem2 index "$REF" 2> "${RESULTS_DIR}/logs/bwa_index.err"
fi

# Index reference for SAMTOOLS/lofreq (.fai) -- REQUIRED by lofreq
if [ ! -f "${REF}.fai" ]; then
    echo "[INFO] Indexing reference for samtools (FAI)..."
    samtools faidx "$REF" 2> "${RESULTS_DIR}/logs/samtools_faidx.err"
fi

job_count=0

# Read the samplesheet directly (skip header without spawning tail)
{
    if ! IFS= read -r _header; then
        echo "ERROR: Samplesheet $SAMPLESHEET is empty or missing a header." >&2
        exit 1
    fi

    while IFS=, read -r SAMPLE GROUP FASTQ_1 FASTQ_2
    do
        (
        echo "[INFO] Processing sample: $SAMPLE (reads: $READ_SOURCE)"

        RAW_FASTQ_1="$FASTQ_1"
        RAW_FASTQ_2="$FASTQ_2"

        case "$READ_SOURCE" in
            trimmed)
                if [[ "$FASTQ_1" == *.trimmed.fastq.gz ]]; then
                    # Samplesheet already points to trimmed reads
                    :
                else
                    base1=$(basename "$FASTQ_1")
                    base2=$(basename "$FASTQ_2")

                    if [[ "$base1" != *_1.fastq.gz || "$base2" != *_2.fastq.gz ]]; then
                        echo "ERROR: Unexpected FASTQ naming for $SAMPLE; expected *_1.fastq.gz and *_2.fastq.gz." >&2
                        exit 1
                    fi

                    prefix=${base1%_1.fastq.gz}
                    expected_mate="${prefix}_2.fastq.gz"

                    if [[ "$base2" != "$expected_mate" ]]; then
                        echo "ERROR: FASTQ pair names for $SAMPLE do not match (${base1} / ${base2})." >&2
                        exit 1
                    fi

                    FASTQ_1="${TRIM_DIR}/${prefix}_1.trimmed.fastq.gz"
                    FASTQ_2="${TRIM_DIR}/${prefix}_2.trimmed.fastq.gz"
                fi
                ;;
            raw)
                :
                ;;
            *)
                echo "ERROR: Unsupported READ_SOURCE '$READ_SOURCE'. Use 'raw' or 'trimmed'." >&2
                exit 1
                ;;
        esac

        # Check that files exist
        if [[ ! -f "$FASTQ_1" || ! -f "$FASTQ_2" ]]; then
            echo "ERROR: FASTQ files for $SAMPLE not found!" >&2
            if [[ "$READ_SOURCE" == "trimmed" ]]; then
                echo "       Expected trimmed files at: $FASTQ_1 and $FASTQ_2" >&2
                echo "       Run scripts/fastqc.sh before scripts/evovariants.sh." >&2
                echo "       Raw files listed in samplesheet: $RAW_FASTQ_1 and $RAW_FASTQ_2" >&2
            fi
            exit 1
        fi
        
        # Alignment -> sorted BAM (skip intermediate BAM conversion)
        bwa-mem2 mem -t "$THREADS" "$REF" "$FASTQ_1" "$FASTQ_2" 2> "${RESULTS_DIR}/logs/${SAMPLE}.bwa2.err" \
          | samtools sort -@ "$THREADS" -O BAM -o "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_sort.err"

        # Add indel qualities (required by LoFreq)
        lofreq indelqual --dindel -f "$REF" \
            -o "${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_indelqual.err"

        samtools index "${RESULTS_DIR}/bam/${SAMPLE}.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_index_indelq.err"

        # Delete sorted BAM to save space (keep indelq only)
        rm -f "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam"

        # Variant calling
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
        ) &
    job_count=$((job_count+1))
    if (( job_count % MAX_JOBS == 0 )); then
        wait
    fi
    done
} < "$SAMPLESHEET"
wait

echo "[INFO] All samples processed successfully."
