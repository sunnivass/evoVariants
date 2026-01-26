#!/bin/bash
set -euo pipefail

# PARAMETERS (default values; overridable via CLI flags)
REF="ref/fasta/masked/R64-1-1_SFS01decoy_PfRPN11.masked.fa"
SAMPLESHEET="data/samplesheet_PfRPN11-2.csv"
RESULTS_DIR="output"
THREADS=18
AF_CUTOFF=0.01
MAX_JOBS=1 # Number of samples to process in parallel

declare -a RUNNING_PIDS=()

usage() {
    cat <<'EOF'
Usage: scripts/evoVariants_lofreq.sh [options]

Options:
  -s, --samplesheet FILE    Path to samplesheet CSV (default: data/samplesheet_PfRPN11-2.csv)
  -n, --sample NAME         Sample name for single-sample mode
  -b, --bam FILE            Deduplicated BAM for single-sample mode (default: output/bam/<sample>.dedup.bam)
  -r, --ref FILE            Reference FASTA path (default: ref/fasta/masked/R64-1-1_SFS01decoy_PfRPN11.masked.fa)
  -t, --threads N           Number of threads to use (default: 18)
  -a, --af-cutoff FLOAT     Allele frequency cutoff for bcftools filter (default: 0.01)
  -h, --help                Show this help and exit

Either provide a samplesheet (default) or use --sample for a single sample.
For samplesheet mode, the script expects input BAMs at output/bam/<sample>.dedup.bam.
EOF
}

SAMPLE_OVERRIDE=""
BAM_OVERRIDE=""
SAMPLESHEET_SET=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samplesheet)
            [[ $# -lt 2 ]] && { echo "ERROR: --samplesheet requires a value." >&2; exit 1; }
            SAMPLESHEET="$2"
            SAMPLESHEET_SET=1
            shift 2
            ;;
        -n|--sample)
            [[ $# -lt 2 ]] && { echo "ERROR: --sample requires a value." >&2; exit 1; }
            SAMPLE_OVERRIDE="$2"
            shift 2
            ;;
        -b|--bam)
            [[ $# -lt 2 ]] && { echo "ERROR: --bam requires a value." >&2; exit 1; }
            BAM_OVERRIDE="$2"
            shift 2
            ;;
        -r|--ref)
            [[ $# -lt 2 ]] && { echo "ERROR: --ref requires a value." >&2; exit 1; }
            REF="$2"
            shift 2
            ;;
        -t|--threads)
            [[ $# -lt 2 ]] && { echo "ERROR: --threads requires a value." >&2; exit 1; }
            THREADS="$2"
            shift 2
            ;;
        -a|--af-cutoff)
            [[ $# -lt 2 ]] && { echo "ERROR: --af-cutoff requires a value." >&2; exit 1; }
            AF_CUTOFF="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "ERROR: Unknown option '$1'." >&2
            usage >&2
            exit 1
            ;;
    esac
done

USE_SAMPLESHEET=1
if [[ -n "$SAMPLE_OVERRIDE" || -n "$BAM_OVERRIDE" ]]; then
    USE_SAMPLESHEET=0
    if [[ -z "$SAMPLE_OVERRIDE" ]]; then
        echo "ERROR: --sample is required for single-sample mode." >&2
        exit 1
    fi
    if (( SAMPLESHEET_SET )); then
        echo "ERROR: Cannot use --samplesheet with --sample/--bam." >&2
        exit 1
    fi
else
    USE_SAMPLESHEET=1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --threads expects a positive integer (got '$THREADS')." >&2
    exit 1
fi

if (( THREADS < 1 )); then
    echo "ERROR: --threads must be >= 1." >&2
    exit 1
fi

if (( USE_SAMPLESHEET )); then
    if [[ ! -f "$SAMPLESHEET" ]]; then
        echo "ERROR: Samplesheet '$SAMPLESHEET' not found." >&2
        exit 1
    fi
else
    if [[ -n "$BAM_OVERRIDE" && ! -f "$BAM_OVERRIDE" ]]; then
        echo "ERROR: BAM file '$BAM_OVERRIDE' not found." >&2
        exit 1
    fi
fi

echo "[INFO] Configuration:"
echo "[INFO]   Reference: $REF"
echo "[INFO]   Threads: $THREADS"
echo "[INFO]   AF cutoff: $AF_CUTOFF"
if (( USE_SAMPLESHEET )); then
    echo "[INFO]   Samplesheet: $SAMPLESHEET"
else
    echo "[INFO]   Single sample: $SAMPLE_OVERRIDE"
    if [[ -n "$BAM_OVERRIDE" ]]; then
        echo "[INFO]     BAM: $BAM_OVERRIDE"
    fi
fi

# ---- sanity checks ----
for exe in lofreq samtools bcftools; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: $exe not found in PATH"; exit 1; }
done

mkdir -p "$RESULTS_DIR"/{bam,vcf,logs}

# Index reference for SAMTOOLS/lofreq (.fai) -- REQUIRED by lofreq
if [ ! -f "${REF}.fai" ]; then
    echo "[INFO] Indexing reference for samtools (FAI)..."
    samtools faidx "$REF" 2> "${RESULTS_DIR}/logs/samtools_faidx.err"
fi

process_sample() {
    local SAMPLE="$1"
    local BAM_IN="$2"

    (
        echo "[INFO] Processing sample: $SAMPLE"

        if [[ -z "$BAM_IN" ]]; then
            BAM_IN="${RESULTS_DIR}/bam/${SAMPLE}.dedup.bam"
        fi

        if [[ ! -f "$BAM_IN" ]]; then
            echo "ERROR: BAM file for $SAMPLE not found: $BAM_IN" >&2
            exit 1
        fi

        local INDELQ_BAM="${RESULTS_DIR}/bam/${SAMPLE}.dedup.indelq.bam"
        if [[ -f "$INDELQ_BAM" ]]; then
            echo "[INFO] Removing existing indelqual BAM for ${SAMPLE}."
            rm -f "$INDELQ_BAM" "${INDELQ_BAM}.bai" "${INDELQ_BAM}.csi"
        fi

        lofreq indelqual --dindel -f "$REF" \
            -o "$INDELQ_BAM" \
            "$BAM_IN" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_indelqual.err"

        samtools index "$INDELQ_BAM" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_index_dedup_indelq.err"

        set +e
        lofreq call-parallel \
            --pp-threads "$THREADS" \
            -f "$REF" \
            -o "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
            "$INDELQ_BAM" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_call_parallel.err"
        rc=$?
        set -e

        if [ $rc -ne 0 ]; then
            echo "[WARN] lofreq call-parallel failed for ${SAMPLE} (see logs). Retrying with single-threaded lofreq call..."
            lofreq call \
                -f "$REF" \
                -o "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
                "$INDELQ_BAM" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_call_single.err"
        fi

        bcftools filter -i "AF>${AF_CUTOFF}" \
            "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
            -o "${RESULTS_DIR}/vcf/${SAMPLE}.filtered.vcf" 2> "${RESULTS_DIR}/logs/${SAMPLE}.bcftools_filter.err"

        rm -f "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf"

        echo "[INFO] Finished $SAMPLE."
    ) &

    local pid=$!
    RUNNING_PIDS+=("$pid")

    while (( ${#RUNNING_PIDS[@]} >= MAX_JOBS )); do
        wait "${RUNNING_PIDS[0]}"
        local rc=$?
        if (( rc != 0 )); then
            exit $rc
        fi
        RUNNING_PIDS=("${RUNNING_PIDS[@]:1}")
    done
}

process_samples_stream() {
    local SOURCE_LABEL="$1"
    if ! IFS= read -r _header; then
        echo "ERROR: ${SOURCE_LABEL} is empty or missing a header." >&2
        exit 1
    fi

    while IFS=, read -r SAMPLE _GROUP _FASTQ_1 _FASTQ_2 || [[ -n "$SAMPLE" ]]; do
        [[ -z "$SAMPLE" ]] && continue
        process_sample "$SAMPLE" ""
    done
}

if (( USE_SAMPLESHEET )); then
    process_samples_stream "$SAMPLESHEET" < "$SAMPLESHEET"
else
    process_sample "$SAMPLE_OVERRIDE" "$BAM_OVERRIDE"
fi

for pid in "${RUNNING_PIDS[@]}"; do
    wait "$pid"
    rc=$?
    if (( rc != 0 )); then
        exit $rc
    fi
done

echo "[INFO] All samples processed successfully."
