#!/bin/bash
set -euo pipefail

# PARAMETERS (default values; overridable via CLI flags)
REF="ref/fasta/sacCer3_SFS01decoy_HsGTaseCeg1.masked2.fa"
SAMPLESHEET="data/samplesheet_HsCeg1.csv"
RESULTS_DIR="output"
TRIM_DIR="data/trimmed"
READ_SOURCE="${READ_SOURCE:-trimmed}" # raw|trimmed (defaults to trimmed output from fastp)
THREADS=18
AF_CUTOFF=0.01
MAX_JOBS=1 # Number of samples to process in parallel

declare -a RUNNING_PIDS=()

usage() {
    cat <<'EOF'
Usage: scripts/Evovariants.sh [options]

Options:
  -s, --samplesheet FILE    Path to samplesheet CSV (default: data/samplesheet_PfRPN11.csv)
  -1, --fastq1 FILE         FASTQ R1 path for single-sample mode (requires --fastq2 and --sample)
  -2, --fastq2 FILE         FASTQ R2 path for single-sample mode (requires --fastq1 and --sample)
  -n, --sample NAME         Sample name for single-sample mode
  -r, --ref FILE            Reference FASTA path (default: ref/fasta/sacCer3_SFS01decoy_PfRpn11.fa)
  -t, --threads N           Number of threads to use (default: 18)
  -a, --af-cutoff FLOAT     Allele frequency cutoff for bcftools filter (default: 0.01)
  -R, --read-source MODE    Read source: raw|trimmed (default: trimmed)
  -h, --help                Show this help and exit

Either provide a samplesheet (default) or use --fastq1/--fastq2/--sample for a single sample.
EOF
}

FASTQ_OVERRIDE_1=""
FASTQ_OVERRIDE_2=""
SAMPLE_OVERRIDE=""
SAMPLESHEET_SET=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samplesheet)
            [[ $# -lt 2 ]] && { echo "ERROR: --samplesheet requires a value." >&2; exit 1; }
            SAMPLESHEET="$2"
            SAMPLESHEET_SET=1
            shift 2
            ;;
        -1|--fastq1)
            [[ $# -lt 2 ]] && { echo "ERROR: --fastq1 requires a value." >&2; exit 1; }
            FASTQ_OVERRIDE_1="$2"
            shift 2
            ;;
        -2|--fastq2)
            [[ $# -lt 2 ]] && { echo "ERROR: --fastq2 requires a value." >&2; exit 1; }
            FASTQ_OVERRIDE_2="$2"
            shift 2
            ;;
        -n|--sample)
            [[ $# -lt 2 ]] && { echo "ERROR: --sample requires a value." >&2; exit 1; }
            SAMPLE_OVERRIDE="$2"
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
        -R|--read-source)
            [[ $# -lt 2 ]] && { echo "ERROR: --read-source requires a value." >&2; exit 1; }
            READ_SOURCE="$2"
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
if [[ -n "$FASTQ_OVERRIDE_1" || -n "$FASTQ_OVERRIDE_2" || -n "$SAMPLE_OVERRIDE" ]]; then
    USE_SAMPLESHEET=0
    if [[ -z "$FASTQ_OVERRIDE_1" || -z "$FASTQ_OVERRIDE_2" || -z "$SAMPLE_OVERRIDE" ]]; then
        echo "ERROR: --fastq1, --fastq2, and --sample must be provided together." >&2
        exit 1
    fi
    if (( SAMPLESHEET_SET )); then
        echo "ERROR: Cannot use --samplesheet with --fastq1/--fastq2/--sample." >&2
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
    for fq in "$FASTQ_OVERRIDE_1" "$FASTQ_OVERRIDE_2"; do
        if [[ ! -f "$fq" ]]; then
            echo "ERROR: FASTQ file '$fq' not found." >&2
            exit 1
        fi
    done
fi

echo "[INFO] Configuration:"
echo "[INFO]   Reference: $REF"
echo "[INFO]   Threads: $THREADS"
echo "[INFO]   AF cutoff: $AF_CUTOFF"
echo "[INFO]   Read source: $READ_SOURCE"
if (( USE_SAMPLESHEET )); then
    echo "[INFO]   Samplesheet: $SAMPLESHEET"
else
    echo "[INFO]   Single sample: $SAMPLE_OVERRIDE"
    echo "[INFO]     FASTQ1: $FASTQ_OVERRIDE_1"
    echo "[INFO]     FASTQ2: $FASTQ_OVERRIDE_2"
fi

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

process_sample() {
    local SAMPLE="$1"
    local GROUP="$2"
    local FASTQ_1="$3"
    local FASTQ_2="$4"

    (
        echo "[INFO] Processing sample: $SAMPLE (reads: $READ_SOURCE)"

        RAW_FASTQ_1="$FASTQ_1"
        RAW_FASTQ_2="$FASTQ_2"

        case "$READ_SOURCE" in
            trimmed)
                if [[ "$FASTQ_1" == *.trimmed.fastq.gz ]]; then
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

        if [[ ! -f "$FASTQ_1" || ! -f "$FASTQ_2" ]]; then
            echo "ERROR: FASTQ files for $SAMPLE not found!" >&2
            if [[ "$READ_SOURCE" == "trimmed" ]]; then
                echo "       Expected trimmed files at: $FASTQ_1 and $FASTQ_2" >&2
                echo "       Run scripts/fastqc.sh before scripts/evovariants.sh." >&2
                echo "       Original FASTQs provided: $RAW_FASTQ_1 and $RAW_FASTQ_2" >&2
            fi
            exit 1
        fi

        bwa-mem2 mem -t "$THREADS" "$REF" "$FASTQ_1" "$FASTQ_2" \
            2> "${RESULTS_DIR}/logs/${SAMPLE}.bwa2.err" \
          | samtools sort -@ "$THREADS" -n -O BAM \
            -o "${RESULTS_DIR}/bam/${SAMPLE}.namesort.bam" \
            2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_sort_namesort.err"

        samtools fixmate -@ "$THREADS" -m \
            "${RESULTS_DIR}/bam/${SAMPLE}.namesort.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.fixmate.bam" \
            2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_fixmate.err"

        # Drop secondary (0x100) and supplementary (0x800) records; they lack ms tags and break markdup.
        samtools view -@ "$THREADS" -b -F 0x900 \
            "${RESULTS_DIR}/bam/${SAMPLE}.fixmate.bam" \
            -o "${RESULTS_DIR}/bam/${SAMPLE}.fixmate.primary.bam" \
            2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_view_primary.err"

        samtools sort -@ "$THREADS" -O BAM \
            -o "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.fixmate.primary.bam" \
            2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_sort_coordinates.err"

        samtools markdup -@ "$THREADS" -r \
            "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.dedup.bam" \
            2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_markdup.err"

        rm -f \
            "${RESULTS_DIR}/bam/${SAMPLE}.namesort.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.fixmate.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.fixmate.primary.bam"

        lofreq indelqual --dindel -f "$REF" \
            -o "${RESULTS_DIR}/bam/${SAMPLE}.dedup.indelq.bam" \
            "${RESULTS_DIR}/bam/${SAMPLE}.dedup.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_indelqual.err"

        samtools index "${RESULTS_DIR}/bam/${SAMPLE}.dedup.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.samtools_index_dedup_indelq.err"

        rm -f "${RESULTS_DIR}/bam/${SAMPLE}.sorted.bam"

        set +e
        lofreq call-parallel \
            --pp-threads "$THREADS" \
            -f "$REF" \
            -o "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
            "${RESULTS_DIR}/bam/${SAMPLE}.dedup.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_call_parallel.err"
        rc=$?
        set -e

        if [ $rc -ne 0 ]; then
            echo "[WARN] lofreq call-parallel failed for ${SAMPLE} (see logs). Retrying with single-threaded lofreq call..."
            lofreq call \
                -f "$REF" \
                -o "${RESULTS_DIR}/vcf/${SAMPLE}.raw.vcf" \
                "${RESULTS_DIR}/bam/${SAMPLE}.dedup.indelq.bam" 2> "${RESULTS_DIR}/logs/${SAMPLE}.lofreq_call_single.err"
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

    while IFS=, read -r SAMPLE GROUP FASTQ_1 FASTQ_2; do
        [[ -z "$SAMPLE" ]] && continue
        process_sample "$SAMPLE" "$GROUP" "$FASTQ_1" "$FASTQ_2"
    done
}

if (( USE_SAMPLESHEET )); then
    process_samples_stream "$SAMPLESHEET" < "$SAMPLESHEET"
else
    process_samples_stream "single-sample input" <<EOF
sample,group,fastq_1,fastq_2
$SAMPLE_OVERRIDE,manual,$FASTQ_OVERRIDE_1,$FASTQ_OVERRIDE_2
EOF
fi

for pid in "${RUNNING_PIDS[@]}"; do
    wait "$pid"
    rc=$?
    if (( rc != 0 )); then
        exit $rc
    fi
done

echo "[INFO] All samples processed successfully."
