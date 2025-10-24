#!/bin/bash
set -euo pipefail

# PARAMETERS
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

RAW_DIR="${PROJECT_ROOT}/data/raw"
TRIM_DIR="${PROJECT_ROOT}/data/trimmed"
QC_DIR="${PROJECT_ROOT}/output/qc"
THREADS=18

usage() {
    cat <<'EOF'
Usage: scripts/fastqc.sh [options]

Options:
  -r, --raw-dir DIR    Directory containing raw FASTQs (default: project data/raw)
  -1, --fastq1 FILE    FASTQ R1 path (requires --fastq2)
  -2, --fastq2 FILE    FASTQ R2 path (requires --fastq1)
  -t, --threads N      Threads for fastqc/fastp (default: 18)
  -h, --help           Show this help message and exit

When --fastq1/--fastq2 are provided, only that pair is processed and RAW_DIR is ignored.
EOF
}

FASTQ_OVERRIDE_1=""
FASTQ_OVERRIDE_2=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--raw-dir)
            [[ $# -lt 2 ]] && { echo "ERROR: --raw-dir requires a value." >&2; exit 1; }
            RAW_DIR="$2"
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
        -t|--threads)
            [[ $# -lt 2 ]] && { echo "ERROR: --threads requires a value." >&2; exit 1; }
            THREADS="$2"
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

if ! [[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -ge 1 ]]; then
    echo "ERROR: --threads expects an integer >=1 (got '$THREADS')." >&2
    exit 1
fi

if [[ -n "$FASTQ_OVERRIDE_1" || -n "$FASTQ_OVERRIDE_2" ]]; then
    if [[ -z "$FASTQ_OVERRIDE_1" || -z "$FASTQ_OVERRIDE_2" ]]; then
        echo "ERROR: --fastq1 and --fastq2 must be provided together." >&2
        exit 1
    fi
    for fq in "$FASTQ_OVERRIDE_1" "$FASTQ_OVERRIDE_2"; do
        if [[ ! -f "$fq" ]]; then
            echo "ERROR: FASTQ file '$fq' not found." >&2
            exit 1
        fi
    done
else
    if [[ ! -d "$RAW_DIR" ]]; then
        echo "ERROR: RAW_DIR '$RAW_DIR' does not exist." >&2
        exit 1
    fi
fi

echo "[INFO] Running fastqc pipeline"
echo "[INFO]   Threads: $THREADS"
if [[ -n "$FASTQ_OVERRIDE_1" ]]; then
    echo "[INFO]   FASTQ1: $FASTQ_OVERRIDE_1"
    echo "[INFO]   FASTQ2: $FASTQ_OVERRIDE_2"
else
    echo "[INFO]   RAW_DIR: $RAW_DIR"
fi

shopt -s nullglob

for exe in fastqc fastp multiqc; do
    if ! command -v "$exe" >/dev/null 2>&1; then
        echo "ERROR: Required executable '$exe' not found in PATH" >&2
        exit 1
    fi
done

mkdir -p "$TRIM_DIR" "$QC_DIR"/raw "$QC_DIR"/trimmed

# Gather raw FASTQs and read pairs
raw_fastqs=()
pair_r1=()
pair_r2=()

if [[ -n "$FASTQ_OVERRIDE_1" ]]; then
    raw_fastqs+=("$FASTQ_OVERRIDE_1" "$FASTQ_OVERRIDE_2")
    pair_r1+=("$FASTQ_OVERRIDE_1")
    pair_r2+=("$FASTQ_OVERRIDE_2")
else
    raw_fastqs=("$RAW_DIR"/*.fastq.gz)
    if [ ${#raw_fastqs[@]} -eq 0 ]; then
        echo "ERROR: No FASTQ files matching *.fastq.gz found in $RAW_DIR" >&2
        exit 1
    fi

    for fq1 in "$RAW_DIR"/*_1.fastq.gz; do
        fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
        if [ ! -f "$fq2" ]; then
            echo "[WARN] Mate file for $fq1 not found (expected $fq2); skipping pair." >&2
            continue
        fi
        pair_r1+=("$fq1")
        pair_r2+=("$fq2")
    done

    if [ ${#pair_r1[@]} -eq 0 ]; then
        echo "ERROR: No FASTQ pairs matching *_1.fastq.gz/_2.fastq.gz found in $RAW_DIR" >&2
        exit 1
    fi
fi

# Run FastQC on raw files
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
trimmed_fastqs=()
for idx in "${!pair_r1[@]}"; do
    fq1="${pair_r1[$idx]}"
    fq2="${pair_r2[$idx]}"

    r1_name=$(basename "$fq1")
    r2_name=$(basename "$fq2")

    html_base=""
    out1=""
    out2=""

    if [[ "$r1_name" == *_1.fastq.gz && "$r2_name" == *_2.fastq.gz ]]; then
        sample_prefix="${r1_name%_1.fastq.gz}"
        html_base="$sample_prefix"
        out1="$TRIM_DIR/${sample_prefix}_1.trimmed.fastq.gz"
        out2="$TRIM_DIR/${sample_prefix}_2.trimmed.fastq.gz"
    elif [[ "$r1_name" == *_R1.fastq.gz && "$r2_name" == *_R2.fastq.gz ]]; then
        sample_prefix="${r1_name%_R1.fastq.gz}"
        html_base="$sample_prefix"
        out1="$TRIM_DIR/${sample_prefix}_R1.trimmed.fastq.gz"
        out2="$TRIM_DIR/${sample_prefix}_R2.trimmed.fastq.gz"
    else
        html_base="${r1_name%.fastq.gz}"
        out1="$TRIM_DIR/${r1_name%.fastq.gz}.trimmed.fastq.gz"
        out2="$TRIM_DIR/${r2_name%.fastq.gz}.trimmed.fastq.gz"
    fi

    fastp -i "$fq1" -I "$fq2" \
          -o "$out1" \
          -O "$out2" \
          -w "$THREADS" \
          --html "$QC_DIR/trimmed/${html_base}_fastp.html"

    trimmed_fastqs+=("$out1" "$out2")
done

# Run FastQC on trimmed files
if [ ${#trimmed_fastqs[@]} -eq 0 ]; then
    echo "[WARN] No trimmed FASTQ files to run FastQC on." >&2
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
