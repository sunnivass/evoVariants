#!/usr/bin/env bash
set -euo pipefail

# Build a custom SnpEff genome from FASTA + GFF3.
# Works with a repo-local snpEff.config + data dir so you never touch conda's global files.
#
# Example:
#   scripts/snpeff_build.sh \
#     --genome-id sacCer3_strainA_v1 \
#     --fasta ref/strainA.fa \
#     --gff ref/strainA.gff3 \
#     --data-dir snpeff_data \
#     --config snpeff_data/snpEff.config \
#     --display-name "sacCer3 strain A (v1)"
#
# Then annotate with:
#   snpEff -c snpeff_data/snpEff.config -dataDir snpeff_data -v sacCer3_strainA_v1 in.vcf.gz | bgzip -c > out.vcf.gz

usage() {
  cat <<EOF
Usage: $0 --genome-id ID --fasta PATH.fa --gff PATH.gff3 [options]

Required:
  --genome-id ID         Unique ID for this build (used by SnpEff). e.g., sacCer3_strainA_v1
  --fasta PATH.fa        Reference FASTA
  --gff PATH.gff3        Matching GFF3 (seqids must match FASTA contigs)

Options:
  --data-dir DIR         SnpEff data directory (default: snpeff_data)
  --config FILE          Path to repo-local snpEff.config (default: snpeff_data/snpEff.config)
  --display-name NAME    Human-readable name in config (default: same as --genome-id)
  --force                Overwrite existing genome folder if it exists
  --no-check             Skip contig-name sanity checks
  -h, --help             Show this help

Notes:
- This script copies FASTA/GFF into: <data-dir>/<genome-id>/{sequences.fa,genes.gff}
- It appends "<genome-id>.genome : <display-name>" to the local snpEff.config (if missing).
- Requires: snpEff, java, md5sum, awk, sort, comm
EOF
}

GENOME_ID=""
FASTA=""
GFF=""
DATA_DIR="snpeff_data"
CONFIG=""
DISPLAY_NAME=""
FORCE="0"
NOCHECK="0"

resolve_path() {
  local target="$1"
  if command -v realpath >/dev/null; then
    if realpath -m "$target" 2>/dev/null; then
      return
    fi
    if realpath "$target" 2>/dev/null; then
      return
    fi
  fi
  if command -v python3 >/dev/null; then
    python3 - "$target" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
    return
  fi
  if command -v python >/dev/null; then
    python - "$target" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
    return
  fi
  if [[ "$target" == /* ]]; then
    echo "$target"
  else
    echo "$(pwd -P)/$target"
  fi
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome-id) GENOME_ID="$2"; shift 2 ;;
    --fasta) FASTA="$2"; shift 2 ;;
    --gff) GFF="$2"; shift 2 ;;
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --display-name) DISPLAY_NAME="$2"; shift 2 ;;
    --force) FORCE="1"; shift ;;
    --no-check) NOCHECK="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "$GENOME_ID" || -z "$FASTA" || -z "$GFF" ]] && { echo "Error: --genome-id, --fasta, --gff are required." >&2; usage; exit 1; }
[[ -f "$FASTA" ]] || { echo "Error: FASTA not found: $FASTA" >&2; exit 1; }
[[ -f "$GFF" ]] || { echo "Error: GFF3 not found: $GFF" >&2; exit 1; }

CONFIG="${CONFIG:-${DATA_DIR}/snpEff.config}"
DISPLAY_NAME="${DISPLAY_NAME:-${GENOME_ID}}"

mkdir -p "$DATA_DIR"

# Create / update local snpEff.config
if [[ ! -f "$CONFIG" ]]; then
  echo "Creating local snpEff.config at: $CONFIG"
  cat > "$CONFIG" <<EOF
# Repo-local SnpEff config
data.dir = ${DATA_DIR}
EOF
fi

# Ensure data.dir points to DATA_DIR (idempotent)
if ! grep -qE "^\s*data\.dir\s*=\s*${DATA_DIR}\s*$" "$CONFIG"; then
  # comment out existing data.dir lines, then append ours
  sed -i.bak 's/^\s*data\.dir\s*=.*/# &/' "$CONFIG"
  rm -f "$CONFIG.bak"
  echo "data.dir = ${DATA_DIR}" >> "$CONFIG"
fi

# Add genome ID mapping if missing
if ! grep -qE "^${GENOME_ID}\.genome\s*:" "$CONFIG"; then
  echo "${GENOME_ID}.genome : ${DISPLAY_NAME}" >> "$CONFIG"
else
  echo "NOTE: ${GENOME_ID}.genome is already present in $CONFIG (leaving as-is)."
fi

# Prepare genome folder
GEN_DIR="${DATA_DIR}/${GENOME_ID}"
if [[ -d "$GEN_DIR" ]]; then
  if [[ "$FORCE" == "1" ]]; then
    echo "WARN: Removing existing folder: $GEN_DIR"
    rm -rf "$GEN_DIR"
  else
    echo "Error: Genome folder already exists: $GEN_DIR (use --force to overwrite)" >&2
    exit 1
  fi
fi
mkdir -p "$GEN_DIR"

# Copy inputs with SnpEff's expected names
cp -f "$FASTA" "${GEN_DIR}/sequences.fa"
cp -f "$GFF"   "${GEN_DIR}/genes.gff"

# Optional: sanity-check contig names
if [[ "$NOCHECK" != "1" ]]; then
  echo "Checking contig names (FASTA vs GFF)..."
  FASTA_CONTIGS=$(awk '/^>/{gsub(/^>/,"",$1); print $1}' "${GEN_DIR}/sequences.fa" | sed 's/\r$//' | sort -u)
  GFF_SEQIDS=$(awk 'BEGIN{FS="\t"} !/^#/ && NF>=1 {print $1}' "${GEN_DIR}/genes.gff" | sed 's/\r$//' | sort -u)

  # Write to temp files
  TMPD=$(mktemp -d)
  echo "$FASTA_CONTIGS" > "$TMPD/f.fa"
  echo "$GFF_SEQIDS"   > "$TMPD/f.gff"
  # normalize blanks
  sed -i '/^$/d' "$TMPD/f.fa" "$TMPD/f.gff"

  ONLY_FASTA=$(comm -23 "$TMPD/f.fa" "$TMPD/f.gff" || true)
  ONLY_GFF=$(comm -13 "$TMPD/f.fa" "$TMPD/f.gff" || true)

  if [[ -n "$ONLY_FASTA" || -n "$ONLY_GFF" ]]; then
    echo "WARNING: Contig name mismatch detected."
    if [[ -n "$ONLY_FASTA" ]]; then
      echo "  Present in FASTA only:"
      echo "$ONLY_FASTA" | sed 's/^/    - /'
    fi
    if [[ -n "$ONLY_GFF" ]]; then
      echo "  Present in GFF only:"
      echo "$ONLY_GFF" | sed 's/^/    - /'
    fi
    echo "If these are expected (e.g., decoy contigs), you can ignore. Otherwise, fix names or use --no-check to skip."
  else
    echo "Contig names look consistent."
  fi
  rm -rf "$TMPD"
fi

# Build the DB
echo "Building SnpEff DB: ${GENOME_ID}"
snpEff build -c "$CONFIG" -dataDir "$DATA_DIR" -gff3 -v "$GENOME_ID"

# Record provenance
SNPEFF_VER=$(snpEff -version 2>/dev/null | awk '{print $NF}')
JAVA_VER=$(java -version 2>&1 | head -n1 | sed 's/\r$//')
FA_MD5=$(md5sum "${GEN_DIR}/sequences.fa" | awk '{print $1}')
GFF_MD5=$(md5sum "${GEN_DIR}/genes.gff"   | awk '{print $1}')
DATE=$(date -Iseconds)
FA_RESOLVED=$(resolve_path "$FASTA")
GFF_RESOLVED=$(resolve_path "$GFF")

cat > "${GEN_DIR}/BUILD_INFO.txt" <<EOF
Date: ${DATE}
Genome ID: ${GENOME_ID}
Display name: ${DISPLAY_NAME}
SnpEff version: ${SNPEFF_VER}
Java: ${JAVA_VER}
FASTA source: ${FA_RESOLVED}
FASTA md5: ${FA_MD5}
GFF3 source: ${GFF_RESOLVED}
GFF3 md5: ${GFF_MD5}
Command: snpEff build -c ${CONFIG} -dataDir ${DATA_DIR} -gff3 -v ${GENOME_ID}
EOF

echo "Done."
echo
echo "Use this genome with:"
echo "  snpEff -c ${CONFIG} -dataDir ${DATA_DIR} -v ${GENOME_ID} input.vcf.gz | bgzip -c > annotated.vcf.gz"
