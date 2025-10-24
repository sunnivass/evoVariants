#!/usr/bin/env bash
set -euo pipefail

# Normalize LoFreq VCFs and annotate with SnpEff.
#
# Example:
#   scripts/normalize_and_annotate.sh \
#     --ref ref/custom.fa \
#     --genome-id sacCer3_strainA_v1 \
#     --cfg snpeff_data/snpeff.config \
#     --datadir snpeff_data \
#     --in-dir vcf_in \
#     --norm-dir work/normalized \
#     --ann-dir work/snpeff \
#     --threads 4 \
#     --heap 8g
#
# (Optional AF filter inline:)
#     ... --af-filter 'AF>=0.01'

usage() {
  cat <<EOF
Usage: $0 --ref REF.fa --genome-id ID [options]

Required:
  --ref REF.fa           Reference FASTA (same as used by LoFreq)
  --genome-id ID         SnpEff genome ID you built (e.g., sacCer3_strainA_v1)

Recommended:
  --cfg FILE             repo-local snpEff.config (default: snpeff_data/snpEff.config)
  --datadir DIR          SnpEff data dir (default: snpeff_data)

I/O:
  --in-dir DIR           Input directory with .vcf or .vcf.gz (default: output/vcf)
  --norm-dir DIR         Output directory for normalized VCFs (default: output/normalized)
  --ann-dir DIR          Output directory for snpEff VCFs (default: output/snpeff)

Performance:
  --threads N            Threads for bcftools/bgzip where supported (default: 4)
  --heap SIZE            Java heap for snpEff, e.g. 8g (default: 8g)

Optional:
  --af-filter EXPR       bcftools -i filter (e.g. 'AF>=0.01'); skipped if not provided
  -h, --help             Show help

Notes:
- Normalization: sort -> index -> left-align/split (-m -any)
- Output files keep sample basename, with .norm.vcf.gz then .snpeff.vcf.gz
EOF
}

REF=""
GENOME_ID=""
CFG="snpeff_data/snpeff.config"
DATADIR="snpeff_data"
IN_DIR="output/vcf"
NORM_DIR="output/normalized"
ANN_DIR="output/snpeff"
THREADS="4"
HEAP="8g"
AF_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --genome-id) GENOME_ID="$2"; shift 2 ;;
    --cfg) CFG="$2"; shift 2 ;;
    --datadir) DATADIR="$2"; shift 2 ;;
    --in-dir) IN_DIR="$2"; shift 2 ;;
    --norm-dir) NORM_DIR="$2"; shift 2 ;;
    --ann-dir) ANN_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --heap) HEAP="$2"; shift 2 ;;
    --af-filter) AF_FILTER="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$REF" || -z "$GENOME_ID" ]] && { echo "Error: --ref and --genome-id are required."; usage; exit 1; }
[[ -f "$REF" ]] || { echo "Error: REF not found: $REF"; exit 1; }
[[ -f "$CFG" ]] || { echo "Error: snpEff config not found: $CFG"; exit 1; }
[[ -d "$DATADIR" ]] || { echo "Error: snpEff data dir not found: $DATADIR"; exit 1; }

# Resolve config/data paths so snpEff never constructs bogus nested paths
if [[ "$CFG" != /* ]]; then
  CFG="$(cd "$(dirname "$CFG")" && pwd)/$(basename "$CFG")"
fi
if [[ "$DATADIR" != /* ]]; then
  DATADIR="$(cd "$DATADIR" && pwd)"
fi

# tool checks
command -v bcftools >/dev/null || { echo "bcftools not found in PATH"; exit 1; }
command -v snpEff >/dev/null   || { echo "snpEff not found in PATH"; exit 1; }
command -v bgzip >/dev/null    || { echo "bgzip (htslib) not found in PATH"; exit 1; }
command -v tabix >/dev/null    || { echo "tabix not found in PATH"; exit 1; }

mkdir -p "$NORM_DIR" "$ANN_DIR"

shopt -s nullglob
inputs=("$IN_DIR"/*.vcf "$IN_DIR"/*.vcf.gz)
if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "No input VCFs found in $IN_DIR (expecting .vcf or .vcf.gz)"
  exit 1
fi

echo ">>> Normalizing ($(date))"
for V in "${inputs[@]}"; do
  base=$(basename "$V")
  out="$NORM_DIR/${base%.vcf.gz}"
  out="${out%.vcf}.norm.vcf.gz"

  echo "  - $base -> $(basename "$out")"
  # sort & compress
  bcftools sort -Oz -o "$out.tmp.gz" "$V"
  tabix -f "$out.tmp.gz"

  # left-align & split multi-allelics
  if [[ -n "$AF_FILTER" ]]; then
    # apply AF filter inline
    bcftools norm --threads "$THREADS" -f "$REF" -m -any "$out.tmp.gz" \
      | bcftools view --threads "$THREADS" -i "$AF_FILTER" -Oz -o "$out"
  else
    bcftools norm --threads "$THREADS" -f "$REF" -m -any "$out.tmp.gz" -Oz -o "$out"
  fi
  tabix -f "$out"
  rm -f "$out.tmp.gz" "$out.tmp.gz".{tbi,csi} 2>/dev/null || true
done

echo ">>> Annotating with SnpEff genome: $GENOME_ID ($(date))"
declare -a sample_tsvs=()
for V in "$NORM_DIR"/*.vcf.gz; do
  b=$(basename "$V" .vcf.gz)
  out="$ANN_DIR/${b}.snpeff.vcf.gz"
  echo "  - $b.vcf.gz -> $(basename "$out")"
  snpEff -Xmx"$HEAP" -c "$CFG" -dataDir "$DATADIR" -v "$GENOME_ID" -noStats "$V" \
    | bgzip -@ "$THREADS" -c > "$out"
  tabix -f "$out"

  sample_name="${b%.norm}"
  sample_tsv="$ANN_DIR/${sample_name}.snpeff.tsv"
  echo "    TSV -> $(basename "$sample_tsv")"
  {
    printf 'SAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tANN\tDP\tAF\n'
    bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/ANN\t%INFO/DP\t%INFO/AF\n' "$out" \
      | awk -v sample="$sample_name" 'BEGIN{OFS="\t"} {print sample, $0}'
  } > "$sample_tsv"
  sample_tsvs+=("$sample_tsv")
done

if [[ ${#sample_tsvs[@]} -gt 0 ]]; then
  merged_tsv="$ANN_DIR/merged.snpeff.tsv"
  echo ">>> Building merged TSV (wide): $(basename "$merged_tsv")"
  python3 - <<'PY' "$merged_tsv" "${sample_tsvs[@]}"
import csv
import os
import sys
from collections import OrderedDict

merged_path = sys.argv[1]
sample_files = sys.argv[2:]

base_cols = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"]
ann_cols = [
    "ANN_EFFECT",
    "ANN_IMPACT",
    "ANN_GENE_NAME",
    "ANN_GENE_ID",
    "ANN_BIOTYPE",
    "ANN_HGVS_P",
]
key_cols = ["CHROM", "POS", "REF", "ALT"]

impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def summarize_ann(raw_ann: str):
    default_summary = {col: "" for col in ann_cols}
    best_summary = default_summary
    best_score = float("inf")

    if not raw_ann:
        return best_summary, best_score

    for entry in raw_ann.split(","):
        parts = entry.split("|")
        if len(parts) < 16:
            parts += [""] * (16 - len(parts))

        annotation = {
            "ANN_EFFECT": parts[1],
            "ANN_IMPACT": parts[2],
            "ANN_GENE_NAME": parts[3],
            "ANN_GENE_ID": parts[4],
            "ANN_BIOTYPE": parts[7],
            "ANN_HGVS_P": parts[10],
        }
        score = impact_rank.get(annotation["ANN_IMPACT"], float("inf"))

        if score < best_score:
            best_summary = annotation
            best_score = score

    return best_summary, best_score


variants = OrderedDict()
sample_order = []


def add_sample(name: str) -> None:
    if name not in sample_order:
        sample_order.append(name)

for sample_file in sample_files:
    sample_name = os.path.basename(sample_file).replace('.snpeff.tsv', '')
    add_sample(sample_name)

    with open(sample_file, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None:
            continue

        for row in reader:
            sample_from_row = row.get('SAMPLE') or sample_name
            add_sample(sample_from_row)

            key = tuple(row.get(col, "") for col in key_cols)
            ann_summary, ann_score = summarize_ann(row.get('ANN', ""))
            entry = variants.get(key)
            if entry is None:
                entry = {
                    'base': {col: row.get(col, "") for col in base_cols},
                    'ann': ann_summary,
                    'ann_score': ann_score,
                    'values': {},
                }
                variants[key] = entry
            else:
                if ann_score < entry.get('ann_score', float('inf')):
                    entry['ann'] = ann_summary
                    entry['ann_score'] = ann_score

            entry['values'][sample_from_row] = (row.get('DP', ''), row.get('AF', ''))

header = base_cols[:]
header.extend(ann_cols)
for sample in sample_order:
    header.append(f"{sample}_DP")
    header.append(f"{sample}_AF")

with open(merged_path, 'w', newline='') as out_handle:
    writer = csv.writer(out_handle, delimiter='\t', lineterminator='\n')
    writer.writerow(header)

    for entry in variants.values():
        row_out = [entry['base'][col] for col in base_cols]
        ann_values = entry['ann'] if 'ann' in entry else {col: '' for col in ann_cols}
        row_out.extend(ann_values.get(col, '') for col in ann_cols)
        for sample in sample_order:
            dp, af = entry['values'].get(sample, ('', ''))
            row_out.extend([dp, af])
        writer.writerow(row_out)
PY
fi

echo "All done. Normalized in: $NORM_DIR ; Annotated in: $ANN_DIR"
