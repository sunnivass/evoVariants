#!/bin/bash

set -euo pipefail

# Ensure required tools are available
if ! command -v bcftools >/dev/null 2>&1; then
    echo "[ERROR] bcftools is not available on PATH." >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "[ERROR] samtools is not available on PATH." >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 is not available on PATH." >&2
    exit 1
fi

# Parameters
REF_FASTA="ref/fasta/sacCer3_SFS01decoy_HsGTaseCeg1.masked2.fa"
GFF3="ref/gff/genomic_HsGTase-CEG1.gff"
SAMPLESHEET="data/samplesheet_HsCeg1.csv"
VCF_DIR="output/vcf"
OUT_DIR="output/annotated_vcf"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Track per-sample outputs for downstream aggregation
declare -a SAMPLE_NAMES=()
declare -A SAMPLE_TSV_PATHS=()

# Index reference if not already
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Indexing reference genome..."
    samtools faidx "$REF_FASTA"
fi

# Annotate with bcftools csq
# Read the CSV file, skip the header - modified to avoid subshell issues
while IFS=, read -r SAMPLE GROUP FASTQ1 FASTQ2; do
    IN_VCF="$VCF_DIR/${SAMPLE}.filtered.vcf"
    OUT_VCF="$OUT_DIR/${SAMPLE}.annotated.vcf.gz"
    OUT_TSV="$OUT_DIR/${SAMPLE}.annotated.tsv"

    echo "[INFO] Annotating variants for ${SAMPLE}..."
    #bcftools reheader -s <(echo "$SAMPLE") "$IN_VCF" | \
    bcftools csq \
        -f "$REF_FASTA" \
        -g "$GFF3" \
        -v \
        -Oz -o "$OUT_VCF" \
        "$IN_VCF"

    echo "[INFO] Adding per-sample FORMAT fields for ${SAMPLE}..."
    TMP_VCF_PLAIN="$OUT_DIR/${SAMPLE}.annotated.with_format.vcf"
    python3 - "$SAMPLE" "$OUT_VCF" "$TMP_VCF_PLAIN" <<'PY'
import gzip
import sys

sample, in_path, out_path = sys.argv[1:]

format_dp_header = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for this sample">\n'
format_af_header = '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency for this sample">\n'

def parse_info(info_field):
    values = {}
    for entry in info_field.split(';'):
        if not entry:
            continue
        if '=' in entry:
            key, value = entry.split('=', 1)
            values[key] = value
        else:
            values[entry] = True
    return values

with gzip.open(in_path, 'rt') as fin, open(out_path, 'w') as fout:
    emitted_format_headers = False
    seen_format_dp = False
    seen_format_af = False

    for line in fin:
        if line.startswith('##'):
            if line.startswith('##FORMAT=<ID=DP'):
                seen_format_dp = True
            elif line.startswith('##FORMAT=<ID=AF'):
                seen_format_af = True
            fout.write(line)
            continue

        if line.startswith('#CHROM'):
            if not seen_format_dp:
                fout.write(format_dp_header)
            if not seen_format_af:
                fout.write(format_af_header)

            header_fields = line.rstrip('\n').split('\t')
            if len(header_fields) <= 8:
                header_line = line.rstrip('\n') + '\tFORMAT\t' + sample
            else:
                header_fields[9] = sample
                header_line = '\t'.join(header_fields)

            fout.write(header_line + '\n')
            emitted_format_headers = True
            continue

        if not emitted_format_headers:
            raise RuntimeError('FORMAT headers not emitted before records')

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 8:
            raise RuntimeError('Unexpected VCF column count: ' + str(len(fields)))

        info = parse_info(fields[7])
        dp_value = info.get('DP', '.')
        af_value = info.get('AF', '.')

        format_spec = 'DP:AF'
        sample_values = f"{dp_value}:{af_value}"

        fout.write('\t'.join(fields[:8]) + f"\t{format_spec}\t{sample_values}\n")
PY

    bgzip -c "$TMP_VCF_PLAIN" > "${OUT_VCF}.tmp"
    mv -f "${OUT_VCF}.tmp" "$OUT_VCF"
    rm -f "$TMP_VCF_PLAIN"

    # Index the compressed VCF (required for merging)
    bcftools index "$OUT_VCF"

    SAMPLE_NAMES+=("$SAMPLE")
    SAMPLE_TSV_PATHS["$SAMPLE"]="$OUT_TSV"

    # Convert to TSV for individual sample (include header)
    printf 'CHROM\tPOS\tREF\tALT\tBCSQ\tDP\tAF\n' > "$OUT_TSV"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\t[%DP]\t[%AF]\n' "$OUT_VCF" >> "$OUT_TSV"

    echo "[INFO] Annotation completed: ${SAMPLE}"
done < <(tail -n +2 "$SAMPLESHEET")

if [ ${#SAMPLE_NAMES[@]} -eq 0 ]; then
    echo "[ERROR] No samples found in $SAMPLESHEET." >&2
    exit 1
fi

# Merge all annotated VCFs
shopt -s nullglob
ANNOTATED_VCFS=("$OUT_DIR"/*.annotated.vcf.gz)
shopt -u nullglob

if [ ${#ANNOTATED_VCFS[@]} -eq 0 ]; then
    echo "[ERROR] No annotated VCFs found in $OUT_DIR; aborting merge." >&2
    exit 1
fi

echo "[INFO] Merging annotated VCFs..."
bcftools merge \
    -O z \
    -o "$OUT_DIR/merged.vcf.gz" \
    "${ANNOTATED_VCFS[@]}"

# Index merged VCF
bcftools index "$OUT_DIR/merged.vcf.gz"

# Convert merged VCF into a combined TSV with per-sample DP/AF columns
MERGED_VARIANTS_TSV="$OUT_DIR/merged.variants.tsv"
{
    printf 'CHROM\tPOS\tREF\tALT\tBCSQ'
    for SAMPLE in "${SAMPLE_NAMES[@]}"; do
        printf '\t%s_DP\t%s_AF' "$SAMPLE" "$SAMPLE"
    done
    printf '\n'
} > "$MERGED_VARIANTS_TSV"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ[\t%DP\t%AF]\n' \
    "$OUT_DIR/merged.vcf.gz" \
    >> "$MERGED_VARIANTS_TSV"

echo "[INFO] Combined TSV written to $MERGED_VARIANTS_TSV"
