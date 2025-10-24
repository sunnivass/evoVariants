#!/usr/bin/env python3
"""Annotate LoFreq VCF files and merge per-sample metrics."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

FORMAT_DP_HEADER = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for this sample\">\n"
FORMAT_AF_HEADER = "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for this sample\">\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate LoFreq-produced VCF files, add per-sample DP/AF fields, and create merged summaries.",
    )
    parser.add_argument(
        "--vcf-dir",
        type=Path,
        default=Path("output/vcf"),
        help="Directory containing per-sample VCF files exported from LoFreq (default: output/vcf).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/annotated_vcf"),
        help="Destination directory for annotated outputs (default: output/annotated_vcf).",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("ref/fasta/sacCer3_SFS01decoy_HsGTaseCeg1.masked2.fa"),
        help="Reference genome FASTA used for annotation.",
    )
    parser.add_argument(
        "--gff",
        type=Path,
        default=Path("ref/gff/genomic_HsGTase-CEG1.gff"),
        help="Gene annotation GFF3 used by bcftools csq.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        help="Optional CSV file (sample,...) to control sample ordering and filtering.",
    )
    parser.add_argument(
        "--vcf-suffix",
        default=".filtered.vcf",
        help="Expected suffix appended to sample names for LoFreq VCFs (default: .filtered.vcf).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to pass to bcftools where supported (default: 1).",
    )
    return parser.parse_args()


def ensure_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        sys.stderr.write(f"[ERROR] Required tool '{tool}' is not on PATH.\n")
        sys.exit(1)


def run_command(cmd: List[str], **kwargs) -> None:
    subprocess.run(cmd, check=True, **kwargs)


def load_sample_order(args: argparse.Namespace) -> List[str]:
    if args.sample_sheet:
        with args.sample_sheet.open(newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
            if header is None:
                raise SystemExit(f"[ERROR] Sample sheet {args.sample_sheet} is empty.")
            samples = [row[0] for row in reader if row]
        if not samples:
            raise SystemExit(f"[ERROR] No samples found in {args.sample_sheet}.")
        return samples

    vcfs = sorted(args.vcf_dir.glob("*.vcf")) + sorted(args.vcf_dir.glob("*.vcf.gz"))
    samples: List[str] = []
    for vcf in vcfs:
        name = vcf.name
        if name.endswith(".vcf.gz"):
            name = name[:-7]
        elif name.endswith(".vcf"):
            name = name[:-4]
        if name.endswith(args.vcf_suffix):
            name = name[: -len(args.vcf_suffix)]
        samples.append(name)
    if not samples:
        raise SystemExit(f"[ERROR] No VCF files found in {args.vcf_dir}.")
    return samples


def locate_vcf(sample: str, args: argparse.Namespace) -> Path:
    candidates = [
        args.vcf_dir / f"{sample}{args.vcf_suffix}",
        args.vcf_dir / f"{sample}{args.vcf_suffix}.gz",
        args.vcf_dir / f"{sample}.vcf",
        args.vcf_dir / f"{sample}.vcf.gz",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise SystemExit(f"[ERROR] Could not locate VCF for sample '{sample}' in {args.vcf_dir}.")


def parse_info_field(info: str) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for entry in info.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            values[key] = value
        else:
            values[entry] = "True"
    return values


def inject_format_fields(sample: str, raw_path: Path, final_path: Path) -> None:
    with raw_path.open("rt") as fin, final_path.open("wt") as fout:
        seen_dp = False
        seen_af = False
        for line in fin:
            if line.startswith("##"):
                if line.startswith("##FORMAT=<ID=DP"):
                    seen_dp = True
                elif line.startswith("##FORMAT=<ID=AF"):
                    seen_af = True
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                if not seen_dp:
                    fout.write(FORMAT_DP_HEADER)
                if not seen_af:
                    fout.write(FORMAT_AF_HEADER)
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= 8:
                    fields.extend(["FORMAT", sample])
                elif len(fields) == 9:
                    fields.append(sample)
                else:
                    fields[9] = sample
                fout.write("\t".join(fields) + "\n")
                break

        for line in fin:
            if not line or line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                fout.write(line)
                continue
            info_values = parse_info_field(fields[7])
            dp_value = info_values.get("DP", ".")
            af_value = info_values.get("AF", ".")
            base = fields[:8]
            base.append("DP:AF")
            base.append(f"{dp_value}:{af_value}")
            fout.write("\t".join(base) + "\n")


def annotate_sample(sample: str, vcf_path: Path, args: argparse.Namespace, out_dir: Path) -> Tuple[Path, Path]:
    raw_path = out_dir / f"{sample}.annotated.raw.vcf"
    final_vcf = out_dir / f"{sample}.annotated.vcf"
    final_gz = out_dir / f"{sample}.annotated.vcf.gz"
    sample_tsv = out_dir / f"{sample}.annotated.tsv"

    cmd = [
        "bcftools",
        "csq",
        "--fasta-ref",
        str(args.reference),
        "--gff-annot",
        str(args.gff),
        "--threads",
        str(args.threads),
        "--output-type",
        "v",
        "--output",
        str(raw_path),
        str(vcf_path),
    ]
    run_command(cmd)

    inject_format_fields(sample, raw_path, final_vcf)
    raw_path.unlink(missing_ok=True)

    run_command(["bgzip", "-f", str(final_vcf)])
    run_command(["bcftools", "index", str(final_gz)])

    with sample_tsv.open("wt") as handle:
        handle.write("CHROM\tPOS\tREF\tALT\tBCSQ\tDP\tAF\n")
        query_fmt = "%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\t[%DP]\t[%AF]\n"
        run_command([
            "bcftools",
            "query",
            "--format",
            query_fmt,
            str(final_gz),
        ], stdout=handle)

    return final_gz, sample_tsv


def ensure_reference_index(reference: Path) -> None:
    fai = reference.with_suffix(reference.suffix + ".fai")
    if not fai.exists():
        run_command(["samtools", "faidx", str(reference)])


def merge_annotations(samples: List[str], annotated_vcfs: List[Path], out_dir: Path, threads: int) -> Tuple[Path, Path]:
    merged_vcf = out_dir / "merged.vcf.gz"
    merged_tsv = out_dir / "merged.variants.tsv"

    merge_cmd = [
        "bcftools",
        "merge",
        "--threads",
        str(threads),
        "--output-type",
        "z",
        "--output",
        str(merged_vcf),
    ] + [str(vcf) for vcf in annotated_vcfs]
    run_command(merge_cmd)
    run_command(["bcftools", "index", str(merged_vcf)])

    with merged_tsv.open("wt") as handle:
        handle.write("CHROM\tPOS\tREF\tALT\tBCSQ")
        for sample in samples:
            handle.write(f"\t{sample}_DP\t{sample}_AF")
        handle.write("\n")
        query_fmt = "%CHROM\t%POS\t%REF\t%ALT\t%BCSQ[\t%DP\t%AF]\n"
        run_command([
            "bcftools",
            "query",
            "--format",
            query_fmt,
            str(merged_vcf),
        ], stdout=handle)

    return merged_vcf, merged_tsv


def main() -> None:
    args = parse_args()
    ensure_tool("bcftools")
    ensure_tool("samtools")
    ensure_tool("bgzip")

    if not args.vcf_dir.is_dir():
        raise SystemExit(f"[ERROR] VCF directory {args.vcf_dir} does not exist or is not a directory.")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not args.reference.exists():
        raise SystemExit(f"[ERROR] Reference FASTA {args.reference} does not exist.")
    if not args.gff.exists():
        raise SystemExit(f"[ERROR] GFF3 file {args.gff} does not exist.")

    ensure_reference_index(args.reference)

    samples = load_sample_order(args)

    annotated_vcfs: List[Path] = []
    for sample in samples:
        vcf_path = locate_vcf(sample, args)
        print(f"[INFO] Annotating {sample} from {vcf_path}")
        annotated_vcf, sample_tsv = annotate_sample(sample, vcf_path, args, args.out_dir)
        annotated_vcfs.append(annotated_vcf)
        print(f"[INFO] Wrote {annotated_vcf} and {sample_tsv}")

    if not annotated_vcfs:
        raise SystemExit("[ERROR] No annotated VCFs generated; aborting merge.")

    print("[INFO] Merging annotated VCFs")
    merged_vcf, merged_tsv = merge_annotations(samples, annotated_vcfs, args.out_dir, args.threads)
    print(f"[INFO] Merged VCF: {merged_vcf}")
    print(f"[INFO] Merged TSV: {merged_tsv}")


if __name__ == "__main__":
    main()
