#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class Thresholds:
    similarity_diff: float
    increase_factor: float
    presence_threshold: float
    max_control_presence: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter SnpEff TSV by treatment-specific or increased-frequency variants "
            "and output SNP and gene-set tables."
        )
    )
    parser.add_argument(
        "--input",
        default="output/snpeff/merged.snpeff.tsv",
        help="Input merged SnpEff TSV.",
    )
    parser.add_argument(
        "--output-snps",
        default="output/snpeff/filtered.snpeff.tsv",
        help="Output SNP TSV.",
    )
    parser.add_argument(
        "--output-genes",
        default="output/snpeff/filtered.snpeff.genes.tsv",
        help="Output gene-set TSV.",
    )
    parser.add_argument(
        "--similarity-diff",
        type=float,
        default=0.05,
        help="Absolute AF difference to control-par treated as similar.",
    )
    parser.add_argument(
        "--increase-factor",
        type=float,
        default=2.0,
        help="Treatment AF must be >= mean control AF * factor to count as increased.",
    )
    parser.add_argument(
        "--presence-threshold",
        type=float,
        default=0.01,
        help="AF threshold for presence/absence logic.",
    )
    parser.add_argument(
        "--max-control-presence",
        type=int,
        default=1,
        help="Max number of control1/2/3 samples allowed for treatment-only variants.",
    )
    parser.add_argument(
        "--min-dp",
        type=int,
        default=30,
        help="Minimum treatment DP required to keep a variant.",
    )
    parser.add_argument(
        "--control-par",
        default="control-par",
        help="Control-par sample name (prefix before .filtered_AF).",
    )
    parser.add_argument(
        "--controls",
        default="control1,control2,control3",
        help="Comma-separated control sample names used for mean AF.",
    )
    return parser.parse_args()


def load_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df


def select_af_columns(df: pd.DataFrame) -> dict[str, str]:
    return {col.replace(".filtered_AF", ""): col for col in df.columns if col.endswith(".filtered_AF")}


def select_dp_columns(df: pd.DataFrame) -> dict[str, str]:
    return {col.replace(".filtered_DP", ""): col for col in df.columns if col.endswith(".filtered_DP")}


def main() -> None:
    args = parse_args()
    thresholds = Thresholds(
        similarity_diff=args.similarity_diff,
        increase_factor=args.increase_factor,
        presence_threshold=args.presence_threshold,
        max_control_presence=args.max_control_presence,
    )

    df = load_table(args.input)
    af_cols = select_af_columns(df)
    dp_cols = select_dp_columns(df)

    control_par_name = args.control_par
    control_names = [name for name in args.controls.split(",") if name]
    treatment_names = [name for name in af_cols if name not in control_names + [control_par_name]]

    missing = [name for name in control_names + [control_par_name] if name not in af_cols]
    if missing:
        raise SystemExit(f"Missing AF columns for: {', '.join(missing)}")

    missing_dp = [name for name in treatment_names if name not in dp_cols]
    if missing_dp:
        raise SystemExit(f"Missing DP columns for treatments: {', '.join(missing_dp)}")

    control_cols = [af_cols[name] for name in control_names]
    control_par_col = af_cols[control_par_name]
    treatment_cols = [af_cols[name] for name in treatment_names]
    treatment_dp_cols = [dp_cols[name] for name in treatment_names]

    df["_max_treatment_af"] = df[treatment_cols].max(axis=1, skipna=True)
    df["_mean_control_af"] = df[control_cols].mean(axis=1, skipna=True)
    df["_control_present_count"] = (df[control_cols] > thresholds.presence_threshold).sum(axis=1)
    df["_max_treatment_dp"] = df[treatment_dp_cols].max(axis=1, skipna=True)

    df["_treatment_present"] = df["_max_treatment_af"] > thresholds.presence_threshold
    df["_treatment_only"] = df["_treatment_present"] & (
        df["_control_present_count"] <= thresholds.max_control_presence
    )

    df["_increased"] = False
    increased_mask = df["_mean_control_af"] > 0
    df.loc[increased_mask, "_increased"] = (
        df.loc[increased_mask, "_max_treatment_af"]
        >= df.loc[increased_mask, "_mean_control_af"] * thresholds.increase_factor
    )

    df["_similar_to_control_par"] = (
        (df[control_par_col].notna())
        & (df["_max_treatment_af"].notna())
        & ((df["_max_treatment_af"] - df[control_par_col]).abs() < thresholds.similarity_diff)
    )

    nonsyn_mask = ~df["ANN_EFFECT"].astype(str).str.contains("synonymous_variant", na=False)
    keep_mask = (
        nonsyn_mask
        & (df["_treatment_only"] | df["_increased"])
        & ~df["_similar_to_control_par"]
        & (df["_max_treatment_dp"] > args.min_dp)
    )

    snp_out = df.loc[keep_mask].copy()
    present_cols = [col for col in snp_out.columns if col.endswith("_PRESENT")]
    if present_cols:
        snp_out = snp_out.drop(columns=present_cols)
    snp_out.to_csv(args.output_snps, sep="\t", index=False)

    gene_cols = ["ANN_GENE_NAME", "ANN_GENE_ID"]
    gene_out = (
        snp_out.groupby(gene_cols, dropna=False)
        .agg(
            variant_count=("POS", "count"),
            max_treatment_af=("_max_treatment_af", "max"),
            mean_control_af=("_mean_control_af", "mean"),
            treatment_only=("_treatment_only", "max"),
            frequency_increased=("_increased", "max"),
            ann_effects=(
                "ANN_EFFECT",
                lambda s: ",".join(sorted({str(x) for x in s.dropna()})),
            ),
        )
        .reset_index()
        .sort_values(["variant_count", "max_treatment_af"], ascending=False)
    )
    treatment_presence = (
        snp_out[treatment_cols]
        .rename(columns=dict(zip(treatment_cols, treatment_names)))
        .gt(thresholds.presence_threshold)
        .apply(lambda row: sorted(row.index[row].tolist()), axis=1)
    )
    treatments_per_gene = (
        treatment_presence.groupby(snp_out[gene_cols].apply(tuple, axis=1))
        .apply(lambda items: ",".join(sorted({name for names in items for name in names})))
    )
    gene_out["treatments_with_mutations"] = (
        gene_out[gene_cols].apply(tuple, axis=1).map(treatments_per_gene).fillna("")
    )
    gene_out.to_csv(args.output_genes, sep="\t", index=False)


if __name__ == "__main__":
    main()
