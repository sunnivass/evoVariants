#!/usr/bin/env python3
"""
Pathway & protein-class enrichment analysis for evoVariants.

Takes a filtered SnpEff CSV (from filter_snpeff_variants.py) and runs:
  1. GO / KEGG / Reactome ORA via g:Profiler REST API
  2. Custom gene-set enrichment for drug-resistance / proteasome /
     ubiquitin / ergosterol classes (Fisher's exact test)
  3. Publication-ready tables (TSV) and plots (PDF/PNG/SVG)

Usage
-----
  python scripts/enrichment_analysis.py \\
      --input output/PfRpn11_filtered-snpAF0-02.snpeff.genes.csv

  # Restrict to HIGH + MODERATE impact only:
  python scripts/enrichment_analysis.py \\
      --input output/PfRpn11_filtered-snpAF0-02.snpeff.genes.csv \\
      --impacts HIGH MODERATE

  # Offline (skip g:Profiler API, only run curated gene-set tests):
  python scripts/enrichment_analysis.py \\
      --input output/PfRpn11_filtered-snpAF0-02.snpeff.genes.csv \\
      --no-gprofiler
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

try:
    from gprofiler import GProfiler

    HAS_GPROFILER = True
except ImportError:
    HAS_GPROFILER = False

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)


THEME_KEYWORDS: dict[str, list[str]] = {
    "ubiquitin": ["ubiquitin", "deubiquitin", "ubp", "ubr", "rsp5"],
    "proteasome": ["proteasome", "rpn", "rpt", "pre", "pup"],
    "multidrug_resistance": ["drug", "multidrug", "xenobiotic", "efflux", "resistance"],
    "transporters": ["transporter", "abc", "mfs", "vba", "azr", "flr", "qdr"],
    "ergosterol": ["ergosterol", "sterol", "erg"],
    "v_atpase": ["v-?atpase", "vacuolar", "vma", "vph1", "stv1", "proton transmembrane"],
}

MANUAL_SET_THEME: dict[str, str] = {
    "ABC_drug_transporters": "transporters",
    "MFS_drug_transporters": "transporters",
    "Drug_resistance_TFs": "multidrug_resistance",
    "Multidrug_resistance_combined": "multidrug_resistance",
    "Proteasome_26S": "proteasome",
    "Ubiquitin_pathway": "ubiquitin",
    "Ergosterol_biosynthesis": "ergosterol",
    "V-ATPase": "v_atpase",
}

# ═══════════════════════════════════════════════════════════════════════
# Curated S. cerevisiae gene sets
# ═══════════════════════════════════════════════════════════════════════

GENE_SETS: dict[str, dict] = {
    "ABC_drug_transporters": {
        "description": "ABC transporter efflux pumps (pleiotropic drug resistance)",
        "genes": {
            "PDR5", "PDR10", "PDR11", "PDR12", "PDR15", "PDR18",
            "SNQ2", "YOR1", "AUS1", "STE6", "YCF1", "VMR1",
            "MDL1", "MDL2", "ATM1", "YBT1", "NFT1",
        },
    },
    "MFS_drug_transporters": {
        "description": "MFS transporters linked to drug efflux / resistance",
        "genes": {
            "FLR1", "AZR1", "ATR1", "TPO1", "TPO2", "TPO3", "TPO4",
            "QDR1", "QDR2", "QDR3", "DTR1", "SGE1", "HOL1",
            "VBA1", "VBA2", "VBA3", "VBA4", "VBA5",
            "FLC1", "FLC2", "FLC3",
        },
    },
    "Drug_resistance_TFs": {
        "description": "Transcription factors regulating multidrug resistance",
        "genes": {
            "PDR1", "PDR3", "YRR1", "YAP1", "STB5", "RDR1",
            "WAR1", "CIN5", "CAD1", "YAP5", "YAP6", "YAP7", "HAL9",
        },
    },
    "Proteasome_26S": {
        "description": "26S proteasome: 19S RP (base + lid) + 20S CP subunits",
        "genes": {
            # 19S base ATPases
            "RPT1", "RPT2", "RPT3", "RPT4", "RPT5", "RPT6",
            # 19S base non-ATPases
            "RPN1", "RPN2", "RPN13",
            # 19S lid
            "RPN3", "RPN5", "RPN6", "RPN7", "RPN8", "RPN9",
            "RPN11", "RPN12", "SEM1",
            # Ubiquitin receptors / shuttle factors
            "RPN4", "RPN10", "DSK2", "RAD23", "DDI1",
            # 20S core particle
            "PRE1", "PRE2", "PRE3", "PRE4", "PRE5", "PRE6",
            "PRE7", "PRE8", "PRE9", "PRE10",
            "PUP1", "PUP2", "PUP3", "SCL1",
            # Assembly / quality control
            "UMP1", "BLM10", "ECM29", "POC4",
        },
    },
    "Ubiquitin_pathway": {
        "description": "E1 / E2 / E3 ligases, DUBs, ubiquitin genes",
        "genes": {
            # Ubiquitin
            "UBI1", "UBI2", "UBI3", "UBI4",
            # E1
            "UBA1",
            # E2
            "UBC1", "UBC2", "UBC3", "UBC4", "UBC5", "UBC6", "UBC7",
            "UBC8", "UBC10", "UBC11", "UBC13",
            # E3
            "RSP5", "UBR1", "UBR2", "TOM1", "HUL4", "HUL5", "SAN1",
            "ITT1", "BUL1", "BUL2", "RKR1", "DOA10", "HRD1", "UFD4",
            # DUBs
            "UBP1", "UBP2", "UBP3", "UBP5", "UBP6", "UBP7", "UBP8",
            "UBP9", "UBP10", "UBP11", "UBP12", "UBP13", "UBP14",
            "UBP15", "UBP16", "DOA4", "OTU1", "OTU2",
            # Cofactors
            "BRE5", "DOA1", "SHP1", "CDC48", "UFD1", "NPL4",
        },
    },
    "Ergosterol_biosynthesis": {
        "description": "Ergosterol / sterol biosynthetic pathway",
        "genes": {
            "ERG1", "ERG2", "ERG3", "ERG4", "ERG5", "ERG6", "ERG7",
            "ERG8", "ERG9", "ERG10", "ERG11", "ERG12", "ERG13",
            "ERG20", "ERG24", "ERG25", "ERG26", "ERG27", "ERG28",
            "HMG1", "HMG2", "MVD1", "IDI1",
        },
    },
    "V-ATPase": {
        "description": "Vacuolar H+-ATPase subunits",
        "genes": {
            "VMA1", "VMA2", "VMA3", "VMA4", "VMA5", "VMA6", "VMA7",
            "VMA8", "VMA9", "VMA10", "VMA11", "VMA13", "VMA16",
            "VMA21", "VMA22", "VPH1", "STV1",
        },
    },
}

# Aggregate MDR superset
GENE_SETS["Multidrug_resistance_combined"] = {
    "description": "Combined: ABC + MFS transporters + drug-resistance TFs",
    "genes": (
        GENE_SETS["ABC_drug_transporters"]["genes"]
        | GENE_SETS["MFS_drug_transporters"]["genes"]
        | GENE_SETS["Drug_resistance_TFs"]["genes"]
    ),
}


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Pathway & protein-class enrichment for evoVariants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input", required=True,
        help="Filtered SnpEff CSV/TSV from filter_snpeff_variants.py.",
    )
    p.add_argument(
        "--output-dir", default=None,
        help="Output directory (default: <input-parent>/enrichment).",
    )
    p.add_argument(
        "--organism", default="scerevisiae",
        help="g:Profiler organism ID (default: scerevisiae).",
    )
    p.add_argument(
        "--impacts", nargs="*", default=None,
        help=(
            "Keep only variants with these ANN_IMPACT values "
            "(e.g. HIGH MODERATE). Default: use all variants."
        ),
    )
    p.add_argument(
        "--background-size", type=int, default=6275,
        help=(
            "Background gene-universe size for Fisher tests "
            "(default: 6275 = S. cerevisiae verified ORFs)."
        ),
    )
    p.add_argument(
        "--background-genes", default=None,
        help="File with one gene per line for custom g:Profiler background.",
    )
    p.add_argument(
        "--fdr", type=float, default=0.05,
        help="FDR threshold for reporting (default: 0.05).",
    )
    p.add_argument(
        "--gprofiler-user-threshold", type=float, default=1.0,
        help=(
            "g:Profiler user_threshold (default: 1.0 to retrieve broader official "
            "term gene-lists; use --fdr for significance interpretation)."
        ),
    )
    p.add_argument(
        "--max-term-size", type=int, default=500,
        help="Skip broad GO terms larger than this (default: 500).",
    )
    p.add_argument(
        "--no-gprofiler", action="store_true",
        help="Skip g:Profiler API (offline mode).",
    )
    p.add_argument(
        "--plot-format", default="pdf", choices=["pdf", "png", "svg"],
        help="Plot file format (default: pdf).",
    )
    p.add_argument(
        "--dpi", type=int, default=300,
        help="Plot DPI (default: 300).",
    )
    p.add_argument(
        "--focus-keywords",
        nargs="*",
        default=[
            "drug", "multidrug", "efflux", "transporter", "xenobiotic",
            "proteasome", "ubiquitin", "deubiquitin", "ergosterol",
            "sterol", "vacuolar", "ATPase",
        ],
        help=(
            "Keywords for selecting focused g:Profiler terms in an extra report "
            "(default: MDR/proteasome/ubiquitin/ergosterol/vATPase terms)."
        ),
    )
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════════

def load_variants(path: str, impacts: list[str] | None = None) -> pd.DataFrame:
    """Load the filtered SnpEff CSV/TSV, optionally subset by impact."""
    sep = "\t" if path.endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)
    if impacts:
        impacts_upper = {i.upper() for i in impacts}
        df = df[df["ANN_IMPACT"].str.upper().isin(impacts_upper)]
    return df


def extract_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicated gene table with name, ID, variant count, best impact, max AF."""
    impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    df = df.copy()
    df["_rank"] = df["ANN_IMPACT"].map(impact_rank).fillna(4)

    af_col = "_max_treatment_af" if "_max_treatment_af" in df.columns else None

    agg_dict: dict = {
        "gene_name": ("ANN_GENE_NAME", "first"),
        "gene_id": ("ANN_GENE_ID", "first"),
        "n_variants": ("POS", "count"),
        "best_impact": ("ANN_IMPACT", "first"),
        "effects": (
            "ANN_EFFECT",
            lambda s: ";".join(sorted({str(x) for x in s.dropna()})),
        ),
    }
    if af_col:
        agg_dict["max_treatment_af"] = (af_col, "max")

    gene_df = (
        df.sort_values("_rank")
        .groupby("ANN_GENE_ID", as_index=False)
        .agg(**agg_dict)
    )
    return gene_df


def canonical_gene_label(gene_name: str | float, gene_id: str | float) -> str:
    """Choose one canonical label per gene: standard name, else systematic ID."""
    if pd.notna(gene_name):
        name = str(gene_name).strip()
        if name and name.lower() != "nan":
            return name
    if pd.notna(gene_id):
        gid = str(gene_id).strip()
        if gid and gid.lower() != "nan":
            return gid
    return ""


def build_hit_gene_sets(gene_df: pd.DataFrame) -> tuple[list[str], set[str], set[str]]:
    """
    Build three gene collections:
      - query_list: ordered canonical list for g:Profiler
      - canonical_set: unique canonical gene labels for Fisher tests
      - alias_set: unique union of gene names and IDs (for diagnostics)
    """
    query_list: list[str] = []
    canonical_set: set[str] = set()
    alias_set: set[str] = set()

    for _, row in gene_df.iterrows():
        canonical = canonical_gene_label(row.get("gene_name"), row.get("gene_id"))
        if canonical:
            query_list.append(canonical)
            canonical_set.add(canonical)

        gene_name = row.get("gene_name")
        gene_id = row.get("gene_id")
        if pd.notna(gene_name) and str(gene_name).strip():
            alias_set.add(str(gene_name).strip())
        if pd.notna(gene_id) and str(gene_id).strip():
            alias_set.add(str(gene_id).strip())

    query_list = list(OrderedDict.fromkeys(query_list))
    return query_list, canonical_set, alias_set


# ═══════════════════════════════════════════════════════════════════════
# g:Profiler ORA
# ═══════════════════════════════════════════════════════════════════════

def run_gprofiler(
    gene_list: list[str],
    organism: str = "scerevisiae",
    background: list[str] | None = None,
    max_term_size: int = 500,
    user_threshold: float = 1.0,
) -> pd.DataFrame:
    """Query the g:Profiler REST API for GO / KEGG / Reactome / WP ORA."""
    if not HAS_GPROFILER:
        print(
            "[WARN] gprofiler-official not installed — "
            "install with:  pip install gprofiler-official\n"
            "       Skipping g:Profiler enrichment.",
            file=sys.stderr,
        )
        return pd.DataFrame()

    gp = GProfiler(return_dataframe=True)
    try:
        result = gp.profile(
            organism=organism,
            query=gene_list,
            sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
            user_threshold=user_threshold,
            significance_threshold_method="fdr",
            no_evidences=False,
            background=background,
        )
    except Exception as exc:
        print(f"[WARN] g:Profiler API call failed: {exc}", file=sys.stderr)
        return pd.DataFrame()

    if result.empty:
        return result

    result = result[result["term_size"] <= max_term_size].copy()

    # Derived columns useful for plotting
    result["gene_ratio"] = result["intersection_size"] / result["query_size"]
    result["enrichment_ratio"] = (
        (result["intersection_size"] / result["query_size"])
        / (result["term_size"] / result["effective_domain_size"])
    )
    result["neg_log10_padj"] = -np.log10(result["p_value"].clip(lower=1e-300))

    return result.sort_values("p_value").reset_index(drop=True)


def build_gprofiler_focus_report(
    gp_df: pd.DataFrame,
    keywords: list[str],
) -> pd.DataFrame:
    """Select terms whose official names match focus keywords and return compact table."""
    if gp_df.empty:
        return pd.DataFrame()
    if not keywords:
        return pd.DataFrame()

    kw = [k.strip().lower() for k in keywords if k.strip()]
    if not kw:
        return pd.DataFrame()

    pattern = re.compile("|".join(re.escape(k) for k in kw), re.IGNORECASE)
    text = gp_df["name"].fillna("") + " | " + gp_df["native"].fillna("")
    mask = text.str.contains(pattern, na=False)
    focus = gp_df.loc[mask].copy()
    if focus.empty:
        return focus

    def matched_keywords(term: str) -> str:
        found = sorted({k for k in kw if re.search(re.escape(k), term, re.IGNORECASE)})
        return ",".join(found)

    focus["matched_keywords"] = focus["name"].fillna("").map(matched_keywords)
    cols = [
        "source",
        "native",
        "name",
        "matched_keywords",
        "p_value",
        "term_size",
        "query_size",
        "intersection_size",
        "effective_domain_size",
        "enrichment_ratio",
        "intersection",
    ]
    cols = [c for c in cols if c in focus.columns]
    return focus.sort_values("p_value")[cols].reset_index(drop=True)


def _theme_from_text(text: str) -> str:
    low = text.lower()
    for theme, kws in THEME_KEYWORDS.items():
        for kw in kws:
            if re.search(kw, low):
                return theme
    return "other"


def build_consensus_evidence(
    fisher_df: pd.DataFrame,
    focus_df: pd.DataFrame,
) -> pd.DataFrame:
    """Combine manual Fisher sets and g:Profiler focus terms into one evidence table."""
    rows: list[dict] = []

    if not fisher_df.empty:
        for _, row in fisher_df.iterrows():
            theme = MANUAL_SET_THEME.get(str(row["gene_set"]), "other")
            rows.append(
                {
                    "theme": theme,
                    "evidence_source": "manual_geneset_fisher",
                    "entity_id": row["gene_set"],
                    "entity_name": row["description"],
                    "p_value": float(row["pvalue"]),
                    "fdr": float(row["fdr"]),
                    "effect_size": float(row["odds_ratio"]),
                    "signal_size": f"{int(row['overlap'])}/{int(row['set_size'])}",
                    "supporting_genes": row["overlapping_genes"],
                }
            )

    if not focus_df.empty:
        for _, row in focus_df.iterrows():
            term_text = f"{row.get('name', '')} {row.get('matched_keywords', '')}"
            theme = _theme_from_text(term_text)
            overlap = row.get("intersection_size", "")
            term_size = row.get("term_size", "")
            rows.append(
                {
                    "theme": theme,
                    "evidence_source": "gprofiler_term",
                    "entity_id": row.get("native", ""),
                    "entity_name": row.get("name", ""),
                    "p_value": float(row.get("p_value", np.nan)),
                    "fdr": float(row.get("p_value", np.nan)),
                    "effect_size": float(row.get("enrichment_ratio", np.nan)),
                    "signal_size": f"{overlap}/{term_size}",
                    "supporting_genes": row.get("intersection", ""),
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                "theme",
                "evidence_source",
                "entity_id",
                "entity_name",
                "p_value",
                "fdr",
                "effect_size",
                "signal_size",
                "supporting_genes",
            ]
        )

    out = pd.DataFrame(rows)
    out["neg_log10_p"] = -np.log10(out["p_value"].clip(lower=1e-300))
    out = out.sort_values(["theme", "evidence_source", "p_value"]).reset_index(drop=True)
    return out


# ═══════════════════════════════════════════════════════════════════════
# Custom Fisher exact enrichment
# ═══════════════════════════════════════════════════════════════════════

def _fisher_one(
    hit_genes: set[str],
    geneset_genes: set[str],
    background_size: int,
) -> dict:
    """One-tailed Fisher exact test for over-representation."""
    a = len(hit_genes & geneset_genes)              # hit ∩ set
    b = len(hit_genes - geneset_genes)               # hit only
    c = len(geneset_genes - hit_genes)               # set only
    d = max(0, background_size - len(hit_genes | geneset_genes))  # neither
    odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "overlap": a,
        "hit_size": len(hit_genes),
        "set_size": len(geneset_genes),
        "background": background_size,
        "odds_ratio": odds_ratio,
        "pvalue": pvalue,
        "overlapping_genes": ", ".join(sorted(hit_genes & geneset_genes)),
    }


def run_fisher_enrichment(
    hit_gene_names: set[str],
    gene_sets: dict,
    background_size: int,
) -> pd.DataFrame:
    """Run Fisher's exact test for every curated gene set; BH-correct."""
    rows = []
    for name, info in gene_sets.items():
        result = _fisher_one(hit_gene_names, info["genes"], background_size)
        result["gene_set"] = name
        result["description"] = info["description"]
        rows.append(result)

    df = pd.DataFrame(rows).sort_values("pvalue").reset_index(drop=True)

    # Benjamini-Hochberg FDR
    n = len(df)
    raw_fdr = df["pvalue"].values * n / (np.arange(n) + 1)
    # enforce monotonicity (step-up)
    for i in range(n - 2, -1, -1):
        raw_fdr[i] = min(raw_fdr[i], raw_fdr[i + 1])
    df["fdr"] = np.clip(raw_fdr, 0, 1)

    cols = [
        "gene_set", "description", "overlap", "set_size", "hit_size",
        "background", "odds_ratio", "pvalue", "fdr", "overlapping_genes",
    ]
    return df[cols]


# ═══════════════════════════════════════════════════════════════════════
# Plotting helpers
# ═══════════════════════════════════════════════════════════════════════

_SOURCE_COLORS = {
    "GO:BP": "#E64B35", "GO:MF": "#4DBBD5", "GO:CC": "#00A087",
    "KEGG": "#3C5488", "REAC": "#F39B7F", "WP": "#8491B4",
}


def _pub_style() -> None:
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 150,
    })


# ── g:Profiler dot plot ──────────────────────────────────────────────

def plot_gprofiler_dotplot(
    df: pd.DataFrame, outpath: Path, top_n: int = 25, dpi: int = 300,
) -> None:
    if df.empty:
        return
    _pub_style()
    plot_df = df.head(top_n).copy()
    plot_df["label"] = plot_df["name"].str[:60]
    plot_df = plot_df.sort_values("gene_ratio")

    fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(plot_df))))
    sc = ax.scatter(
        plot_df["gene_ratio"],
        range(len(plot_df)),
        s=plot_df["intersection_size"] * 30 + 20,
        c=plot_df["neg_log10_padj"],
        cmap="viridis_r",
        edgecolors="black", linewidths=0.5, zorder=3,
    )
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["label"])
    ax.set_xlabel("Gene ratio (hits / query)")
    ax.set_title("g:Profiler ORA – enriched terms")
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label(r"$-\log_{10}$(adj. p)")

    # size legend
    for sz in (1, 3, 5, 10):
        ax.scatter([], [], s=sz * 30 + 20, c="grey",
                   edgecolors="black", linewidths=0.5, label=str(sz))
    ax.legend(title="Gene count", loc="lower right", frameon=True, fontsize=8)

    plt.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {outpath}")


# ── g:Profiler bar plot ──────────────────────────────────────────────

def plot_gprofiler_barplot(
    df: pd.DataFrame, outpath: Path, top_n: int = 20, dpi: int = 300,
) -> None:
    if df.empty:
        return
    _pub_style()
    plot_df = df.head(top_n).copy()
    plot_df["label"] = plot_df["name"].str[:55]
    plot_df = plot_df.sort_values("neg_log10_padj")

    colors = [_SOURCE_COLORS.get(s, "#999999") for s in plot_df["source"]]

    fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(plot_df))))
    ax.barh(range(len(plot_df)), plot_df["neg_log10_padj"],
            color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["label"])
    ax.set_xlabel(r"$-\log_{10}$(adj. p-value)")
    ax.set_title("Top enriched pathways / GO terms")

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, edgecolor="black", label=s)
               for s, c in _SOURCE_COLORS.items()
               if s in plot_df["source"].values]
    ax.legend(handles=handles, loc="lower right", fontsize=8, frameon=True)

    plt.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {outpath}")


# ── Custom gene-set bar plot ────────────────────────────────────────

def plot_geneset_barplot(
    df: pd.DataFrame, outpath: Path, dpi: int = 300,
) -> None:
    if df.empty:
        return
    _pub_style()
    plot_df = df.copy()
    plot_df["neg_log10_p"] = -np.log10(plot_df["pvalue"].clip(lower=1e-300))
    plot_df = plot_df.sort_values("neg_log10_p")

    colors = [
        "#E64B35" if fdr < 0.05 else "#4DBBD5" if fdr < 0.1 else "#AAAAAA"
        for fdr in plot_df["fdr"]
    ]
    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.5 * len(plot_df))))
    ax.barh(range(len(plot_df)), plot_df["neg_log10_p"],
            color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(plot_df)))
    labels = [
        f"{row.gene_set}  ({row.overlap}/{row.set_size})"
        for _, row in plot_df.iterrows()
    ]
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$-\log_{10}$(p-value)")
    ax.set_title("Curated gene-set enrichment (Fisher's exact test)")
    ax.axvline(-np.log10(0.05), color="grey", ls="--", lw=0.8)

    from matplotlib.patches import Patch
    legend_el = [
        Patch(facecolor="#E64B35", edgecolor="black", label="FDR < 0.05"),
        Patch(facecolor="#4DBBD5", edgecolor="black", label="FDR < 0.10"),
        Patch(facecolor="#AAAAAA", edgecolor="black", label="n.s."),
    ]
    ax.legend(handles=legend_el, loc="lower right", fontsize=8, frameon=True)

    plt.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {outpath}")


# ── Resistance-gene heatmap ─────────────────────────────────────────

def plot_variant_heatmap(
    variant_df: pd.DataFrame,
    gene_sets: dict,
    outpath: Path,
    dpi: int = 300,
) -> None:
    """Heatmap of allele frequencies for curated-set genes across treatments."""
    af_cols = [c for c in variant_df.columns if c.endswith(".filtered_AF")]
    control_names = {"control1", "control2", "control3", "control-par"}
    treat_af = [c for c in af_cols
                if c.replace(".filtered_AF", "") not in control_names]
    if not treat_af:
        return

    all_set_genes: set[str] = set()
    for name, info in gene_sets.items():
        if name == "Multidrug_resistance_combined":
            continue
        all_set_genes |= info["genes"]

    mask = (
        variant_df["ANN_GENE_NAME"].isin(all_set_genes)
        | variant_df["ANN_GENE_ID"].isin(all_set_genes)
    )
    subset = variant_df[mask].copy()
    if subset.empty:
        return

    subset["gene_label"] = subset["ANN_GENE_NAME"]
    heat = subset.groupby("gene_label")[treat_af].max().fillna(0)
    heat.columns = [c.replace(".filtered_AF", "") for c in heat.columns]

    # Sort by gene-set membership
    set_of = {}
    for name, info in gene_sets.items():
        if name == "Multidrug_resistance_combined":
            continue
        for g in info["genes"]:
            if g in heat.index and g not in set_of:
                set_of[g] = name
    order = sorted(heat.index, key=lambda g: (set_of.get(g, "ZZZ"), g))
    heat = heat.loc[order]
    if heat.shape[0] == 0:
        return

    _pub_style()
    fig_h = max(3, 0.4 * heat.shape[0])
    fig_w = max(5, 1.2 * heat.shape[1] + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        heat, ax=ax, cmap="YlOrRd", linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Allele frequency", "shrink": 0.6},
        vmin=0, vmax=max(0.2, heat.values.max()),
    )

    # Add gene-set color strip on left
    set_color_map = {
        "ABC_drug_transporters": "#E64B35",
        "MFS_drug_transporters": "#4DBBD5",
        "Drug_resistance_TFs": "#00A087",
        "Proteasome_26S": "#3C5488",
        "Ubiquitin_pathway": "#F39B7F",
        "Ergosterol_biosynthesis": "#8491B4",
        "V-ATPase": "#91D1C2",
    }
    for i, gene in enumerate(order):
        c = set_color_map.get(set_of.get(gene, ""), "#CCCCCC")
        ax.add_patch(plt.Rectangle((-0.6, i), 0.5, 1, color=c,
                                   clip_on=False, transform=ax.transData))

    ax.set_title("Variant AF in curated gene sets across treatments")
    ax.set_ylabel("")
    ax.set_xlabel("")

    # Gene-set legend
    from matplotlib.patches import Patch
    used_sets = sorted({set_of[g] for g in order if g in set_of})
    handles = [Patch(color=set_color_map.get(s, "#CCC"), label=s.replace("_", " "))
               for s in used_sets]
    if handles:
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 0), loc="lower left",
                  fontsize=7, frameon=True, title="Gene set")

    plt.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {outpath}")


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main() -> None:
    args = parse_args()
    inpath = Path(args.input)
    if not inpath.exists():
        sys.exit(f"Error: input not found: {inpath}")

    outdir = Path(args.output_dir) if args.output_dir else inpath.parent / "enrichment"
    outdir.mkdir(parents=True, exist_ok=True)
    fmt = args.plot_format

    print(f"Input : {inpath}")
    print(f"Output: {outdir}/")
    if args.impacts:
        print(f"Impact filter: {', '.join(args.impacts)}")
    print()

    # ── 1. Load data ──
    print("═══ Loading variants ═══")
    variant_df = load_variants(str(inpath), impacts=args.impacts)
    gene_df = extract_genes(variant_df)
    gene_query, canonical_hit_genes, alias_gene_set = build_hit_gene_sets(gene_df)

    print(f"  Variants loaded : {len(variant_df)}")
    print(f"  Unique genes    : {len(gene_df)}")
    print(f"  Canonical hit genes (for Fisher): {len(canonical_hit_genes)}")
    print(f"  Alias labels (name+ID union): {len(alias_gene_set)}")
    for imp, cnt in variant_df["ANN_IMPACT"].value_counts().items():
        print(f"    {imp}: {cnt}")

    # Quick scan – flag curated-set genes already present
    print()
    print("  ── Quick scan: curated-set genes in hit list ──")
    for set_name, info in GENE_SETS.items():
        overlap = canonical_hit_genes & info["genes"]
        if overlap:
            print(f"    {set_name}: {', '.join(sorted(overlap))}")

    # Save gene summary
    gene_summary_path = outdir / "hit_genes_summary.tsv"
    gene_df.to_csv(gene_summary_path, sep="\t", index=False)
    print(f"\n  Saved: {gene_summary_path}")

    # ── 2. g:Profiler ORA ──
    gp_df = pd.DataFrame()
    focus_df = pd.DataFrame()
    if not args.no_gprofiler:
        print("\n═══ g:Profiler ORA ═══")
        background = None
        if args.background_genes:
            bg_path = Path(args.background_genes)
            if bg_path.exists():
                background = [
                    line.strip() for line in bg_path.read_text().splitlines()
                    if line.strip()
                ]
                print(f"  Custom background: {len(background)} genes")

        gp_df = run_gprofiler(
            gene_query,
            organism=args.organism,
            background=background,
            max_term_size=args.max_term_size,
            user_threshold=args.gprofiler_user_threshold,
        )
        if not gp_df.empty:
            gp_path = outdir / "gprofiler_results.tsv"
            gp_df.to_csv(gp_path, sep="\t", index=False)
            sig = gp_df[gp_df["p_value"] < args.fdr]
            print(f"  Total enriched terms : {len(gp_df)}")
            print(f"  Significant (FDR<{args.fdr}): {len(sig)}")
            for src in sorted(sig["source"].unique()):
                print(f"    {src}: {(sig['source'] == src).sum()}")
            print(f"  Saved: {gp_path}")

            focus_df = build_gprofiler_focus_report(gp_df, args.focus_keywords)
            focus_path = outdir / "gprofiler_focus_terms.tsv"
            focus_df.to_csv(focus_path, sep="\t", index=False)
            print(
                "  Focus report (official g:Profiler term gene-lists/intersections): "
                f"{len(focus_df)} terms"
            )
            print(f"  Saved: {focus_path}")
        else:
            print("  No enriched terms returned.")
    else:
        print("\n═══ Skipping g:Profiler (--no-gprofiler) ═══")

    # ── 3. Custom gene-set enrichment (Fisher) ──
    print("\n═══ Custom gene-set enrichment (Fisher's exact) ═══")
    print(f"  Hit genes      : {len(canonical_hit_genes)} (canonical labels)")
    print(f"  Background size: {args.background_size}")
    fisher_df = run_fisher_enrichment(
        canonical_hit_genes, GENE_SETS, args.background_size,
    )
    fisher_path = outdir / "geneset_enrichment.tsv"
    fisher_df.to_csv(fisher_path, sep="\t", index=False)
    print()
    for _, row in fisher_df.iterrows():
        stars = (
            "***" if row.fdr < 0.001 else
            "** " if row.fdr < 0.01 else
            "*  " if row.fdr < 0.05 else
            "   "
        )
        print(
            f"  {stars} {row.gene_set:32s}  "
            f"overlap={row.overlap:2d}/{row.set_size:2d}  "
            f"OR={row.odds_ratio:7.2f}  "
            f"p={row.pvalue:.2e}  FDR={row.fdr:.2e}"
        )
        if row.overlap > 0:
            print(f"        genes: {row.overlapping_genes}")
    print(f"\n  Saved: {fisher_path}")

    consensus_df = build_consensus_evidence(fisher_df, focus_df)
    consensus_path = outdir / "consensus_evidence.tsv"
    consensus_df.to_csv(consensus_path, sep="\t", index=False)
    print(f"  Saved: {consensus_path}")

    if not consensus_df.empty:
        print("  Top consensus signals by theme:")
        for theme, block in consensus_df.groupby("theme", sort=True):
            best = block.sort_values("p_value").iloc[0]
            print(
                f"    {theme:20s} {best.evidence_source:22s} "
                f"p={best.p_value:.2e}  effect={best.effect_size:.2f}"
            )

    # ── 4. Plots ──
    print("\n═══ Generating plots ═══")
    if not gp_df.empty:
        plot_gprofiler_dotplot(gp_df, outdir / f"gprofiler_dotplot.{fmt}", dpi=args.dpi)
        plot_gprofiler_barplot(gp_df, outdir / f"gprofiler_barplot.{fmt}", dpi=args.dpi)
    plot_geneset_barplot(fisher_df, outdir / f"geneset_enrichment.{fmt}", dpi=args.dpi)
    plot_variant_heatmap(
        variant_df, GENE_SETS, outdir / f"resistance_heatmap.{fmt}", dpi=args.dpi,
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
