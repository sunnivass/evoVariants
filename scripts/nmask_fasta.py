#!/usr/bin/env python3
"""
N-mask a reference FASTA using a BED of regions, accepting flexible chromosome identifiers.

Examples
--------
# Basic
python nmask_fasta.py --fasta ref.fa --bed mask.bed --out ref.masked.fa

# With gzip
python nmask_fasta.py --fasta ref.fa.gz --bed mask.bed.gz --out - > ref.masked.fa

# Treat BED coordinates as 1-based closed (if your BED isn't standard)
python nmask_fasta.py --fasta ref.fa --bed mask.tsv --bed-coords 1-based-closed

What it does
------------
- Parses FASTA headers like:
  ">NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence"
- Builds a many-to-one alias map so that BED chrom names like "II", "2", "chrII",
  "chr2", and the RefSeq accession (e.g., "NC_001134.8") all resolve to the
  correct sequence.
- Applies 0-based half-open BED intervals by default (like bcftools consensus -m).
- Writes wrapped FASTA (60 bp/line) unless --no-wrap is given.

Notes
-----
- Aliasing is heuristic but robust for common genomes. It prefers extracting a
  chromosome token from the description (e.g. "chromosome I"). You may also
  provide an explicit alias map via --alias TSV if needed.
- Intervals extending beyond sequence ends are clipped. Overlaps are merged.

Author: you
License: MIT
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple

# ---------------------------- IO helpers ---------------------------- #

def _open_maybe_gzip(path: str, mode: str = "rt"):
    if path == "-":
        # stdin/stdout
        if "r" in mode:
            return io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8")
        else:
            return io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

# ------------------------- FASTA parsing/writing -------------------- #

def read_fasta(path: str) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) for each record. Header excludes leading '>'."""
    with _open_maybe_gzip(path, "rt") as fh:
        header = None
        chunks: List[str] = []
        for line in fh:
            if not line:
                break
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n\r")
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)

def wrap_fasta(seq: str, width: int = 60) -> str:
    if width <= 0:
        return seq + "\n"
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width)) + "\n"

# --------------------------- Alias utilities ------------------------ #

_ROMAN_MAP = {
    "I":1, "V":5, "X":10, "L":50, "C":100, "D":500, "M":1000
}

def roman_to_int(s: str) -> int | None:
    s = s.upper()
    if not s or any(ch not in _ROMAN_MAP for ch in s):
        return None
    total = 0
    prev = 0
    for ch in reversed(s):
        val = _ROMAN_MAP[ch]
        if val < prev:
            total -= val
        else:
            total += val
            prev = val
    return total

_acc_re = re.compile(r"^(?P<acc>[A-Z]{2}_[0-9]+\.[0-9]+)")
_chr_token_re = re.compile(r"chromosome\s+([ivxlcdm]+|\d+)", re.I)
_mito_re = re.compile(r"mito(chondr(ial|ion))?|chr\s*m(t)?\b", re.I)
_plasmid_re = re.compile(r"2[- ]?micron|2μm|2um", re.I)


def normalize_chrom_name(name: str) -> str:
    """Normalize a chromosome-like string to a canonical key for aliasing.

    Strategy:
    - strip leading 'chr' and common punctuation
    - map roman numerals to integer
    - unify mitochondrial names to 'MT'
    - uppercase canonical form
    """
    raw = name.strip()
    n = re.sub(r"^chr", "", raw, flags=re.I)
    n = re.sub(r"[_\s-]", "", n)
    up = n.upper()

    # Mito shortcuts
    if up in {"M", "MT", "MITO", "MITOCHONDRIA", "MITOCHONDRION"}:
        return "MT"

    # Roman numerals or integer
    r = roman_to_int(up)
    if r is not None:
        return f"{r}"
    if up.isdigit():
        return str(int(up))

    # Keep as-is for accessions or other tokens
    return up


def build_aliases(header: str) -> List[str]:
    """Return a list of alias keys for a FASTA header.

    Includes: accession, parsed chromosome tokens, with/without 'chr',
    arabic and roman variants, and special handling for MT and 2-micron.
    """
    aliases = set()
    h = header

    # 1) Accession at start
    m = _acc_re.match(h)
    if m:
        aliases.add(m.group("acc").upper())

    # 2) Chromosome token from description
    m2 = _chr_token_re.search(h)
    if m2:
        tok = m2.group(1)
        norm = normalize_chrom_name(tok)
        if norm.isdigit():
            n = int(norm)
            roman = int_to_roman(n)
            for base in {str(n), roman}:
                aliases.update({
                    base.upper(),
                    ("CHR" + base).upper(),
                })
        else:
            # Uncommon path
            aliases.update({norm, ("CHR" + norm)})

    # 3) Mitochondrial names
    if _mito_re.search(h):
        aliases.update({"MT", "CHRMT"})

    # 4) 2-micron plasmid common in S. cerevisiae
    if _plasmid_re.search(h):
        aliases.update({"2MICRON", "CHR2MICRON"})

    # 5) Also add the full defline token before first space as an alias (common tool behavior)
    first_tok = h.split()[0]
    aliases.add(first_tok.upper())

    return sorted(aliases)


def int_to_roman(num: int) -> str:
    if num <= 0:
        return str(num)
    vals = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
    ]
    res = []
    n = num
    for v, sym in vals:
        while n >= v:
            res.append(sym)
            n -= v
    return "".join(res)

# ------------------------------ BED parsing ------------------------ #

def parse_bed(path: str, coord_mode: str = "0-based-half-open") -> Dict[str, List[Tuple[int,int]]]:
    """Return mapping chrom -> list of (start,end) intervals in 0-based half-open.

    coord_mode:
        - '0-based-half-open' (default, standard BED)
        - '1-based-closed'    (will be converted to 0-based half-open)
    """
    conv = coord_mode
    out: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    with _open_maybe_gzip(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 3:
                sys.stderr.write(f"[WARN] Skipping line {ln}: not enough columns\n")
                continue
            chrom, a, b = parts[:3]
            try:
                start = int(a)
                end = int(b)
            except ValueError:
                sys.stderr.write(f"[WARN] Skipping line {ln}: start/end not integers\n")
                continue

            if conv == "1-based-closed":
                # Convert [start,end] (1-based, inclusive) -> [start-1, end) 0-based half-open
                start = start - 1
                # end remains inclusive, so +1 to make half-open
                end = end
            elif conv == "0-based-half-open":
                pass
            else:
                raise ValueError(f"Unknown coord_mode: {coord_mode}")

            if end < start:
                start, end = end, start

            out[chrom].append((start, end))
    return out


def merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged

# ----------------------------- Masking core ------------------------ #

def apply_masks(seq: str, ivals: List[Tuple[int,int]]) -> str:
    if not ivals:
        return seq
    arr = bytearray(seq.encode("ascii"))
    n = len(arr)
    for s, e in ivals:
        s = max(0, s)
        e = min(n, e)
        if s >= e:
            continue
        arr[s:e] = b"N" * (e - s)
    return arr.decode("ascii")

# ------------------------------ Main logic ------------------------- #

def load_alias_map(fasta_path: str) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """Return (alias->record_id, records)
    records list contains tuples of (header, seq). The record_id = defline token (without '>').
    """
    alias_to_id: Dict[str, str] = {}
    records: List[Tuple[str, str]] = []

    for header, seq in read_fasta(fasta_path):
        record_id = header.split()[0]
        records.append((header, seq))
        aliases = build_aliases(header)
        for a in aliases:
            alias_to_id[a] = record_id

    return alias_to_id, records


def load_user_aliases(path: str | None) -> Dict[str, str]:
    """Optional TSV with two columns: alias\trecord_id"""
    if not path:
        return {}
    mapping: Dict[str, str] = {}
    with _open_maybe_gzip(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                sys.stderr.write(f"[WARN] alias TSV line {ln} has <2 columns\n")
                continue
            alias, rec = parts[0].strip(), parts[1].strip()
            mapping[alias.upper()] = rec
    return mapping


def main():
    p = argparse.ArgumentParser(description="N-mask a FASTA using a BED, with smart chromosome aliasing.")
    p.add_argument("--fasta", required=True, help="Reference FASTA (can be .gz). Use '-' for stdin.")
    p.add_argument("--bed", required=True, help="BED (can be .gz). Use '-' for stdin.")
    p.add_argument("--out", required=True, help="Output FASTA (use '-' for stdout)")
    p.add_argument("--bed-coords", default="0-based-half-open", choices=["0-based-half-open", "1-based-closed"],
                   help="BED coordinate convention; default: standard BED (0-based half-open)")
    p.add_argument("--alias-tsv", default=None, help="Optional TSV: alias\trecord_id to supplement/override inference")
    p.add_argument("--no-wrap", action="store_true", help="Do not wrap FASTA lines (default wraps at 60bp)")

    args = p.parse_args()

    # Load FASTA and alias map
    alias_map, records = load_alias_map(args.fasta)

    # Merge in user-provided aliases (override built ones)
    user_alias = load_user_aliases(args.alias_tsv)
    alias_map.update(user_alias)

    # Build a set of valid IDs for quick membership
    valid_ids = {h.split()[0] for h, _ in records}

    # Parse BED and remap chromosomes via aliasing
    bed_raw = parse_bed(args.bed, coord_mode=args.bed_coords)

    per_id_intervals: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    unknown_chroms = set()

    for chrom, ivals in bed_raw.items():
        # Try direct accession/id match first
        key_direct = chrom.upper()
        rec_id = alias_map.get(key_direct)
        if rec_id is None:
            # Try normalized forms
            norm_key = normalize_chrom_name(chrom)
            # Accept both raw and with CHR prefix as aliases
            rec_id = alias_map.get(norm_key)
            if rec_id is None:
                rec_id = alias_map.get(("CHR" + norm_key))
        if rec_id is None or rec_id not in valid_ids:
            unknown_chroms.add(chrom)
            continue
        per_id_intervals[rec_id].extend(ivals)

    # Merge overlaps per record
    for k in list(per_id_intervals.keys()):
        per_id_intervals[k] = merge_intervals(per_id_intervals[k])

    # Mask and write
    width = 0 if args.no_wrap else 60
    out_fh = _open_maybe_gzip(args.out, "wt")
    with out_fh as oh:
        for header, seq in records:
            rec_id = header.split()[0]
            masked = apply_masks(seq, per_id_intervals.get(rec_id, []))
            oh.write(f">{header}\n")
            oh.write(wrap_fasta(masked, width))

    if unknown_chroms:
        sys.stderr.write("[WARN] BED contained chromosomes not found via aliasing: " + ", ".join(sorted(unknown_chroms)) + "\n")


if __name__ == "__main__":
    main()
