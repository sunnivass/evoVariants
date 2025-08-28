#!/usr/bin/env python3
import sys, argparse, re
from collections import OrderedDict

DNA_ALPHABET = set("ACGTNacgtn")

# --- FASTA I/O ---
def load_fasta(path):
    ref = OrderedDict()
    h_token, h_full, chunks = None, None, []
    fh = sys.stdin if path == "-" else open(path, "r")
    with fh:
        for line in fh:
            if line.startswith(">"):
                if h_token is not None:
                    ref[h_token] = {"full": h_full, "seq": "".join(chunks)}
                h_full = line[1:].rstrip("\n").strip()
                h_token = h_full.split()[0]  # accession/token
                chunks = []
            else:
                chunks.append(line.strip())
        if h_token is not None:
            ref[h_token] = {"full": h_full, "seq": "".join(chunks)}
    if not ref:
        sys.exit("ERROR: no sequences found in FASTA.")
    return ref

def write_fasta(ref, path, width=60):
    out = sys.stdout if path == "-" else open(path, "w")
    with out:
        for token, rec in ref.items():
            out.write(f">{rec['full']}\n")
            s = rec["seq"]
            for i in range(0, len(s), width):
                out.write(s[i:i+width] + "\n")

# --- Utils ---
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

def read_new_seq(literal=None, fasta_path=None):
    if literal is not None:
        return re.sub(r"\s+", "", literal)
    if fasta_path is None:
        sys.exit("ERROR: provide --new-seq or --new-seq-fasta.")
    seq, started = [], False
    fh = sys.stdin if fasta_path == "-" else open(fasta_path, "r")
    with fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            seq.append(line.strip())
    s = "".join(seq)
    if not s:
        sys.exit("ERROR: replacement sequence empty.")
    return s

def replace_slice(seq, start1, end1, insert_seq):
    if start1 < 1 or start1 > end1:
        sys.exit("ERROR: start must be >=1 and <= end.")
    if end1 > len(seq):
        sys.exit(f"ERROR: end ({end1}) exceeds sequence length ({len(seq)}).")
    s0, e0 = start1 - 1, end1
    return seq[:s0] + insert_seq + seq[e0:]

# --- Roman/Arabic helpers ---
_roman_map = [
    (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
    (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
    (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
]
def int_to_roman(n: int) -> str:
    if n <= 0:
        raise ValueError("Roman numerals require positive integers")
    out = []
    for val, sym in _roman_map:
        while n >= val:
            out.append(sym); n -= val
    return "".join(out)

def roman_to_int(s: str):
    s = s.upper()
    vals = {"I":1,"V":5,"X":10,"L":50,"C":100,"D":500,"M":1000}
    total = 0; prev = 0
    for ch in reversed(s):
        if ch not in vals: return None
        v = vals[ch]
        if v < prev: total -= v
        else: total += v; prev = v
    # basic validity check: re-encode and compare
    try:
        if int_to_roman(total) != s: return None
    except Exception:
        return None
    return total

def normalize_numeral_forms(q: str):
    """
    Given 'II', '2', 'chr2', 'chrII', return a set of candidate numerals to try in headers:
      e.g., {'II','2'}
    """
    x = q.strip()
    if x.lower().startswith("chr"):
        x = x[3:].strip()
    candidates = set()
    # raw
    if x:
        candidates.add(x.upper())
    # if arabic int -> add roman
    if re.fullmatch(r"\d+", x):
        try:
            candidates.add(int_to_roman(int(x)))
        except Exception:
            pass
    else:
        # if roman -> add arabic
        val = roman_to_int(x)
        if val is not None:
            candidates.add(str(val))
    return candidates

# --- Flexible contig finder ---
def find_target_contig(query, ref):
    """
    Matching priority:
      1) Exact accession token (e.g. NC_001134.8)
      2) If query includes the word 'chromosome' -> substring match in full header (case-insensitive)
      3) Numeral-aware regex on 'chromosome <NUM>' trying both Arabic and Roman forms
    """
    q = query.strip()
    # 1) exact token
    if q in ref:
        return q

    q_lower = q.lower()

    # 2) substring only if clearly descriptive (avoid tiny queries like '2')
    if "chromosome" in q_lower:
        for token, rec in ref.items():
            if q_lower in rec["full"].lower():
                return token

    # 3) numerals (II / 2 / chrII / chr2)
    numerals = normalize_numeral_forms(q)
    for token, rec in ref.items():
        full = rec["full"]
        for num in numerals:
            # try 'chromosome <NUM>' exact word match
            patt = re.compile(r"\bchromosome\s+" + re.escape(num) + r"\b", re.IGNORECASE)
            if patt.search(full):
                return token
            # also allow headers that literally have 'chrNUM' tokens (common in custom FASTAs)
            patt2 = re.compile(r"\bchr" + re.escape(num) + r"\b", re.IGNORECASE)
            if patt2.search(full):
                return token

    return None

# --- Main ---
def main():
    ap = argparse.ArgumentParser(
        description="Replace a contiguous region in a reference FASTA using contig, start, end. "
                    "Contig may be an accession (e.g., NC_001134.8), 'chromosome II', 'II', '2', 'chrII', or 'chr2'."
    )
    ap.add_argument("--ref", required=True, help="Reference FASTA (use - for stdin)")
    ap.add_argument("--out", required=True, help="Output FASTA (use - for stdout)")
    ap.add_argument("--contig", required=True, help="Accession, 'chromosome II', 'II', '2', 'chrII', or 'chr2'")
    ap.add_argument("--start", type=int, required=True, help="1-based start (inclusive)")
    ap.add_argument("--end", type=int, required=True, help="1-based end (inclusive)")
    ap.add_argument("--strand", choices=["+", "-"], default="+",
                    help="If '-', reverse-complement the provided sequence when --new-seq-is-coding-sense is set (default: +)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--new-seq", help="Replacement sequence literal (ACGTN...)")
    g.add_argument("--new-seq-fasta", help="Single-record FASTA with replacement (use - for stdin)")
    ap.add_argument("--new-seq-is-coding-sense", action="store_true",
                    help="If set and --strand is '-', reverse-complement the provided sequence before insertion.")
    ap.add_argument("--allow-non-dna", action="store_true",
                    help="Skip alphabet check for replacement sequence.")
    args = ap.parse_args()

    ref = load_fasta(args.ref)
    target_token = find_target_contig(args.contig, ref)
    if not target_token:
        examples = [ref[k]["full"] for k in list(ref.keys())[:5]]
        sys.exit(f"ERROR: Could not match contig '{args.contig}' to any FASTA header.\n"
                 f"Examples: {examples}")

    new_seq = read_new_seq(args.new_seq, args.new_seq_fasta)

    if not args.allow_non_dna:
        bad = set(new_seq) - DNA_ALPHABET
        if bad:
            sys.exit(f"ERROR: replacement contains non-DNA characters: {''.join(sorted(bad))}")

    if args.new_seq_is_coding_sense and args.strand == "-":
        new_seq = reverse_complement(new_seq)

    orig = ref[target_token]["seq"]
    replaced = replace_slice(orig, args.start, args.end, new_seq)
    ref[target_token]["seq"] = replaced

    write_fasta(ref, args.out)

    orig_len = args.end - args.start + 1
    delta = len(new_seq) - orig_len
    sys.stderr.write(
        f"Replaced {target_token}:{args.start}-{args.end} (len {orig_len}, strand {args.strand}) "
        f"with {len(new_seq)} bp (Δ {delta:+}). New contig length {len(replaced)} bp.\n"
    )

if __name__ == "__main__":
    main()


#How to use:
#bash
#python replace_region_in_fasta.py \
#  --ref genome.fa \
#  --out edited.fa \
#  --contig chr3 \
#  --start 1250034 \
#  --end 1254567 \
#  --new-seq ACGTACGT...ACGT