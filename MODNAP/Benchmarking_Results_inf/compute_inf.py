#!/usr/bin/env python3
"""
compute_inf.py
==============
Compute INF_all and INF_WC between AlphaFold3 (.cif) and experimental (.pdb)
RNA/DNA structures.

Directory layout expected
-------------------------
MODNAP_AF3_models/
  <CCD>/
    <PDB>/
      af_output/
        <subfolder>/
          *.cif          ← AlphaFold3 prediction (first .cif found is used)

MODNAP/
  <CCD>/
    <PDB>/
      *.pdb              ← experimental reference (first .pdb found is used)

INF formula (Parisien et al.)
-----------------------------
  INF = sqrt(PPV * sensitivity)
  PPV         = TP / (TP + FP)
  sensitivity = TP / (TP + FN)
  → INF = sqrt( TP / (TP + FP) * TP / (TP + FN) )

Base-pair annotation
--------------------
Uses DSSR (x3dna-dssr) when available; falls back to a pure-Python
geometric annotator based on biopython (slower, but no extra binary needed).

Usage
-----
python compute_inf.py --af3_root  ../MODNAP_AF3_models --exp_root  ../MODNAP  --output    results_inf.csv --verbose


      [--tool     dssr|biopython]   # default: auto-detect dssr, else biopython
      [--dssr_bin /path/to/x3dna-dssr]
      [--verbose]
"""

import argparse
import csv
import glob
import json
import pandas as pd
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# load CCD mapping

def load_ccd_mapping(xlsx_file):
    df = pd.read_excel(xlsx_file)

    mapping = {}

    for _, row in df.iterrows():
        ccd = str(row["CCD_ID"]).strip().upper()
        base = str(row["Nucleobase"]).strip().upper()

        if base in {"A", "G", "C", "U", "T"}:
            mapping[ccd] = base

    return mapping

CCD_TO_PARENT = load_ccd_mapping(
    "../../list-of-modified-nucleotides-in-the-PDB/all_nucleotide_CCDs_PDB.xlsx"
)

# ===========================================================================
# 1.  DSSR-based base-pair extraction
# ===========================================================================

def run_dssr(structure_path: str, dssr_bin: str = "x3dna-dssr") -> dict:
    """
    Run DSSR on a PDB or CIF file and return parsed JSON output.
    Raises RuntimeError if DSSR is not found or fails.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [dssr_bin, f"-i={structure_path}", "--json", "--quiet",
               f"--prefix={os.path.join(tmpdir, 'out')}"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        json_path = os.path.join(tmpdir, "out.json")
        if not os.path.exists(json_path):
            raise RuntimeError(
                f"DSSR did not produce JSON for {structure_path}.\n"
                f"stderr: {result.stderr[:500]}"
            )
        with open(json_path) as fh:
            return json.load(fh)


def dssr_base_pairs(structure_path: str, dssr_bin: str = "x3dna-dssr"):
    """
    Return two sets of (chain, resi, ins_code, chain, resi, ins_code) tuples:
      wc_pairs   – Watson-Crick / canonical pairs
      all_pairs  – all base pairs (WC + non-canonical + wobble)

    Each pair is stored as a *sorted* 2-tuple of residue keys so that
    (A, B) == (B, A).
    """
    data = run_dssr(structure_path, dssr_bin)

    wc_pairs  = set()
    all_pairs = set()

    pairs_section = data.get("pairs", []) or []
    for bp in pairs_section:
        # nt1 / nt2 have format like "A.G5" or "A.5"
        nt1 = bp.get("nt1", "")
        nt2 = bp.get("nt2", "")
        if not nt1 or not nt2:
            continue
        key = tuple(sorted([_dssr_nt_key(nt1), _dssr_nt_key(nt2)]))
        all_pairs.add(key)

        lw = bp.get("LW", "") or ""
        # WC / canonical: cWW (cis Watson-Crick/Watson-Crick)
        if "WC" in bp.get("name", "") or lw.upper() == "CWW":
            wc_pairs.add(key)

    return wc_pairs, all_pairs


def _dssr_nt_key(nt_id: str):
    """Convert DSSR nucleotide id (e.g. 'A.G14') to a comparable key."""
    return nt_id.strip()


# ===========================================================================
# 2.  Pure-Python / biopython fallback annotator
# ===========================================================================

def biopython_base_pairs(structure_path: str):
    """
    Detect Watson-Crick base pairs from 3D coordinates using biopython.

    Strategy: for each nucleotide pair, check whether the N1/N9-to-N1/N9
    distance and C1'-to-C1' distance are in the canonical range.  This is
    a fast geometric heuristic — good enough for INF calculation.

    Returns (wc_pairs, all_pairs).  all_pairs == wc_pairs here because the
    simple geometric criterion only covers Watson-Crick interactions.
    To detect non-canonical pairs use DSSR.
    """
    try:
        from Bio import PDB as bpdb
        from Bio.PDB import MMCIFParser, PDBParser
    except ImportError:
        raise RuntimeError(
            "biopython is not installed.  "
            "Run:  pip install biopython  or use --tool dssr"
        )

    ext = Path(structure_path).suffix.lower()
    if ext == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("s", structure_path)
    # Use first model only
    model = next(structure.get_models())

    nucleotides = []  # list of (chain_id, res_id, atoms_dict)
    for chain in model:
        for residue in chain:
            res_name = residue.get_resname().strip().upper()
            if res_name not in RNA_DNA_RESIDUES:
                continue
            atoms = {a.get_name().strip(): a.get_vector()
                     for a in residue.get_atoms()}
            nucleotides.append((chain.id, residue.get_id(), res_name, atoms))

    wc_pairs  = set()
    all_pairs = set()

    for i in range(len(nucleotides)):
        ci, ri, ni, ai = nucleotides[i]
        for j in range(i + 1, len(nucleotides)):
            cj, rj, nj, aj = nucleotides[j]

            # C1'–C1' distance filter (fast reject)
            c1i = ai.get("C1'") or ai.get("C1*")
            c1j = aj.get("C1'") or aj.get("C1*")
            if c1i is None or c1j is None:
                continue
            d_c1 = (c1i - c1j).norm()
            if not (8.0 < d_c1 < 13.5):
                continue

            # Glycosidic N distance
            gi = ai.get("N9") or ai.get("N1") or ai.get("N9*") or ai.get("N1*")
            gj = aj.get("N9") or aj.get("N1") or aj.get("N9*") or aj.get("N1*")
            if gi is None or gj is None:
                continue
            d_g = (gi - gj).norm()
            if not (7.5 < d_g < 13.5):
                continue

            pair = tuple(sorted([_bp_key(ci, ri), _bp_key(cj, rj)]))

            # Canonical WC check: complementary bases
            if _is_wc(ni, nj):
                wc_pairs.add(pair)
                all_pairs.add(pair)
            else:
                # Could still be non-canonical; keep in all_pairs
                all_pairs.add(pair)

    return wc_pairs, all_pairs


RNA_DNA_RESIDUES = {
    "A", "U", "G", "C",        # RNA
    "DA", "DT", "DG", "DC",    # DNA
    "ADE", "URA", "GUA", "CYT",
    "rA", "rU", "rG", "rC",
}

WC_PAIRS = {
    frozenset({"A", "U"}),
    frozenset({"G", "C"}),
    frozenset({"DA", "DT"}),
    frozenset({"DG", "DC"}),
    frozenset({"A",  "T"}),     # sometimes DNA residues are named without D
    frozenset({"G",  "C"}),
    frozenset({"G",  "U"}),     # wobble — included in WC for INF_WC
}


def _is_wc(r1: str, r2: str) -> bool:
    return frozenset({r1, r2}) in WC_PAIRS


def _bp_key(chain_id, res_id):
    """Create a hashable residue key from biopython res_id tuple."""
    # res_id = (hetflag, seqnum, icode)
    return (chain_id, res_id[1], res_id[2].strip())


# ===========================================================================
# 3.  Sequence alignment-based residue mapping
# ===========================================================================

def extract_sequence(structure_path: str) -> list:
    """
    Extract one-letter nucleotide sequence (with residue keys) from a
    PDB or CIF file.
    Returns list of (key, one_letter) sorted by chain + seqnum.
    """
    try:
        from Bio import PDB as bpdb
        from Bio.PDB import MMCIFParser, PDBParser
        from Bio.SeqUtils import seq1
    except ImportError:
        raise RuntimeError("biopython required")

    ext = Path(structure_path).suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("s", structure_path)
    model = next(structure.get_models())
    
    records = []

    for chain in model:
        for residue in chain:

            rname = residue.get_resname().strip().upper()

            if rname in {"A","U","G","C","DA","DT","DG","DC"}:
                base = NUCLEOTIDE_3TO1.get(rname, None)
            elif rname in CCD_TO_PARENT:
                base = CCD_TO_PARENT[rname]
                #print(rname, "->", base)
            else:
                continue

            key = _bp_key(chain.id, residue.get_id())

            records.append(
                (
                    chain.id,
                    residue.get_id()[1],
                    key,
                    base
                )
            )
    
    print(structure_path, "nucleotides:", len(records))
    
    records.sort(key=lambda x: (x[0], x[1]))
    return [(r[2], r[3]) for r in records]


NUCLEOTIDE_3TO1 = {
    "A": "A", "U": "U", "G": "G", "C": "C",
    "DA": "A", "DT": "T", "DG": "G", "DC": "C",
    "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
}


def align_sequences(seq_exp: list, seq_af3: list):
    """
    Pairwise global alignment of two nucleotide sequences.
    Returns two lists of (key | None) of equal length representing the
    aligned positions.

    Uses a simple Needleman-Wunsch implementation to avoid extra dependencies.
    """
    letters_exp = [s[1] for s in seq_exp]
    letters_af3 = [s[1] for s in seq_af3]

    keys_exp = [s[0] for s in seq_exp]
    keys_af3 = [s[0] for s in seq_af3]

    M, N = len(letters_exp), len(letters_af3)
    GAP   = -2
    MATCH = 1
    MISM  = -1

    # Score matrix
    dp = [[0] * (N + 1) for _ in range(M + 1)]
    for i in range(M + 1):
        dp[i][0] = i * GAP
    for j in range(N + 1):
        dp[0][j] = j * GAP

    for i in range(1, M + 1):
        for j in range(1, N + 1):
            match = MATCH if letters_exp[i-1] == letters_af3[j-1] else MISM
            dp[i][j] = max(
                dp[i-1][j-1] + match,
                dp[i-1][j]   + GAP,
                dp[i][j-1]   + GAP,
            )

    # Traceback
    aligned_exp = []
    aligned_af3 = []
    i, j = M, N
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            match = MATCH if letters_exp[i-1] == letters_af3[j-1] else MISM
            if dp[i][j] == dp[i-1][j-1] + match:
                aligned_exp.append(keys_exp[i-1])
                aligned_af3.append(keys_af3[j-1])
                i -= 1; j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + GAP:
            aligned_exp.append(keys_exp[i-1])
            aligned_af3.append(None)
            i -= 1
        else:
            aligned_exp.append(None)
            aligned_af3.append(keys_af3[j-1])
            j -= 1

    aligned_exp.reverse()
    aligned_af3.reverse()
    return aligned_exp, aligned_af3


def build_key_mapping(aligned_exp, aligned_af3):
    """
    Return dict: af3_key → exp_key for aligned (non-gap) positions.
    """
    mapping = {}
    for ke, ka in zip(aligned_exp, aligned_af3):
        if ke is not None and ka is not None:
            mapping[ka] = ke
    return mapping


def translate_pairs(pairs: set, key_map: dict) -> set:
    """
    Re-key AF3 base-pair set into experimental residue keys.
    Pairs where either residue has no mapping are dropped.
    """
    translated = set()
    for (ka1, ka2) in pairs:
        ke1 = key_map.get(ka1)
        ke2 = key_map.get(ka2)
        if ke1 is not None and ke2 is not None:
            translated.add(tuple(sorted([ke1, ke2])))
    return translated


# ===========================================================================
# 4.  INF calculation
# ===========================================================================

def inf_score(ref_pairs: set, pred_pairs: set) -> float:
    """
    INF = sqrt(PPV * sensitivity)
    Returns 0.0 when ref or pred is empty.
    """
    if not ref_pairs or not pred_pairs:
        return 0.0
    tp = len(ref_pairs & pred_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)
    if tp == 0:
        return 0.0
    ppv  = tp / (tp + fp)
    sens = tp / (tp + fn)
    return math.sqrt(ppv * sens)


# ===========================================================================
# 5.  File discovery
# ===========================================================================

def find_af3_cif(af3_root: str, ccd: str, pdb: str) -> str | None:
    pattern = os.path.join(af3_root, ccd, pdb, "af_output", "**", "*.cif")
    hits = sorted(glob.glob(pattern, recursive=True))
    return hits[0] if hits else None


def find_exp_pdb(exp_root: str, ccd: str, pdb: str) -> str | None:
    pattern = os.path.join(exp_root, ccd, pdb, "*.pdb")
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def iter_pairs(af3_root: str, exp_root: str):
    """
    Yield (ccd, pdb, af3_path, exp_path) for every matched pair.
    """
    for ccd_dir in sorted(Path(af3_root).iterdir()):
        if not ccd_dir.is_dir():
            continue
        ccd = ccd_dir.name
        for pdb_dir in sorted(ccd_dir.iterdir()):
            if not pdb_dir.is_dir():
                continue
            pdb = pdb_dir.name
            af3_path = find_af3_cif(af3_root, ccd, pdb)
            exp_path = find_exp_pdb(exp_root, ccd, pdb)
            if af3_path and exp_path:
                yield ccd, pdb, af3_path, exp_path
            else:
                if not af3_path:
                    log.debug(f"No AF3 .cif for {ccd}/{pdb}")
                if not exp_path:
                    log.debug(f"No experimental .pdb for {ccd}/{pdb}")


# ===========================================================================
# 6.  Per-structure processing
# ===========================================================================

def process_pair(
    ccd: str,
    pdb: str,
    af3_path: str,
    exp_path: str,
    tool: str,
    dssr_bin: str,
    verbose: bool,
) -> dict:
    """
    Compute INF_WC and INF_all for one (AF3, experimental) pair.
    Returns a result dict.
    """
    result = {
        "CCD": ccd, "PDB": pdb,
        "af3_path": af3_path, "exp_path": exp_path,
        "INF_WC": "", "INF_all": "",
        "TP_WC": "", "FP_WC": "", "FN_WC": "",
        "TP_all": "", "FP_all": "", "FN_all": "",
        "note": "",
    }

    try:
        # --- annotate base pairs ---
        if tool == "dssr":
            exp_wc,  exp_all  = dssr_base_pairs(exp_path,  dssr_bin)
            af3_wc,  af3_all  = dssr_base_pairs(af3_path,  dssr_bin)
        else:
            exp_wc,  exp_all  = biopython_base_pairs(exp_path)
            af3_wc,  af3_all  = biopython_base_pairs(af3_path)
            
        if verbose:
            print(f"\n[DEBUG {ccd}/{pdb}] base-pair counts")
            print(f"  EXP WC : {len(exp_wc)}")
            print(f"  AF3 WC : {len(af3_wc)}")
            print(f"  EXP ALL: {len(exp_all)}")
            print(f"  AF3 ALL: {len(af3_all)}\n")

        # --- sequence alignment → residue mapping ---
        seq_exp = extract_sequence(exp_path)
        seq_af3 = extract_sequence(af3_path)
        
        if verbose:
            print("EXP length:", len(seq_exp))
            print("AF3 length:", len(seq_af3))

        if not seq_exp or not seq_af3:
            result["note"] = "empty sequence"
            return result

        al_exp, al_af3 = align_sequences(seq_exp, seq_af3)
        key_map = build_key_mapping(al_exp, al_af3)

        if len(key_map) < 3:
            result["note"] = "alignment too short"
            return result

        # Translate AF3 keys into experimental residue space
        af3_wc_t  = translate_pairs(af3_wc,  key_map)
        af3_all_t = translate_pairs(af3_all, key_map)

        # --- INF ---
        inf_wc  = inf_score(exp_wc,  af3_wc_t)
        inf_all = inf_score(exp_all, af3_all_t)

        def counts(ref, pred):
            tp = len(ref & pred)
            fp = len(pred - ref)
            fn = len(ref - pred)
            return tp, fp, fn

        tp_wc, fp_wc, fn_wc   = counts(exp_wc,  af3_wc_t)
        tp_all, fp_all, fn_all = counts(exp_all, af3_all_t)

        result.update({
            "INF_WC":  round(inf_wc,  4),
            "INF_all": round(inf_all, 4),
            "TP_WC": tp_wc, "FP_WC": fp_wc, "FN_WC": fn_wc,
            "TP_all": tp_all, "FP_all": fp_all, "FN_all": fn_all,
        })

        if verbose:
            log.info(
                f"{ccd}/{pdb}  INF_WC={inf_wc:.4f}  INF_all={inf_all:.4f}  "
                f"aligned={len(key_map)} residues"
            )

    except Exception as exc:
        result["note"] = str(exc)
        log.warning(f"{ccd}/{pdb}: {exc}")

    return result


# ===========================================================================
# 7.  Main
# ===========================================================================

CSV_FIELDS = [
    "CCD", "PDB",
    "INF_WC", "INF_all",
    "TP_WC", "FP_WC", "FN_WC",
    "TP_all", "FP_all", "FN_all",
    "note", "af3_path", "exp_path",
]


def main():
    parser = argparse.ArgumentParser(
        description="Compute INF_WC and INF_all for AF3 vs experimental structures."
    )
    parser.add_argument("--af3_root", required=True,
                        help="Root of MODNAP_AF3_models directory")
    parser.add_argument("--exp_root", required=True,
                        help="Root of MODNAP directory")
    parser.add_argument("--output", default="results_inf.csv",
                        help="Output CSV path (default: results_inf.csv)")
    parser.add_argument("--tool", choices=["dssr", "biopython", "auto"],
                        default="auto",
                        help="Base-pair annotation tool (default: auto)")
    parser.add_argument("--dssr_bin", default="x3dna-dssr",
                        help="Path to DSSR binary (default: x3dna-dssr)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print per-structure results to console")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only first N pairs (for testing)")
    args = parser.parse_args()

    # Resolve tool
    tool = args.tool
    if tool == "auto":
        tool = "dssr" if shutil.which(args.dssr_bin) else "biopython"
        log.info(f"Tool auto-detected: {tool}")
    else:
        if tool == "dssr" and not shutil.which(args.dssr_bin):
            log.warning(
                f"DSSR binary '{args.dssr_bin}' not found on PATH; "
                "falling back to biopython"
            )
            tool = "biopython"

    if not args.verbose:
        log.setLevel(logging.WARNING)

    log.info(f"AF3 root : {args.af3_root}")
    log.info(f"Exp root : {args.exp_root}")
    log.info(f"Output   : {args.output}")
    log.info(f"Tool     : {tool}")

    pairs = list(iter_pairs(args.af3_root, args.exp_root))
    if not pairs:
        sys.exit("No matched (AF3, experimental) pairs found. "
                 "Check --af3_root and --exp_root.")

    if args.limit:
        pairs = pairs[: args.limit]

    log.warning(f"Processing {len(pairs)} structure pairs…")

    rows = []
    for i, (ccd, pdb, af3_path, exp_path) in enumerate(pairs, 1):
        if args.verbose:
            log.warning(f"[{i}/{len(pairs)}] {ccd}/{pdb}")
        row = process_pair(ccd, pdb, af3_path, exp_path,
                           tool, args.dssr_bin, args.verbose)
        rows.append(row)

    # Write CSV
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    # Summary stats
    inf_wc_vals  = [r["INF_WC"]  for r in rows if isinstance(r["INF_WC"],  float)]
    inf_all_vals = [r["INF_all"] for r in rows if isinstance(r["INF_all"], float)]
    n_err = sum(1 for r in rows if r["note"])

    def _mean(vals):
        return round(sum(vals) / len(vals), 4) if vals else "n/a"

    print(f"\n{'='*55}")
    print(f"Structures processed : {len(rows)}")
    print(f"Errors / skipped     : {n_err}")
    print(f"Mean INF_WC          : {_mean(inf_wc_vals)}")
    print(f"Mean INF_all         : {_mean(inf_all_vals)}")
    print(f"Results written to   : {args.output}")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()
