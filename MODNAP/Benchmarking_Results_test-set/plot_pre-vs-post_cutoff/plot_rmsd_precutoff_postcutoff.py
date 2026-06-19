#!/usr/bin/env python3
"""
useage: python3 plot_rmsd_precutoff_postcutoff.py

Reads:
    pymol_rmsd_summary_mod.txt           — RMSD + Mod RMSD values per CCD/PDB
    pre_cutoff.txt                       — PDBs released on or before 2021-09-30
    post_cutoff.txt                      — PDBs released after 2021-09-30
    nuc_type_count_per_pdb.txt           — category classification per CCD/PDB

Produces 4 figures (PNG + PDF):
    1. full_rmsd_prepost.png   — whole-structure RMSD, % below 2.0 Å
    2. mod_rmsd_prepost.png    — modification RMSD, % below 2.0 Å

Each figure has 6 groups (DNA, RNA, DNA-protein, RNA-protein, DNA-RNA, DNA-RNA-protein).
Each group has 2 bars: pre-cutoff and post-cutoff.
Y-axis = % of PDBs with RMSD < 2.0 Å.
Bar labels show: XX.X% (n/total)
"""

from pathlib import Path
from collections import defaultdict, Counter
import re

import matplotlib
import matplotlib.font_manager as fm

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['pdf.use14corefonts'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'

arial_exists = any("Arial" in f for f in fm.findSystemFonts())
matplotlib.rcParams['font.family'] = 'Arial' if arial_exists else 'DejaVu Sans'

import numpy as np
import matplotlib.pyplot as plt

# --- Paths ---
RMSD_FILE   = Path("../rmsd_modification/pymol_rmsd_summary_mod.txt")
PRE_FILE    = Path("../../MODNAP_features_pre-vs-post-cutoff/pre_cutoff.txt")
POST_FILE   = Path("../../MODNAP_features_pre-vs-post-cutoff/post_cutoff.txt")
NUC_FILE    = Path("../../../MODNAP/MODNAP_features/nuc_type_count_per_pdb.txt")
OUT_DIR     = Path(".")

CATEGORIES  = ["DNA", "RNA", "DNA-protein", "RNA-protein", "DNA-RNA", "DNA-RNA-protein"]
RMSD_CUTOFF = 2.0  # Å

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_rmsd_file(path: Path) -> list[dict]:
    """
    Parse pymol_rmsd_summary_mod.txt.
    Expected header: CCD  PDB  RMSD  N_atoms  Mod_RMSD  Mod_N_atoms
    Returns list of dicts with keys: ccd, pdb, rmsd, mod_rmsd.
    Skips header, comment lines, and rows where values are not numeric.
    """
    records = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            if parts[0].strip().upper() == "CCD":
                continue
            try:
                records.append({
                    "ccd":      parts[0].strip().upper(),
                    "pdb":      parts[1].strip().upper(),
                    "rmsd":     float(parts[2].strip()),
                    "mod_rmsd": float(parts[4].strip()),
                })
            except ValueError:
                continue
    return records


def parse_cutoff_file(path: Path) -> set[str]:
    """
    Parse pre_cutoff.txt or post_cutoff.txt.
    Returns a set of PDB IDs (upper-cased).
    Skips comment lines (#) and the header row.
    """
    pdbs = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            if parts[0].strip().upper() == "CCD":
                continue
            pdb = parts[1].strip().upper()
            pdbs.add(pdb)
    return pdbs


def parse_nuc_type_file(path: Path) -> dict:
    """
    Parse nuc_type_count_per_pdb.txt into a mapping (ccd, pdb) -> category_label.
    Taken directly from plot_openstructure_cdfs.py.
    """
    mapping = {}
    if not path.exists():
        print(f"[WARN] nuc_type_count_per_pdb.txt not found: {path}")
        return mapping

    text = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    if not text:
        return mapping

    header = text[0]
    delim = '\t' if '\t' in header else None
    cols = header.split(delim) if delim else header.split()

    try:
        idx_ccd      = cols.index("CCD")
        idx_pdb      = cols.index("PDB")
        idx_category = cols.index("Category")
    except ValueError:
        lower = [c.lower() for c in cols]
        try:
            idx_ccd      = lower.index("ccd")
            idx_pdb      = lower.index("pdb")
            idx_category = lower.index("category")
        except ValueError:
            print(f"[WARN] Could not find CCD/PDB/Category columns in: {path}")
            return mapping

    for line in text[1:]:
        if not line.strip():
            continue
        parts = line.split(delim) if delim else line.split()
        if len(parts) <= max(idx_ccd, idx_pdb, idx_category):
            continue

        ccd = parts[idx_ccd].strip()
        if ccd.lower().endswith("_ok"):
            ccd = ccd[:-3]
        ccd = ccd.upper()

        pdb = parts[idx_pdb].strip().upper()

        category = parts[idx_category].strip()
        # normalize DNA+RNA variants to match CATEGORIES list
        if category == "DNA+RNA":
            category = "DNA-RNA"
        elif category == "DNA+RNA-protein":
            category = "DNA-RNA-protein"

        mapping[(ccd, pdb)] = category

    return mapping


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_category_map(records: list[dict], nuc_map: dict) -> dict:
    """
    Map (ccd, pdb) -> category by looking up nuc_map.
    Tries several key variants to tolerate '_ok' suffixes.
    """
    category_map = {}
    for rec in records:
        ccd, pdb = rec["ccd"], rec["pdb"]
        cat = None
        for key in [(ccd, pdb), (ccd.rstrip("_ok"), pdb), (ccd, pdb.rstrip("_ok"))]:
            if key in nuc_map:
                cat = nuc_map[key]
                break
        if cat is None:
            print(f"  [SKIP] No category found for {ccd}/{pdb}")
            continue
        category_map[(ccd, pdb)] = cat
    return category_map


def compute_percentages(records, category_map, pre_pdbs, post_pdbs, rmsd_key):
    """
    Count total PDBs per category (pre/post) and how many are below RMSD_CUTOFF.
    Returns percentage dicts and raw count dicts.
    """
    pre_total = defaultdict(int)
    post_total = defaultdict(int)
    pre_below  = defaultdict(int)
    post_below = defaultdict(int)

    for rec in records:
        key = (rec["ccd"], rec["pdb"])
        cat = category_map.get(key)
        if cat is None:
            continue

        pdb_upper = rec["pdb"].upper()
        val = rec[rmsd_key]

        if pdb_upper in pre_pdbs:
            pre_total[cat] += 1
            if val <= RMSD_CUTOFF:
                pre_below[cat] += 1
        elif pdb_upper in post_pdbs:
            post_total[cat] += 1
            if val <= RMSD_CUTOFF:
                post_below[cat] += 1

    def pct(below, total):
        return 100.0 * below / total if total > 0 else 0.0

    pre_pcts  = {cat: pct(pre_below[cat],  pre_total[cat])  for cat in CATEGORIES}
    post_pcts = {cat: pct(post_below[cat], post_total[cat]) for cat in CATEGORIES}

    return pre_pcts, post_pcts, pre_total, post_total, pre_below, post_below


def plot_figure(pre_pcts, post_pcts, pre_total, post_total, pre_below, post_below, title, out_stem):
    """
    Grouped bar chart: 6 categories x 2 bars (pre/post cutoff).
    Y-axis = % of PDBs with RMSD <= RMSD_CUTOFF.
    Bar labels: XX.X% \n below/total
    """
    bar_width = 0.35
    x = np.arange(len(CATEGORIES))

    pre_color  = "#4C9BE8"
    post_color = "#F4A460"

    fig, ax = plt.subplots(figsize=(13, 6))

    bars_pre  = ax.bar(x - bar_width / 2,
                       [pre_pcts[cat]  for cat in CATEGORIES],
                       bar_width, label="Pre-cutoff (≤ 2021-09-30)",
                       color=pre_color, alpha=0.85, edgecolor="white")

    bars_post = ax.bar(x + bar_width / 2,
                       [post_pcts[cat] for cat in CATEGORIES],
                       bar_width, label="Post-cutoff (> 2021-09-30)",
                       color=post_color, alpha=0.85, edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(CATEGORIES, fontsize=10)
    ax.set_ylabel(f"% of PDBs with RMSD ≤ {RMSD_CUTOFF} Å", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_ylim(0, 115)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    for bars, totals, belows in [
        (bars_pre,  pre_total,  pre_below),
        (bars_post, post_total, post_below),
    ]:
        for bar, cat in zip(bars, CATEGORIES):
            h = bar.get_height()
            n = totals[cat]
            b = belows[cat]
            label = f"{h:.1f}%\n{b}/{n}"
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                h + 1.5,
                label,
                ha="center", va="bottom", fontsize=8
            )

    ax.legend(fontsize=9)
    plt.tight_layout()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUT_DIR / f"{out_stem}.png"
    pdf_path = OUT_DIR / f"{out_stem}.pdf"
    plt.savefig(str(png_path), dpi=200)
    plt.savefig(str(pdf_path))
    plt.close()
    print(f"  Saved: {png_path}")
    print(f"  Saved: {pdf_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    for path in [RMSD_FILE, PRE_FILE, POST_FILE, NUC_FILE]:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    print("Parsing RMSD file...", flush=True)
    records = parse_rmsd_file(RMSD_FILE)
    print(f"  Total records: {len(records)}", flush=True)

    print("Parsing pre/post cutoff files...", flush=True)
    pre_pdbs  = parse_cutoff_file(PRE_FILE)
    post_pdbs = parse_cutoff_file(POST_FILE)
    print(f"  Pre-cutoff PDBs:  {len(pre_pdbs)}", flush=True)
    print(f"  Post-cutoff PDBs: {len(post_pdbs)}", flush=True)

    print("\nLoading category map from nuc_type_count_per_pdb.txt...", flush=True)
    nuc_map = parse_nuc_type_file(NUC_FILE)
    print(f"  Entries in nuc_type map: {len(nuc_map)}", flush=True)

    category_map = build_category_map(records, nuc_map)
    print(f"  Categorized: {len(category_map)} / {len(records)} records", flush=True)

    cat_summary = Counter(category_map.values())
    for cat in CATEGORIES:
        print(f"    {cat}: {cat_summary.get(cat, 0)}", flush=True)

    # --- Figure 1: Full structure RMSD ---
    print("\nComputing full structure RMSD percentages...", flush=True)
    pre_pcts, post_pcts, pre_total, post_total, pre_below, post_below = compute_percentages(
        records, category_map, pre_pdbs, post_pdbs, "rmsd"
    )
    plot_figure(
        pre_pcts, post_pcts, pre_total, post_total, pre_below, post_below,
        title=f"Full Structure RMSD ≤ {RMSD_CUTOFF} Å by Category and Release Period",
        out_stem="full_rmsd_prepost"
    )

    # --- Figure 2: Modification RMSD ---
    print("\nComputing modification RMSD percentages...", flush=True)
    pre_pcts, post_pcts, pre_total, post_total, pre_below, post_below = compute_percentages(
        records, category_map, pre_pdbs, post_pdbs, "mod_rmsd"
    )
    plot_figure(
        pre_pcts, post_pcts, pre_total, post_total, pre_below, post_below,
        title=f"Modification RMSD ≤ {RMSD_CUTOFF} Å by Category and Release Period",
        out_stem="mod_rmsd_prepost"
    )

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
