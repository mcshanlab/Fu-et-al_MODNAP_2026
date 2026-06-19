#!/usr/bin/env python3
"""
usage: python3 plot_rmsd_by_category.py

Reads pymol_rmsd_summary_mod.txt and categorizes each PDB
using nuc_type_count_per_pdb.txt into:
DNA, RNA, DNA-protein, RNA-protein, DNA-RNA, DNA-RNA-protein.

Produces 2 figures (PNG + PDF):
  1. full_rmsd_by_category.png  — whole-structure RMSD per category
  2. mod_rmsd_by_category.png   — modification RMSD per category

Each figure shows grouped bars for the 6 categories.
Each bar = % of PDBs in that category with RMSD below the cutoff.
Cutoffs: 2.0, 2.5, 3.0 Å
"""

from pathlib import Path
from collections import Counter

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
NUC_FILE    = Path("../../../MODNAP/MODNAP_features/nuc_type_count_per_pdb.txt")
OUT_DIR   = Path(".")

CATEGORIES = ["DNA", "RNA", "DNA-protein", "RNA-protein", "DNA-RNA", "DNA-RNA-protein"]
CUTOFFS    = [2.0, 2.5, 3.0]


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_rmsd_file(path: Path):
    """
    Parse pymol_rmsd_summary_mod.txt.
    Expected header: CCD  PDB  RMSD  N_atoms  Mod_RMSD  Mod_N_atoms
    Returns list of dicts with keys: ccd, pdb, rmsd, n_atoms, mod_rmsd, mod_n_atoms.
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
                    "ccd":         parts[0].strip().upper(),
                    "pdb":         parts[1].strip().upper(),
                    "rmsd":        float(parts[2].strip()),
                    "n_atoms":     parts[3].strip(),
                    "mod_rmsd":    float(parts[4].strip()),
                    "mod_n_atoms": parts[5].strip() if len(parts) > 5 else "",
                })
            except ValueError:
                continue
    return records


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
        for key in [
            (ccd, pdb),
            (ccd[:-3] if ccd.endswith("_OK") else ccd, pdb),
            (ccd, pdb[:-3] if pdb.endswith("_OK") else pdb),
        ]:
            if key in nuc_map:
                cat = nuc_map[key]
                break
        if cat is None:
            print(f"  [SKIP] No category found for {ccd}/{pdb}")
            continue
        category_map[(ccd, pdb)] = cat
    return category_map


def compute_percentages(records, category_map, rmsd_key, cutoffs):
    """
    For each category and cutoff, compute % of PDBs with RMSD <= cutoff.
    Returns dict: category -> list of percentages (one per cutoff),
    and total counts per category.
    """
    cat_vals = {cat: [] for cat in CATEGORIES}
    for rec in records:
        key = (rec["ccd"], rec["pdb"])
        cat = category_map.get(key)
        if cat is None:
            continue
        cat_vals[cat].append(rec[rmsd_key])

    cat_pcts   = {}
    cat_counts = {}
    for cat in CATEGORIES:
        vals  = cat_vals[cat]
        total = len(vals)
        cat_counts[cat] = total
        if total == 0:
            cat_pcts[cat] = [0.0] * len(cutoffs)
        else:
            cat_pcts[cat] = [100.0 * sum(v <= c for v in vals) / total for c in cutoffs]

    return cat_pcts, cat_counts


def plot_figure(cat_pcts, cat_counts, cutoffs, title, out_stem):
    """
    Grouped bar chart: x = categories, grouped bars = cutoffs.
    """
    n_cutoffs = len(cutoffs)
    bar_width = 0.22
    x = np.arange(len(CATEGORIES))

    colors = ["#4C9BE8", "#F4A460", "#6DBF7E"]

    fig, ax = plt.subplots(figsize=(12, 5.5))

    bars_all = []
    for i, (cutoff, color) in enumerate(zip(cutoffs, colors)):
        pcts   = [cat_pcts[cat][i] for cat in CATEGORIES]
        offset = (i - (n_cutoffs - 1) / 2.0) * bar_width
        bars   = ax.bar(x + offset, pcts, bar_width,
                        label=f"RMSD ≤ {cutoff} Å",
                        color=color, alpha=0.85, edgecolor="white", linewidth=0.5)
        bars_all.append(bars)

    x_labels = [f"{cat}\n(n={cat_counts[cat]})" for cat in CATEGORIES]
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_ylabel("% of PDBs below RMSD cutoff", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_ylim(0, 110)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    for bars in bars_all:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2.0, h + 1.5,
                        f"{h:.1f}%", ha="center", va="bottom", fontsize=7.5)

    ax.legend(title="Cutoff", fontsize=9, title_fontsize=9)
    plt.tight_layout()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUT_DIR / f"{out_stem}.png"
    pdf_path = OUT_DIR / f"{out_stem}.pdf"
    plt.savefig(str(png_path), dpi=200)
    plt.savefig(str(pdf_path))
    plt.close()
    print(f"  Saved: {png_path}")
    print(f"  Saved: {pdf_path}")
    return png_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    for path in [RMSD_FILE, NUC_FILE]:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    print("Reading RMSD summary file...", flush=True)
    records = parse_rmsd_file(RMSD_FILE)
    print(f"  Total records: {len(records)}", flush=True)

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
    cat_pcts_full, cat_counts_full = compute_percentages(
        records, category_map, "rmsd", CUTOFFS)
    plot_figure(
        cat_pcts_full, cat_counts_full, CUTOFFS,
        title="Full Structure RMSD by Category",
        out_stem="full_rmsd_by_category"
    )

    # --- Figure 2: Modification RMSD ---
    print("\nComputing modification RMSD percentages...", flush=True)
    cat_pcts_mod, cat_counts_mod = compute_percentages(
        records, category_map, "mod_rmsd", CUTOFFS)
    plot_figure(
        cat_pcts_mod, cat_counts_mod, CUTOFFS,
        title="Modification RMSD by Category",
        out_stem="mod_rmsd_by_category"
    )

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
