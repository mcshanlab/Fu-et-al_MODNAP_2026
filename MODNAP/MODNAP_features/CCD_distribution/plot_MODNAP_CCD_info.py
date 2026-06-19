#!/usr/bin/env python3
"""
Plot a compact horizontal bar chart of PDB count per CCD,
sorted descending, designed to fit on a standard page.

Usage:
python3 plot_MODNAP_CCD_info.py --base ../../MODNAP --outdir .
"""

from pathlib import Path
import argparse
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']        = False
matplotlib.rcParams['svg.fonttype']       = 'none'
matplotlib.rcParams['pdf.use14corefonts'] = False

import matplotlib.font_manager as fm
arial_exists = any("Arial" in f for f in fm.findSystemFonts())
matplotlib.rcParams['font.family'] = 'Arial' if arial_exists else 'DejaVu Sans'

import matplotlib.pyplot as plt
import numpy as np


def save_fig(outpath: Path, dpi: int = 200):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outpath), dpi=dpi, bbox_inches="tight")
    try:
        plt.savefig(str(outpath.with_suffix(".pdf")),
                    format="pdf", bbox_inches="tight")
    except Exception as e:
        print(f"[WARN] PDF save failed: {e}")
    plt.close()
    print(f"[INFO] Saved: {outpath}  (+ .pdf)")


def count_pdbs_per_ccd(base: Path) -> dict:
    counts = {}
    for ccd_dir in sorted(base.iterdir()):
        if not ccd_dir.is_dir():
            continue
        n = sum(1 for p in ccd_dir.iterdir() if p.is_dir())
        if n > 0:
            counts[ccd_dir.name] = n
    return counts


def plot_horizontal_bar(counts: dict, outpath: Path):
    if not counts:
        print("[WARN] No CCD counts found.")
        return

    # sort ascending so the longest bar is at top when rendered
    sorted_items = sorted(counts.items(), key=lambda x: x[1])
    labels = [item[0] for item in sorted_items]
    values = [item[1] for item in sorted_items]

    n_ccds   = len(labels)
    max_v    = max(values)

    # --- page-friendly dimensions ---
    # 0.18 in per row fits ~129 CCDs on a standard letter/A4 height
    row_h    = 0.18
    fig_h    = max(8, n_ccds * row_h)
    fig_w    = 7          # fits within a standard page width

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap   = plt.cm.Blues
    colors = [cmap(0.30 + 0.70 * (v / max_v)) for v in values]

    y = np.arange(n_ccds)
    ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.3, height=0.8)

    # CCD labels on the y-axis
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)

    # count label to the right of each bar
    x_off = max_v * 0.01
    for i, v in enumerate(values):
        ax.text(v + x_off, i, str(v),
                va='center', ha='left', fontsize=5.5)

    # add a little right-margin so the count labels aren't clipped
    ax.set_xlim(right=max_v * 1.15)

    ax.set_xlabel("Number of PDB structures", fontsize=10)
    ax.set_title(
        f"PDB count per CCD  (n = {n_ccds} CCDs, "
        f"{sum(values)} total PDBs)",
        fontsize=11, pad=8,
    )

    ax.xaxis.grid(True, linestyle=':', linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_fig(outpath)

    # companion TSV
    tsv_path = outpath.with_suffix(".tsv")
    with tsv_path.open("w", encoding="utf-8") as fh:
        fh.write("CCD\tPDB_count\n")
        for ccd, n in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            fh.write(f"{ccd}\t{n}\n")
    print(f"[INFO] Table saved: {tsv_path}")


def main():
    default_base = (Path.home() / "Desktop" / "ML_aptamer_project"
                    / "step12_benchmark" / "MODNAP")
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", "-b", default=str(default_base),
                        help="MODNAP base directory")
    parser.add_argument("--outdir", "-o", default=None,
                        help="Output directory (default: current dir)")
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    if not base.exists() or not base.is_dir():
        print(f"[ERROR] Directory not found: {base}")
        return

    outdir = (Path(args.outdir).expanduser().resolve()
              if args.outdir else Path(".").resolve())

    print(f"[INFO] Counting PDBs per CCD in: {base}")
    counts = count_pdbs_per_ccd(base)
    print(f"[INFO] {len(counts)} CCDs, {sum(counts.values())} total PDBs")

    plot_horizontal_bar(counts, outdir / "pdb_count_per_ccd.png")
    print("[DONE]")


if __name__ == "__main__":
    main()
