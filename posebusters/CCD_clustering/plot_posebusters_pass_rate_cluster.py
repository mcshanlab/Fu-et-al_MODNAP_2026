#!/usr/bin/env python3
"""
plot_posebusters_pass_rate_cluster.py

Plots Chemical Validity pass/fail rates grouped by CCD_ID.

Each CCD gets ONE stacked horizontal bar:
    green = % passing PDBs
    red   = % failing PDBs

Example:
    CCD 5CM
        7 passing PDBs
        3 failing PDBs

    -> plotted as:
        70% green
        30% red

Usage
-----
python3 plot_posebusters_pass_rate_cluster.py \
    --csv ../combined_posebusters.csv \
    --outdir .

"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
    "svg.fonttype": "none",
    "pdf.use14corefonts": False,
})

import matplotlib.font_manager as fm

matplotlib.rcParams["font.family"] = (
    "Arial"
    if any("Arial" in f for f in fm.findSystemFonts())
    else "DejaVu Sans"
)

import matplotlib.pyplot as plt
import numpy as np


# ─────────────────────────────────────────────────────────────
# Chemical Validity parameters
# ─────────────────────────────────────────────────────────────

CHEM_VALIDITY_PARAMS = [
    "mol_pred_loaded",
    "sanitization",
    "inchi_convertible",
    "no_radicals",
]


# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

def _is_true(v):
    return str(v).strip().upper() == "TRUE"


def normalise_header(rows):
    """
    Map canonical parameter names to actual CSV headers.
    """

    actual = list(rows[0].keys())
    mapping = {}

    for p in CHEM_VALIDITY_PARAMS:
        for h in actual:
            if h.strip().lower() == p.lower():
                mapping[p] = h
                break

    return mapping


# ─────────────────────────────────────────────────────────────
# Core computation
# ─────────────────────────────────────────────────────────────

def compute_ccd_pass_fail(csv_path):
    """
    Returns:
        [(ccd, pass_count, fail_count), ...]

    grouped by CCD_ID column.
    """

    with open(csv_path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        raise ValueError(f"CSV is empty: {csv_path}")

    header_map = normalise_header(rows)

    counts = defaultdict(lambda: {"pass": 0, "fail": 0})

    for row in rows:

        # USE CCD_ID DIRECTLY
        ccd = str(row.get("CCD_ID", "")).strip().upper()

        if not ccd:
            continue

        passed = all(
            _is_true(row.get(header_map[p], ""))
            for p in CHEM_VALIDITY_PARAMS
        )

        if passed:
            counts[ccd]["pass"] += 1
        else:
            counts[ccd]["fail"] += 1

    # sort by total descending
    sorted_items = sorted(
        counts.items(),
        key=lambda x: x[1]["pass"] + x[1]["fail"],
        reverse=True,
    )

    return [
        (ccd, d["pass"], d["fail"])
        for ccd, d in sorted_items
    ]


# ─────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────

def plot_ccd_pass_fail(data, outdir):

    # reverse order so largest appears at top
    data_plot = list(reversed(data))

    ccds = [d[0] for d in data_plot]
    pass_counts = [d[1] for d in data_plot]
    fail_counts = [d[2] for d in data_plot]

    totals = [
        p + f
        for p, f in zip(pass_counts, fail_counts)
    ]

    pass_pcts = [
        (p / t) * 100 if t > 0 else 0
        for p, t in zip(pass_counts, totals)
    ]

    fail_pcts = [
        (f / t) * 100 if t > 0 else 0
        for f, t in zip(fail_counts, totals)
    ]

    n_ccds = len(ccds)

    # figure size
    row_h = 0.22
    fig_h = max(8, n_ccds * row_h)
    fig_w = 7

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    COLOR_PASS = "#7EC8A4"
    COLOR_FAIL = "#F4A7A3"

    y = np.arange(n_ccds)

    # PASS
    ax.barh(
        y,
        pass_pcts,
        color=COLOR_PASS,
        edgecolor="white",
        linewidth=0.3,
        height=0.8,
        label="Pass",
    )

    # FAIL
    ax.barh(
        y,
        fail_pcts,
        left=pass_pcts,
        color=COLOR_FAIL,
        edgecolor="white",
        linewidth=0.3,
        height=0.8,
        label="Fail",
    )

    # y-axis labels = CCD IDs
    ax.set_yticks(y)
    ax.set_yticklabels(ccds, fontsize=6)

    # annotations
    for i, (pct, total) in enumerate(zip(pass_pcts, totals)):

        ax.text(
            101,
            i,
            f"{pct:.0f}% (n={total})",
            va="center",
            ha="left",
            fontsize=5.5,
            color="#333333",
        )

    ax.set_xlim(0, 135)

    ax.set_xlabel(
        "Percentage of PDBs within each CCD (%)",
        fontsize=10,
    )

    ax.axvline(
        100,
        color="#888888",
        linewidth=0.6,
        linestyle="--",
    )

    total_pdbs = sum(totals)
    total_pass = sum(pass_counts)

    overall_pct = (
        total_pass / total_pdbs * 100
        if total_pdbs else 0
    )

    ax.set_title(
        f"Chemical Validity pass/fail grouped by CCD\n"
        f"({n_ccds} CCDs · "
        f"{total_pdbs} PDBs · "
        f"overall pass {overall_pct:.1f}%)",
        fontsize=10,
        pad=8,
    )

    ax.xaxis.grid(
        True,
        linestyle=":",
        linewidth=0.5,
        alpha=0.6,
    )

    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1,
                          color=COLOR_PASS,
                          label="Pass"),
            plt.Rectangle((0, 0), 1, 1,
                          color=COLOR_FAIL,
                          label="Fail"),
        ],
        loc="lower right",
        frameon=True,
        fontsize=8,
    )

    plt.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)

    png_path = outdir / "posebusters_ccd_pass_fail.png"
    pdf_path = outdir / "posebusters_ccd_pass_fail.pdf"

    fig.savefig(
        png_path,
        dpi=200,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"[INFO] Saved: {png_path}")
    print(f"[INFO] Saved: {pdf_path}")

    # TSV table
    tsv_path = outdir / "posebusters_ccd_pass_fail.tsv"

    with open(tsv_path, "w", encoding="utf-8") as fh:

        fh.write("CCD\tpass\tfail\ttotal\tpass_pct\n")

        for ccd, p, f in data:

            total = p + f

            pct = (p / total) * 100 if total else 0

            fh.write(
                f"{ccd}\t{p}\t{f}\t{total}\t{pct:.1f}\n"
            )

    print(f"[INFO] Saved table: {tsv_path}")


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Plot Chemical Validity pass/fail "
            "rates grouped by CCD_ID."
        )
    )

    parser.add_argument(
        "--csv",
        "-c",
        required=True,
        help="Input CSV file",
    )

    parser.add_argument(
        "--outdir",
        "-o",
        default=".",
        help="Output directory",
    )

    args = parser.parse_args()

    csv_path = Path(args.csv).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    if not csv_path.exists():
        print(f"[ERROR] CSV not found: {csv_path}")
        return

    print(f"[INFO] Reading: {csv_path}")

    data = compute_ccd_pass_fail(csv_path)

    print(f"[INFO] {len(data)} CCD IDs found")

    print("\nTop 10 CCDs:")

    for ccd, p, f in data[:10]:

        total = p + f

        pct = (p / total) * 100 if total else 0

        print(
            f"{ccd:8s} "
            f"pass={p:4d} "
            f"fail={f:4d} "
            f"total={total:4d} "
            f"pass%={pct:5.1f}"
        )

    plot_ccd_pass_fail(data, outdir)


if __name__ == "__main__":
    main()
