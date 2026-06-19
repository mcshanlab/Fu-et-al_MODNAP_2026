#!/usr/bin/env python3
"""
usage: python3 plot_posebusters_pass_rate.py

Reads a PoseBusters results CSV and produces:
  1. A column (bar) figure showing the passing percentage for each
     quality-check parameter.
  2. A scatter/dot figure showing Chemical Validity and Intramolecular
     Validity pass rates broken down by nucleotide/protein category
     (RNA, DNA, RNA-protein, DNA-protein, DNA+RNA, DNA+RNA-protein),
     derived by walking the AF3 JSON files under openstructure_last/.

Rules:
  - TRUE  -> pass
  - FALSE -> fail
  - empty -> fail (PoseBusters could not evaluate; treated as failure)

  For the category figure:
  - Chemical Validity   : a PDB passes iff ALL 8 chem-validity params are TRUE
  - Intramolecular Val. : a PDB passes iff ALL 9 intramol-validity params are TRUE
  - JSON folder         : /Volumes/My_Passport/Clary/step19_posebusters/openstructure_last
                          structure: <CCD_dir>/<PDB_dir>/<single .json>
  - Linking             : the 'file_name' (or 'mol_pred'/'id') column in the CSV
                          is matched case-insensitively against PDB folder names.

Output:
  - posebusters_pass_rate.png / .pdf
  - posebusters_by_nuc_category.png / .pdf
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

matplotlib.rcParams['pdf.fonttype']        = 42
matplotlib.rcParams['ps.fonttype']         = 42
matplotlib.rcParams['text.usetex']         = False
matplotlib.rcParams['svg.fonttype']        = 'none'
matplotlib.rcParams['pdf.use14corefonts'] = False

import matplotlib.font_manager as fm
_arial = any("Arial" in f for f in fm.findSystemFonts())
matplotlib.rcParams['font.family'] = 'Arial' if _arial else 'DejaVu Sans'

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ── Parameter definitions ─────────────────────────────────────────────────────

PARAMS = [
    "mol_pred_loaded",
    "sanitization",
    "inchi_convertible",
    "no_radicals",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "internal_energy",
]

PARAM_LABELS = [
    "mol_pred\nloaded",
    "sanitization",
    "inchi\nconvertible",
    "no_radicals",
    "bond\nlengths",
    "bond\nangles",
    "internal\nsteric_clash",
    "aromatic_ring\nflatness",
    "non-aromatic\nring non-flat",
    "double_bond\nflatness",
    "internal\nenergy",
]

CHEM_VALIDITY_PARAMS = [
    "mol_pred_loaded",
    "sanitization",
    "inchi_convertible",
    "no_radicals",
]

INTRAMOL_VALIDITY_PARAMS = [
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "internal_energy",
]

NUC_CATEGORIES = [
    "RNA",
    "DNA",
    "RNA-protein",
    "DNA-protein",
    "DNA+RNA",
    "DNA+RNA-protein",
]

NUC_CATEGORY_LABELS = [
    "RNA",
    "DNA",
    "RNA-\nprotein",
    "DNA-\nprotein",
    "DNA+RNA\nhybrid",
    "DNA+RNA\n+protein",
]

PASTEL_COLORS = [
    "#FFB3B3", "#FFDAB9", "#FFE599", "#C9E4CA", "#A8D8EA", "#D0C4F7",
    "#F7C6D4", "#B5EAD7", "#FFDFD3", "#C7CEEA", "#E2F0CB", "#FFD6E0",
    "#BDE0FE", "#FFCBA4", "#CDE5D9", "#F6C7B3", "#D0F0FF", "#E8D5B7",
]

# ── MODNAP classification helpers ───────────────────────────────────────────

def build_pdb_to_category_map(nuc_type_txt: Path) -> dict:
    """
    Reads nuc_type_count_per_pdb.txt

    Expected columns:
      CCD    PDB    Nucleotide    Has_Protein    Category    JSON

    Returns:
      {
          "1r3o": "RNA",
          "1mq2": "DNA-protein",
          ...
      }
    """

    mapping = {}

    with nuc_type_txt.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            pdb = str(row.get("PDB") or "").strip().lower()
            category = str(row.get("Category") or "").strip()

            if not pdb or not category:
                continue

            mapping[pdb] = category

    return mapping

# ── CSV helpers ─────────────────────────────────────────────────────────────
def normalize_pdb_key(x: str) -> str:
    x = str(x).strip().lower()

    # remove file path if present
    x = Path(x).name

    # remove extension if present
    x = Path(x).stem

    # strip known suffixes (IMPORTANT ORDER FIX)
    for suffix in ["_ligand", "_pred", "_true", "_mol", "_ok"]:
        if x.endswith(suffix):
            x = x[: -len(suffix)]

    return x
    
def _normalise_headers(rows):
    if not rows:
        return {}
    actual = list(rows[0].keys())
    mapping = {}
    for p in PARAMS:
        for h in actual:
            if h.strip().lower() == p.lower():
                mapping[p] = h
                break
    return mapping


def _is_true(v: str) -> bool:
    return str(v).strip().upper() == "TRUE"

def compute_category_pass_rates(rows, header_map, id_col, pdb_to_category):

    buckets = {c: [] for c in NUC_CATEGORIES}
    buckets["unknown"] = []

    for row in rows:
        raw = row.get(id_col, "")

        key = normalize_pdb_key(raw)

        category = pdb_to_category.get(key)

        if category is None:
            print(f"[DEBUG unmatched] raw='{raw}' → key='{key}'")
            buckets["unknown"].append(row)
            continue

        if category in buckets:
            buckets[category].append(row)
        else:
            buckets["unknown"].append(row)

    results = {}

    for cat in NUC_CATEGORIES:
        group = buckets[cat]
        n = len(group)

        if n == 0:
            results[cat] = {
                "n": 0,
                "chem_pct": None,
                "intramol_pct": None
            }
            continue

        n_chem = sum(
            all(_is_true(r.get(header_map.get(p, ""), ""))
                for p in CHEM_VALIDITY_PARAMS)
            for r in group
        )

        n_intr = sum(
            all(_is_true(r.get(header_map.get(p, ""), ""))
                for p in INTRAMOL_VALIDITY_PARAMS)
            for r in group
        )

        results[cat] = {
            "n": n,
            "chem_pct": n_chem / n * 100,
            "intramol_pct": n_intr / n * 100,
        }

    return results, buckets["unknown"]


# ── FIGURE 1 unchanged ──────────────────────────────────────────────────────
def plot_pass_rates(pass_rates, n_total, outdir):
    values = [pass_rates[p] for p in PARAMS]
    colors = [PASTEL_COLORS[i % len(PASTEL_COLORS)] for i in range(len(PARAMS))]

    fig, ax = plt.subplots(figsize=(14, 6))
    x = range(len(PARAMS))

    ax.bar(x, values, color=colors, edgecolor='black', linewidth=0.6)

    ax.set_xticks(list(x))
    ax.set_xticklabels(PARAM_LABELS, fontsize=7.5)
    ax.set_ylabel("Passing percentage (%)")
    ax.set_ylim(0, 115)
    ax.axhline(100, color='#888888', linewidth=0.8, linestyle='--')

    for i, v in enumerate(values):
        ax.text(i, v + 1.5, f"{v:.1f}%", ha='center', fontsize=7)

    plt.tight_layout()
    outdir.mkdir(exist_ok=True)

    fig.savefig(outdir / "posebusters_pass_rate.png", dpi=200)
    fig.savefig(outdir / "posebusters_pass_rate.pdf", dpi=300)
    plt.close(fig)

# ── FIGURE 2 unchanged ──────────────────────────────────────────────────────
def plot_category_pass_rates(cat_results: dict, n_total: int, outdir: Path):
    COLOR_CHEM = "#C9A8D4"
    COLOR_INTR = "#B8CCE4"

    fig, ax = plt.subplots(figsize=(9, 5))

    x_positions = list(range(len(NUC_CATEGORIES)))

    chem_plotted = False
    intramol_plotted = False

    for xi, cat in enumerate(NUC_CATEGORIES):
        res = cat_results.get(cat, {})
        n = res.get("n", 0)

        chem = res.get("chem_pct", None)
        intr = res.get("intramol_pct", None)

        if n == 0 or chem is None or intr is None:
            ax.scatter(xi, 0, marker='o', s=60, color='#cccccc', alpha=0.5, zorder=3)
            ax.scatter(xi, 0, marker='s', s=60, color='#cccccc', alpha=0.5, zorder=3)
            continue

        ax.scatter(
            xi, chem,
            marker='o', s=100,
            color=COLOR_CHEM,
            edgecolors='#9b7db5', linewidths=0.8,
            zorder=4, alpha=0.92,
            label="Chemical Validity" if not chem_plotted else ""
        )
        ax.scatter(
            xi, intr,
            marker='s', s=90,
            color=COLOR_INTR,
            edgecolors='#7a9bbf', linewidths=0.8,
            zorder=4, alpha=0.92,
            label="Intramolecular Validity" if not intramol_plotted else ""
        )

        ax.text(xi, chem + 2, f"{chem:.1f}%", ha='center', fontsize=7, color='#333333')
        ax.text(xi, intr - 4, f"{intr:.1f}%", ha='center', fontsize=7, color='#333333')

        chem_plotted = True
        intramol_plotted = True

    # ---- X labels (MATCH STYLE) ----
    x_labels = []
    for cat, lab in zip(NUC_CATEGORIES, NUC_CATEGORY_LABELS):
        n = cat_results.get(cat, {}).get("n", 0)
        x_labels.append(f"{lab}\n(n={n})")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=7.5, ha='center')

    # ---- Axes styling (MATCH FIGURE 1) ----
    ax.set_ylabel("Passing percentage (%)", fontsize=10)
    ax.set_ylim(0, 115)
    ax.set_xlim(-0.6, len(NUC_CATEGORIES) - 0.4)

    ax.set_title(
        f"PoseBusters pass rates by nucleotide/protein category  (n = {n_total} PDB structures)",
        fontsize=11, pad=10,
    )

    ax.axhline(100, color='#888888', linewidth=0.8, linestyle='--')

    ax.yaxis.grid(True, linestyle=':', linewidth=0.5, alpha=0.7, color='#aaaaaa')
    ax.set_axisbelow(True)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ---- Legend (MATCH STYLE) ----
    legend_handles = [
        mlines.Line2D([], [], marker='o', color='w',
                      markerfacecolor=COLOR_CHEM,
                      markeredgecolor='#9b7db5',
                      markersize=9,
                      label="Chemical Validity"),
        mlines.Line2D([], [], marker='s', color='w',
                      markerfacecolor=COLOR_INTR,
                      markeredgecolor='#7a9bbf',
                      markersize=9,
                      label="Intramolecular Validity"),
    ]

    ax.legend(
        handles=legend_handles,
        loc='upper right',
        frameon=False,
        fontsize=8.5,
    )

    plt.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)

    fig.savefig(outdir / "posebusters_by_nuc_category.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / "posebusters_by_nuc_category.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_pass_rates(csv_path: Path):
    with csv_path.open() as fh:
        rows = list(csv.DictReader(fh))

    header_map = _normalise_headers(rows)
    id_col = "PDB_ID"

    n_total = len(rows)
    pass_rates = {}

    print("[DEBUG CSV columns]")
    print(list(rows[0].keys()))
    
    for p in PARAMS:
        col = header_map.get(p)
        if not col:
            pass_rates[p] = 0
            continue

        pass_rates[p] = (
            sum(_is_true(r.get(col, "")) for r in rows)
            / n_total * 100
        )

    return pass_rates, n_total, rows, header_map, id_col

# ── MAIN unchanged ──────────────────────────────────────────────────────────

def main():
    csv_path = Path("combined_posebusters.csv")
    nuc_type_txt = Path("../MODNAP/MODNAP_features/nuc_type_count_per_pdb.txt")

    pass_rates, n_total, rows, header_map, id_col = compute_pass_rates(csv_path)

    print(f"[INFO] Reading CSV: {csv_path}")
    print(f"[INFO] Total rows: {n_total}")
    print(f"[INFO] ID column: {id_col}")

    for p in PARAMS:
        print(f"  {p:<35s} {pass_rates[p]:.2f}%")

    outdir = Path(".")

    plot_pass_rates(pass_rates, n_total, outdir)

    pdb_to_category = build_pdb_to_category_map(nuc_type_txt)
    print(f"[INFO] PDB classifications loaded: {len(pdb_to_category)}")
    print("[DEBUG] example PDB keys:", list(pdb_to_category.keys())[:10])
    
    # 👇 ADD THIS BLOCK HERE
    print("\n[DEBUG] EXAMPLE MATCH vs MISMATCH CHECK")

    matched_example = None
    unmatched_example = None

    for row in rows:
        raw = row.get(id_col, "")
        key = normalize_pdb_key(raw)

        if key in pdb_to_category and matched_example is None:
            matched_example = (raw, key, pdb_to_category[key])

        if key not in pdb_to_category and unmatched_example is None:
            unmatched_example = (raw, key)

        if matched_example and unmatched_example:
            break

    print("\n[EXAMPLE MATCHED]")
    if matched_example:
        raw, key, path = matched_example
        print(f"RAW: {raw}")
        print(f"KEY: {key}")
        print(f"CATEGORY: {path}")
    else:
        print("No matched example found")

    print("\n[EXAMPLE UNMATCHED]")
    if unmatched_example:
        raw, key = unmatched_example
        print(f"RAW: {raw}")
        print(f"KEY: {key}")
        print(f"IN MAP?: {key in pdb_to_category}")
    else:
        print("No unmatched example found")

# 👆 END BLOCK

    cat_results, unknown_rows = compute_category_pass_rates(
        rows, header_map, id_col, pdb_to_category
    )
    
    print("[INFO] Category results:")
    
    for c, r in cat_results.items():
        if "n" in r:
            print(
                f"  {c:<20s} n={r['n']:4d} "
                f"chem={r.get('chem_pct')} "
                f"intramol={r.get('intramol_pct')}"
            )
            
    plot_category_pass_rates(cat_results, n_total, outdir)

    print("\n[INFO] UNMATCHED / UNKNOWN ENTRIES (these did not map to nuc_type_count_per_pdb.txt):")
    print(f"Total unmatched: {len(unknown_rows)}\n")

    for r in unknown_rows:
        raw = r.get(id_col, "")
        key = normalize_pdb_key(raw)
        print(f"RAW: {raw} | KEY: {key}")
    
if __name__ == "__main__":
    main()
