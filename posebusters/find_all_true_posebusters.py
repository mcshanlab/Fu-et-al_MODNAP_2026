#!/usr/bin/env python3
"""
Filter PoseBusters results: find PDB IDs (grouped by CCD ID) that pass ALL
parameters, ignoring the rmsd column.

Parameters checked:
    mol_pred_loaded, mol_true_loaded, sanitization, inchi_convertible,
    all_atoms_connected, no_radicals, molecular_formula, molecular_bonds,
    double_bond_stereochemistry, tetrahedral_chirality, bond_lengths,
    bond_angles, internal_steric_clash, aromatic_ring_flatness,
    non-aromatic_ring_non-flatness, double_bond_flatness, internal_energy

Usage:

python3 find_all_true_posebusters.py [input_csv] [output_txt]

python3 find_all_true_posebusters.py combined_posebusters.csv passing_pdbs.txt

Defaults:
    input_csv  → combined_posebusters_results.csv
    output_txt → passing_pdbs.txt
"""

import sys
import pandas as pd
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
input_csv  = sys.argv[1] if len(sys.argv) > 1 else "combined_posebusters.csv"
output_txt = sys.argv[2] if len(sys.argv) > 2 else "passing_pdbs.txt"

# ── parameters to check (rmsd excluded) ─────────────────────────────────────
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

# ── load ─────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PoseBusters Pass-Filter")
print(f"{'='*60}")
print(f"Input  : {input_csv}")
print(f"Output : {output_txt}")
print(f"{'='*60}\n")

df = pd.read_csv(input_csv)
print(f"Total rows loaded : {len(df)}")

# ── verify all param columns exist ───────────────────────────────────────────
missing = [c for c in PARAMS if c not in df.columns]
if missing:
    print(f"\n[ERROR] Missing columns in CSV: {missing}")
    sys.exit(1)

# ── filter: all PARAMS must be True ──────────────────────────────────────────
mask = df[PARAMS].eq(True).all(axis=1)
passing = df[mask][["CCD_ID", "PDB_ID"]].drop_duplicates()

print(f"Rows passing all parameters: {len(passing)}\n")

# ── group by CCD_ID ──────────────────────────────────────────────────────────
grouped = passing.groupby("CCD_ID")["PDB_ID"].apply(list)

# ── print to terminal ─────────────────────────────────────────────────────────
print(f"{'─'*60}")
print(f"{'CCD ID':<12}  PDB IDs")
print(f"{'─'*60}")
for ccd_id, pdbs in grouped.items():
    pdbs_str = ", ".join(sorted(set(pdbs)))
    print(f"{str(ccd_id):<12}  {pdbs_str}")
print(f"{'─'*60}")
print(f"\nTotal CCD IDs with ≥1 passing PDB : {len(grouped)}")
print(f"Total passing PDB entries          : {sum(len(v) for v in grouped)}\n")

# ── write output txt ─────────────────────────────────────────────────────────
out_path = Path(output_txt)
with out_path.open("w") as fh:
    fh.write("PoseBusters – PDBs passing ALL parameters (rmsd excluded)\n")
    fh.write("="*60 + "\n\n")
    fh.write(f"{'CCD_ID':<12}  PDB_IDs\n")
    fh.write("-"*60 + "\n")
    for ccd_id, pdbs in grouped.items():
        pdbs_str = ", ".join(sorted(set(pdbs)))
        fh.write(f"{str(ccd_id):<12}  {pdbs_str}\n")
    fh.write("-"*60 + "\n")
    fh.write(f"\nTotal CCD IDs with ≥1 passing PDB : {len(grouped)}\n")
    fh.write(f"Total passing PDB entries          : {sum(len(v) for v in grouped)}\n")

print(f"Results saved to: {out_path.resolve()}")
print(f"{'='*60}\n")
