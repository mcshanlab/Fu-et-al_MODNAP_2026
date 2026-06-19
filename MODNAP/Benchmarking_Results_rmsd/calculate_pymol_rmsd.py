#!/usr/bin/env python3

"""
usage: python3 calculate_pymol_rmsd.py

Calculates RMSD between AF3 CIF models and experimental PDBs using PyMOL in headless mode.

Outputs:
pymol_rmsd_summary.txt

Format:
CCD    PDB    RMSD    N_atoms
"""

import glob
from pathlib import Path
import sys

# ============================================================
# PYMOL HEADLESS
# ============================================================

try:
    import pymol
    from pymol import cmd

    pymol.finish_launching(["pymol", "-cq"])

except Exception as e:
    print("ERROR: Failed to launch PyMOL headless")
    print(e)
    sys.exit(1)

# ============================================================
# PATHS
# ============================================================

AF3_BASE = Path("../MODNAP_AF3_models")
EXP_BASE = Path("../MODNAP")

OUT_TXT = Path("pymol_rmsd_summary.txt")

# ============================================================
# FIND CIF FILES
# ============================================================

cif_files = glob.glob(
    str(AF3_BASE / "*" / "*" / "af_output" / "*" / "*.cif")
)

print(f"Found {len(cif_files)} CIF files")

results = []

# ============================================================
# MAIN LOOP
# ============================================================

for cif in cif_files:

    cif = Path(cif)

    model_name = cif.stem

    # Require *_model.cif
    if not model_name.endswith("_model"):
        continue

    base = model_name.replace("_model", "")

    try:
        ccd, pdb_id = base.split("_", 1)

    except Exception:

        print(f"[SKIP] could not parse {model_name}")

        results.append((
            "NA",
            model_name,
            "NA",
            "NA"
        ))

        continue

    # ========================================================
    # REFERENCE PDB
    # ========================================================

    ref = (
        EXP_BASE
        / ccd.upper()
        / pdb_id.upper()
        / f"{ccd.upper()}_{pdb_id.upper()}.pdb"
    )

    if not ref.exists():

        print(f"[MISSING_REF] {base}")

        results.append((
            ccd.upper(),
            pdb_id.upper(),
            "NA",
            "NA"
        ))

        continue

    print(f"[RUNNING] {base}")

    try:

        # ====================================================
        # CLEAN SESSION
        # ====================================================

        cmd.reinitialize()

        # avoid reserved keywords
        ref_obj = "ref_obj"
        model_obj = "pred_obj"

        # ====================================================
        # LOAD STRUCTURES
        # ====================================================

        cmd.load(str(ref), ref_obj)
        cmd.load(str(cif), model_obj)

        # ====================================================
        # SANITY CHECK
        # ====================================================

        n_ref = cmd.count_atoms(ref_obj)
        n_model = cmd.count_atoms(model_obj)

        if n_ref == 0 or n_model == 0:

            print(f"[ZERO_ATOMS] {base}")

            results.append((
                ccd.upper(),
                pdb_id.upper(),
                "NA",
                "NA"
            ))

            continue

        # ====================================================
        # ALIGN + RMSD
        # ====================================================

        aln = cmd.align(model_obj, ref_obj)

        rmsd = aln[0]
        n_atoms = aln[1]

        # ====================================================
        # STORE RESULT
        # ====================================================

        results.append((
            ccd.upper(),
            pdb_id.upper(),
            rmsd,
            n_atoms
        ))

        print(
            f"[OK] {base} "
            f"RMSD={rmsd:.4f} "
            f"aligned_atoms={n_atoms}"
        )

    except Exception as e:

        print(f"[FAIL] {base}")
        print(e)

        results.append((
            ccd.upper(),
            pdb_id.upper(),
            "NA",
            "NA"
        ))

# ============================================================
# SORT RESULTS (MINIMAL ADDITION)
# ============================================================

results.sort(key=lambda x: (x[0], x[1]))

# ============================================================
# WRITE OUTPUT
# ============================================================

with OUT_TXT.open("w") as fh:

    fh.write("CCD\tPDB\tRMSD\tN_atoms\n")

    for ccd, pdb, rmsd, n_atoms in results:

        if isinstance(rmsd, float):

            fh.write(
                f"{ccd}\t{pdb}\t{rmsd:.6f}\t{n_atoms}\n"
            )

        else:

            fh.write(
                f"{ccd}\t{pdb}\t{rmsd}\t{n_atoms}\n"
            )

# ============================================================
# DONE
# ============================================================

print(f"\nDONE → {OUT_TXT}")
