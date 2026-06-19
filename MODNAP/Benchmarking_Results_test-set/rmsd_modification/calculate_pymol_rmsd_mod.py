#!/usr/bin/env python3

"""
usage: python3 calculate_pymol_rmsd_mod.py

Calculates full-structure RMSD and modification (CCD ligand) RMSD
between AF3 CIF models and experimental PDBs using PyMOL in headless mode.

Outputs:
pymol_rmsd_summary_mod.txt

Format:
CCD    PDB    RMSD    N_atoms    Mod_RMSD    Mod_N_atoms
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

AF3_BASE = Path("../../MODNAP_AF3_models")
EXP_BASE = Path("../../MODNAP")

OUT_TXT = Path("pymol_rmsd_summary_mod.txt")

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
            "NA", model_name, "NA", "NA", "NA", "NA"
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
            ccd.upper(), pdb_id.upper(), "NA", "NA", "NA", "NA"
        ))

        continue

    print(f"[RUNNING] {base}")

    try:

        # ====================================================
        # CLEAN SESSION
        # ====================================================

        cmd.reinitialize()

        ref_obj   = "ref_obj"
        model_obj = "pred_obj"

        # ====================================================
        # LOAD STRUCTURES
        # ====================================================

        cmd.load(str(ref), ref_obj)
        cmd.load(str(cif), model_obj)

        # ====================================================
        # SANITY CHECK
        # ====================================================

        n_ref   = cmd.count_atoms(ref_obj)
        n_model = cmd.count_atoms(model_obj)

        if n_ref == 0 or n_model == 0:

            print(f"[ZERO_ATOMS] {base}")

            results.append((
                ccd.upper(), pdb_id.upper(), "NA", "NA", "NA", "NA"
            ))

            continue

        # ====================================================
        # FULL STRUCTURE RMSD via cmd.align
        # (handles atom count mismatches via sequence alignment)
        # model_obj is superimposed onto ref_obj after this call
        # ====================================================

        aln    = cmd.align(model_obj, ref_obj)
        rmsd   = aln[0]
        n_atoms = aln[1]

        # ====================================================
        # MOD RMSD: CCD ligand only
        # Use first matching residue in each structure
        # rms_cur with matchmaker=4: sequence-based atom matching,
        # no extra fitting (structures already aligned above)
        # ====================================================

        sel_ref = f"{ref_obj} and resn {ccd}"
        sel_tgt = f"{model_obj} and resn {ccd}"

        if cmd.count_atoms(sel_ref) == 0 or cmd.count_atoms(sel_tgt) == 0:

            print(f"[NO_LIGAND] {base} — resn {ccd} not found in one or both structures")

            results.append((
                ccd.upper(), pdb_id.upper(), rmsd, n_atoms, "NA", "NA"
            ))

            continue

        resi_ref = cmd.get_model(sel_ref).atom[0].resi
        resi_tgt = cmd.get_model(sel_tgt).atom[0].resi

        cmd.select("ref_lig", f"{ref_obj} and resn {ccd} and resi {resi_ref}")
        cmd.select("tgt_lig", f"{model_obj} and resn {ccd} and resi {resi_tgt}")

        n_mod    = cmd.count_atoms("ref_lig")
        rmsd_mod = cmd.rms_cur("tgt_lig", "ref_lig", matchmaker=4)

        # ====================================================
        # STORE RESULT
        # ====================================================

        results.append((
            ccd.upper(), pdb_id.upper(), rmsd, n_atoms, rmsd_mod, n_mod
        ))

        print(
            f"[OK] {base} "
            f"RMSD={rmsd:.4f} aligned_atoms={n_atoms} "
            f"Mod_RMSD={rmsd_mod:.4f} Mod_N_atoms={n_mod}"
        )

    except Exception as e:

        print(f"[FAIL] {base}")
        print(e)

        results.append((
            ccd.upper(), pdb_id.upper(), "NA", "NA", "NA", "NA"
        ))

# ============================================================
# SORT RESULTS
# ============================================================

results.sort(key=lambda x: (x[0], x[1]))

# ============================================================
# WRITE OUTPUT
# ============================================================

with OUT_TXT.open("w") as fh:

    fh.write("CCD\tPDB\tRMSD\tN_atoms\tMod_RMSD\tMod_N_atoms\n")

    for ccd, pdb, rmsd, n_atoms, rmsd_mod, n_mod in results:

        rmsd_str    = f"{rmsd:.6f}"    if isinstance(rmsd, float)    else rmsd
        mod_str     = f"{rmsd_mod:.6f}" if isinstance(rmsd_mod, float) else rmsd_mod
        n_mod_str   = str(n_mod)       if isinstance(n_mod, int)      else n_mod

        fh.write(f"{ccd}\t{pdb}\t{rmsd_str}\t{n_atoms}\t{mod_str}\t{n_mod_str}\n")

# ============================================================
# DONE
# ============================================================

print(f"\nDONE → {OUT_TXT}")
