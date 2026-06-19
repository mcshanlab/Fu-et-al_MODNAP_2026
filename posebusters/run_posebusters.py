#!/usr/bin/env python3
"""
python3 run_posebusters.py

For each PDB folder, this script:
  1) loads the AF3 predicted CIF in PyMOL
  2) extracts only the CCD residue from each structure
  3) saves:
       PDBID_CCDID_model.pdb
     inside the same PDB folder
  4) runs PoseBusters:
       python -m posebusters -l PDBID_CCDID.pdb --outfmt csv --output pred_posebusters_report.csv
  5) saves .csv inside the same PDB folder

Requirements:
  - PyMOL executable available as `pymol` in PATH
  - PoseBusters installed in the same Python environment as this script
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path


DEFAULT_BASE_DIR = Path("/Volumes/GigiMurin/Clary_ML_paper/step19_posebusters/openstructure_last")


def strip_ok_suffix(name: str) -> str:
    """Remove a trailing _ok / _OK / _Ok / etc. if present."""
    lower = name.lower()
    if lower.endswith("_ok"):
        return name[:-3]
    return name


def find_case_insensitive_exact_file(folder: Path, filename: str) -> Path | None:
    """Find a file in folder matching filename case-insensitively."""
    target = filename.lower()
    for p in folder.iterdir():
        if p.is_file() and p.name.lower() == target:
            return p
    return None


def run_pymol_extract(
    pred_cif: Path,
    ccd_id: str,
    out_pred_pdb: Path,
) -> None:
    """
    Use headless PyMOL to extract the residue named ccd_id from both inputs.
    """
    pymol_exe = shutil.which("pymol")
    if not pymol_exe:
        raise RuntimeError(
            "Could not find the PyMOL executable 'pymol' in PATH. "
            "Please activate the environment that has PyMOL installed."
        )

    
    ccd_sel = ccd_id.upper()

    pml = textwrap.dedent(
        f"""
        reinitialize

        load {pred_cif.as_posix()}, pred
        select mod_pred, resn {ccd_sel}
        save {out_pred_pdb.as_posix()}, mod_pred, format=pdb

        quit
        """
    ).strip() + "\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        pml_path = Path(tmpdir) / "extract_mod.pml"
        pml_path.write_text(pml, encoding="utf-8")

        subprocess.run(
            [pymol_exe, "-cq", str(pml_path)],
            check=True,
        )

    # Verify files were actually created
    if not out_pred_pdb.exists() or out_pred_pdb.stat().st_size == 0:
        raise RuntimeError(f"PyMOL failed to create {out_pred_pdb}")

def run_posebusters(pred_pdb: Path, out_csv: Path) -> None:
    """
    Run PoseBusters on a predicted model only.
    """
    cmd = [
        sys.executable,
        "-m",
        "posebusters",
        str(pred_pdb),
        "--outfmt",
        "csv",
        "--output",
        str(out_csv),
    ]
    subprocess.run(cmd, check=True)


def process_one_pdb_folder(pdb_dir: Path) -> None:
    """
    Process one PDB folder.
    """
    ccd_dir = pdb_dir.parent
    ccd_id = strip_ok_suffix(ccd_dir.name)
    pdb_id = strip_ok_suffix(pdb_dir.name)

    # Expected predicted file: af_output/ccdid_pdbid/ccdid_pdbid_model.cif
    model_rel_dir = Path("af_output") / f"{ccd_id.lower()}_{pdb_id.lower()}"
    pred_cif = pdb_dir / model_rel_dir / f"{ccd_id.lower()}_{pdb_id.lower()}_model.cif"
    if not pred_cif.exists():
        raise FileNotFoundError(f"Could not find predicted CIF file {pred_cif}")

    run_dir = Path.cwd() / "posebusters_outputs" / ccd_id / pdb_id
    run_dir.mkdir(parents=True, exist_ok=True)

    out_pred_pdb = run_dir / f"{pdb_id}_{ccd_id}_model.pdb"
    out_csv = run_dir / "pred_posebusters_report.csv"

    print(f"\n[INFO] Processing: {ccd_id} / {pdb_id}")
    print(f"       predicted:    {pred_cif}")

    run_pymol_extract(pred_cif, ccd_id, out_pred_pdb)
    print(f"[OK]   Saved extracted files:")
    print(f"       {out_pred_pdb.name}")

    run_posebusters(out_pred_pdb, out_csv)
    print(f"[OK]   Saved PoseBusters CSV: {out_csv.name}")


def iter_pdb_folders(base_dir: Path):
    """
    Yield likely PDB folders under CCD folders.
    """
    for ccd_dir in sorted(base_dir.iterdir()):
        if not ccd_dir.is_dir():
            continue
        if ccd_dir.name.startswith("."):
            continue

        for pdb_dir in sorted(ccd_dir.iterdir()):
            if not pdb_dir.is_dir():
                continue
            if pdb_dir.name.startswith("."):
                continue
            # Only consider folders that appear to be PDB folders.
            yield pdb_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-extract CCD residues with PyMOL and run PoseBusters CSV reports."
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help=f"Root directory to scan (default: {DEFAULT_BASE_DIR})",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop immediately if one folder fails instead of continuing.",
    )
    args = parser.parse_args()

    base_dir: Path = args.base_dir
    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory does not exist: {base_dir}")

    print(f"[INFO] Scanning: {base_dir}")

    n_total = 0
    n_ok = 0
    n_fail = 0

    for pdb_dir in iter_pdb_folders(base_dir):
        n_total += 1
        try:
            process_one_pdb_folder(pdb_dir)
            n_ok += 1
        except Exception as exc:
            n_fail += 1
            print(f"[ERROR] {pdb_dir}: {exc}")
            if args.stop_on_error:
                raise

    print("\n[DONE]")
    print(f"  total folders seen: {n_total}")
    print(f"  successful:         {n_ok}")
    print(f"  failed:              {n_fail}")


if __name__ == "__main__":
    main()
