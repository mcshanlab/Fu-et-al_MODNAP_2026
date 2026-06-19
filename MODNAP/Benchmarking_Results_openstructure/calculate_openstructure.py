#!/usr/bin/env python3

# usage: python3 calculate_openstructure.py

import csv
import glob
import subprocess
from pathlib import Path
import json

# ============================================================
# PATHS
# ============================================================

AF3_BASE = Path("../MODNAP_AF3_models")
EXP_BASE = Path("../MODNAP")

OUT_CSV = Path("OST_lddt_TM_results.csv")
LOG_DIR = Path("ost_logs")
LOG_DIR.mkdir(exist_ok=True)

# ============================================================
# FILES
# ============================================================

cif_files = glob.glob(str(AF3_BASE / "*" / "*" / "af_output" / "*" / "*.cif"))
print(f"Found {len(cif_files)} CIF files")

# ============================================================
# CSV
# ============================================================

f_out = open(OUT_CSV, "w", newline="")
writer = csv.DictWriter(
    f_out,
    fieldnames=[
        "ccd",
        "pdb_id",
        "model_file",
        "ref_file",
        "lddt",
        "tm_score",
        "status"
    ]
)
writer.writeheader()

# ============================================================
# OST PARSER (UNCHANGED)
# ============================================================

def parse_ost(json_path):
    try:
        with open(json_path) as f:
            data = json.load(f)

        lddt = data.get("lddt") or data.get("scores", {}).get("lddt")
        tm = data.get("tm_score") or data.get("scores", {}).get("tm_score")

        return lddt, tm

    except Exception:
        return None, None


# ============================================================
# MAIN LOOP
# ============================================================

for cif in cif_files:

    cif = Path(cif)

    model_name = cif.stem
    if not model_name.endswith("_model"):
        continue

    base = model_name.replace("_model", "")
    ccd, pdb_id = base.split("_", 1)

    ref = EXP_BASE / ccd / pdb_id / f"{base}.pdb"

    if not ref.exists():
        writer.writerow({
            "ccd": ccd,
            "pdb_id": pdb_id,
            "model_file": str(cif),
            "ref_file": str(ref),
            "lddt": None,
            "tm_score": None,
            "status": "MISSING_REF"
        })
        continue

    # ========================================================
    # OUTPUT DIR
    # ========================================================

    out_dir = LOG_DIR / base
    out_dir.mkdir(parents=True, exist_ok=True)

    json_out = out_dir / "ost.json"
    log_out = out_dir / "ost.log"

    # ========================================================
    # OST (DO NOT TOUCH)
    # ========================================================

    cmd = [
        "ost", "compare-structures",
        "-m", str(cif),
        "-r", str(ref),
        "--lddt",
        "--tm-score",
        "--fault-tolerant",
        "-o", str(json_out)
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    log_out.write_text(proc.stdout)

    lddt, tm = parse_ost(json_out)
    
    if lddt is None:
        print(f"  [DEBUG] {base} returncode={proc.returncode}")
        print(f"  [DEBUG] {base} json_exists={json_out.exists()}")
        print(f"  [DEBUG] {base} stdout={proc.stdout[:200]}")
        print(f"  [DEBUG] {base} stderr={proc.stderr[:200]}")
        if json_out.exists():
            try:
                with open(json_out) as f:
                    raw = json.load(f)
                print(f"  [DEBUG] {base} keys: {list(raw.keys())}")
                print(f"  [DEBUG] {base} JSON: {json.dumps(raw, indent=2)[:300]}")
            except Exception as e:
                print(f"  [DEBUG] {base} JSON parse error: {e}")

    # ========================================================
    # WRITE OUTPUT (NO RMSD)
    # ========================================================

    writer.writerow({
        "ccd": ccd,
        "pdb_id": pdb_id,
        "model_file": str(cif),
        "ref_file": str(ref),
        "lddt": lddt,
        "tm_score": tm,
        "status": "OK" if proc.returncode == 0 else "OST_FAIL"
    })

    f_out.flush()

    print(f"[OK] {base} lDDT={lddt} TM={tm}")

# ============================================================
# DONE
# ============================================================

f_out.close()
print("DONE →", OUT_CSV)
