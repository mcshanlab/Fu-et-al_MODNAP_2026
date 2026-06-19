#!/usr/bin/env python3
"""
summarize_openstructure.py

Walk a tree of CCD -> PDB -> af_output to extract:
  - lDDT and TM-score (from OST_lddt_TM_results.csv)
  - pTM and ipTM (from af_output/*/*_summary_confidences.json)

Outputs printed table to terminal and writes a text file:
  ./openstructure_summary.txt   (in current working directory)

Usage:

python3 summarize_openstructure.py --base ../MODNAP_AF3_models --csv ./OST_lddt_TM_results.csv --out ./openstructure_summary.txt


"""
from pathlib import Path
import json
import re
import csv
import argparse
from datetime import datetime
import sys

FLOAT_RE = re.compile(r"([-+]?\d*\.\d+|\d+(\.\d*)?)([eE][-+]?\d+)?")


# ---------------------------------------------------------------------------
# CSV parsing
# ---------------------------------------------------------------------------

def load_ost_csv(csv_path: Path):
    """
    Parse OST_lddt_TM_results.csv and return a dict keyed by (ccd_upper, pdb_upper)
    with values {'lddt': str, 'tm': str, 'status': str}.
    Columns expected (case-insensitive): ccd, pdb_id, lddt, tm_score, status
    """
    results = {}
    try:
        with open(csv_path, newline='', encoding='utf-8') as fh:
            # detect delimiter — try tab first, then comma
            sample = fh.read(4096)
            fh.seek(0)
            dialect = 'excel-tab' if '\t' in sample else 'excel'
            reader = csv.DictReader(fh, dialect=dialect)
            # normalise header names
            for row in reader:
                norm = {k.strip().lower(): v.strip() for k, v in row.items() if k}
                ccd = norm.get('ccd', '').strip()
                pdb = norm.get('pdb_id', '').strip()
                if not ccd or not pdb:
                    continue
                lddt = norm.get('lddt')
                tm = norm.get('tm_score')
                status = norm.get('status', '')
                key = (ccd.upper(), pdb.upper())
                results[key] = {
                    'lddt': lddt if lddt not in (None, '') else None,
                    'tm':   tm   if tm   not in (None, '') else None,
                    'status': status,
                }
    except Exception as e:
        print(f"WARNING: could not parse CSV {csv_path}: {e}", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# AF3 summary confidences
# ---------------------------------------------------------------------------

def find_in_dict(d, key_patterns):
    """Recursively search dict/list for keys matching any compiled regex."""
    if isinstance(d, dict):
        for k, v in d.items():
            for patt in key_patterns:
                if patt.search(str(k)):
                    return v
            found = find_in_dict(v, key_patterns)
            if found is not None:
                return found
    elif isinstance(d, list):
        for item in d:
            found = find_in_dict(item, key_patterns)
            if found is not None:
                return found
    return None


def normalize_scalar(x):
    if x is None:
        return None
    if isinstance(x, (int, float, str)):
        return str(x)
    txt = json.dumps(x)
    m = FLOAT_RE.search(txt)
    return m.group(0) if m else str(x)


def find_af_summary_json(af_output_dir: Path, ccd: str, pdb: str):
    """
    Given an af_output directory, find *_summary_confidences.json.
    Search strategy:
      1. Subdir whose name contains both ccd and pdb (case-insensitive)
      2. Constructed name {ccd_low}_{pdb_low}
      3. Any *summary_confidences.json recursively
    """
    if not af_output_dir.exists() or not af_output_dir.is_dir():
        return None

    ccd_low = ccd.lower()
    pdb_low = pdb.lower()

    subdirs = [d for d in sorted(af_output_dir.iterdir()) if d.is_dir()]

    # 1) exact constructed name
    for cand_name in [f"{ccd_low}_{pdb_low}", f"{ccd}_{pdb}"]:
        candp = af_output_dir / cand_name
        if candp.is_dir():
            for match in candp.glob("*summary_confidences.json"):
                return match

    # 2) heuristic: subdir name contains both ccd and pdb
    for d in subdirs:
        dn = d.name.lower()
        if ccd_low in dn and pdb_low in dn:
            for match in d.glob("*summary_confidences.json"):
                return match

    # 3) recursive fallback
    for match in af_output_dir.rglob("*summary_confidences.json"):
        parent = match.parent.name.lower()
        if ccd_low in parent or pdb_low in parent:
            return match

    return None


def extract_ptm_iptm(json_path: Path):
    """Return (ptm, iptm, errstr) from a summary_confidences.json."""
    try:
        with open(json_path, 'r', encoding='utf-8') as fh:
            data = json.load(fh)
    except Exception as e:
        return (None, None, f"af_json_error:{e}")

    ptm_patts  = [re.compile(p, re.I) for p in [r"\bpTM\b",  r"\bptm\b",  r"ptm_score",  r"\bptm_confidence\b"]]
    iptm_patts = [re.compile(p, re.I) for p in [r"\bipTM\b", r"\biptm\b", r"iptm_score", r"\bip_tm\b"]]

    ptm  = normalize_scalar(find_in_dict(data, ptm_patts))
    iptm = normalize_scalar(find_in_dict(data, iptm_patts))
    return (ptm, iptm, None)


# ---------------------------------------------------------------------------
# Main walk
# ---------------------------------------------------------------------------

def walk_and_extract(base_dir: Path, ost_data: dict):
    """
    Walk base_dir where immediate children are CCD folders, their children are PDB folders.
    For each PDB folder:
      - look up lDDT and TM from ost_data
      - look up pTM / ipTM from af_output/*_summary_confidences.json
    """
    results = []

    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory does not exist: {base_dir}")
    if not base_dir.is_dir():
        raise NotADirectoryError(f"Base path is not a directory: {base_dir}")

    for ccd_dir in sorted(base_dir.iterdir()):
        if not ccd_dir.is_dir():
            continue
        for pdb_dir in sorted(ccd_dir.iterdir()):
            if not pdb_dir.is_dir():
                continue

            ccd_name = ccd_dir.name
            pdb_name = pdb_dir.name
            note_parts = []

            # lDDT + TM from CSV
            key = (ccd_name.upper(), pdb_name.upper())
            csv_row = ost_data.get(key, {})
            tm   = csv_row.get('tm')
            lddt = csv_row.get('lddt')
            if not csv_row:
                note_parts.append("csv: no match")
            elif csv_row.get('status', '').upper() not in ('', 'OK'):
                note_parts.append(f"csv_status:{csv_row['status']}")

            # pTM + ipTM from af_output
            af_output_dir = pdb_dir / "af_output"
            af_json_path  = find_af_summary_json(af_output_dir, ccd_name, pdb_name)
            if af_json_path is None:
                note_parts.append("af_output: summary_confidences.json none")
                ptm = iptm = None
            else:
                ptm, iptm, af_err = extract_ptm_iptm(af_json_path)
                if af_err:
                    note_parts.append(af_err)

            results.append({
                "ccd":  ccd_name,
                "pdb":  pdb_name,
                "tm":   tm,
                "lddt": lddt,
                "ptm":  ptm,
                "iptm": iptm,
                "note": " ; ".join(note_parts) if note_parts else None,
            })

    return results


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def format_table(results):
    cols   = ["CCD", "PDB", "TM", "LDDT", "pTM", "ipTM", "NOTE"]
    widths = {c: len(c) for c in cols}
    for r in results:
        widths["CCD"]  = max(widths["CCD"],  len(str(r["ccd"])))
        widths["PDB"]  = max(widths["PDB"],  len(str(r["pdb"])))
        widths["TM"]   = max(widths["TM"],   len(str(r["tm"]   or "none")))
        widths["LDDT"] = max(widths["LDDT"], len(str(r["lddt"] or "none")))
        widths["pTM"]  = max(widths["pTM"],  len(str(r["ptm"]  or "none")))
        widths["ipTM"] = max(widths["ipTM"], len(str(r["iptm"] or "none")))
        widths["NOTE"] = max(widths["NOTE"], len(str(r["note"] or "")))

    header = "  ".join(c.ljust(widths[c]) for c in cols)
    lines  = [header, "-" * len(header)]
    for r in results:
        line = "  ".join([
            str(r["ccd"]).ljust(widths["CCD"]),
            str(r["pdb"]).ljust(widths["PDB"]),
            str(r["tm"]   or "none").ljust(widths["TM"]),
            str(r["lddt"] or "none").ljust(widths["LDDT"]),
            str(r["ptm"]  or "none").ljust(widths["pTM"]),
            str(r["iptm"] or "none").ljust(widths["ipTM"]),
            str(r["note"] or "").ljust(widths["NOTE"]),
        ])
        lines.append(line)
    return "\n".join(lines)


def write_summary_file(out_path: Path, results):
    now    = datetime.now().isoformat(sep=' ', timespec='seconds')
    header = f"OpenStructure summary generated: {now}\nBase scan produced {len(results)} PDB entries.\n\n"
    table  = format_table(results)
    with out_path.open("w", encoding='utf-8') as fh:
        fh.write(header)
        fh.write(table)
        fh.write("\n\n# End of report\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    default_base = Path.home() / "Desktop" / "ML_aptamer_project" / "MODNAP_AF3_models"
    default_csv  = Path.home() / "Desktop" / "ML_aptamer_project" / "OST_lddt_TM_results.csv"

    parser = argparse.ArgumentParser(
        description="Summarize lDDT/TM from OST CSV and pTM/ipTM from AF3 summary_confidences.json."
    )
    parser.add_argument("--base", "-b", default=str(default_base),
                        help=f"Base AF3 models folder (default: {default_base})")
    parser.add_argument("--csv",  "-c", default=str(default_csv),
                        help=f"OST lDDT/TM CSV file (default: {default_csv})")
    parser.add_argument("--out",  "-o", default="openstructure_summary.txt",
                        help="Output summary text file (default: ./openstructure_summary.txt)")
    args = parser.parse_args()

    base_dir = Path(args.base).expanduser().resolve()
    csv_path = Path(args.csv).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()

    # Load CSV data
    if not csv_path.exists():
        print(f"WARNING: CSV file not found: {csv_path}", file=sys.stderr)
        ost_data = {}
    else:
        ost_data = load_ost_csv(csv_path)
        print(f"Loaded {len(ost_data)} rows from CSV: {csv_path}")

    # Walk directory tree
    try:
        results = walk_and_extract(base_dir, ost_data)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    if not results:
        print(f"No PDB entries found under: {base_dir}")
    else:
        print(format_table(results))

    # Write output file
    try:
        write_summary_file(out_path, results)
        print(f"\nResults written to: {out_path}")
    except Exception as e:
        print(f"Failed to write output file {out_path}: {e}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
