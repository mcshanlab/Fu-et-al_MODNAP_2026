#!/usr/bin/env python3
"""
usage: python3 plot_openstructure_cdfs.py --outdir .

"""
from pathlib import Path
import argparse
import re
import csv
import math
import numpy as np
import matplotlib
import matplotlib.font_manager as fm
import json

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['pdf.use14corefonts'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'

arial_exists = any("Arial" in f for f in fm.findSystemFonts())
matplotlib.rcParams['font.family'] = 'Arial' if arial_exists else 'DejaVu Sans'

import matplotlib.pyplot as plt

SPLIT_RE = re.compile(r"\s{2,}")  # split on two or more spaces

def parse_summary_file(path: Path):
    text = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    header_idx = None
    for i, line in enumerate(text):
        if "CCD" in line and "PDB" in line and "TM" in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Could not find table header in {path!s}")
    data_start = header_idx + 1
    if data_start < len(text) and set(text[data_start].strip()) <= set("- "):
        data_start += 1

    records = []
    for line in text[data_start:]:
        if not line.strip():
            continue
        if line.startswith("# End of report"):
            break
        parts = SPLIT_RE.split(line.strip())
        if len(parts) < 6:
            parts += [""] * (6 - len(parts))
        ccd = parts[0].strip()
        pdb = parts[1].strip() if len(parts) > 1 else ""
        tm = parts[2].strip() if len(parts) > 2 else ""
        # extract LDDT (new column in the summary), then shift pTM/iPTM/RMSD accordingly
        lddt = parts[3].strip() if len(parts) > 3 else ""
        ptm = parts[4].strip() if len(parts) > 4 else ""
        iptm = parts[5].strip() if len(parts) > 5 else ""
        rmsd = parts[6].strip() if len(parts) > 6 else ""
        note = parts[7].strip() if len(parts) > 7 else ""
        records.append({
            "ccd": ccd,
            "pdb": pdb,
            "tm": tm,
            "lddt": lddt,
            "ptm": ptm,
            "iptm": iptm,
            "rmsd": rmsd,
            "note": note
        })
    return records

import csv

def parse_inf_csv(path: Path):
    """
    Parse results_inf.csv and return:
        records = [
            {
                "ccd": ...,
                "pdb": ...,
                "inf_wc": float,
                "inf_all": float
            },
            ...
        ]
    """
    records = []

    if not path.exists():
        print(f"[WARN] INF file not found: {path}")
        return records

    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)

        for row in reader:

            inf_wc = to_float(row.get("INF_WC"))
            inf_all = to_float(row.get("INF_all"))

            records.append({
                "ccd": row["CCD"].strip().upper(),
                "pdb": row["PDB"].strip().upper(),
                "inf_wc": inf_wc,
                "inf_all": inf_all,
            })

    return records

def to_float(x):
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    if s.lower().startswith("none"):
        return None
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    m = float_re.search(s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except Exception:
        return None

def parse_pymol_rmsd(path: Path):
    """Read RMSD values from a pymol_rmsd_summary.txt TSV-like file.

    Expected header: CCD\tPDB\tRMSD\tN_atoms
    Lines: CCD\tPDB\t<RMSD>\t<N_atoms>
    Returns: list of floats (or empty list if none) — kept for compatibility
    with your original workflow.
    """
    if not path.exists():
        print(f"[ERROR] pymol RMSD file not found: {path}")
        return None
    text = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    vals = []
    for line in text:
        line = line.strip()
        if not line:
            continue
        # skip header if present
        if line.lower().startswith('ccd') or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        rmsd_str = parts[2]
        v = to_float(rmsd_str)
        if v is not None:
            vals.append(v)
    return vals

def parse_pymol_rmsd_map(path: Path):
    """Parse pymol_rmsd_summary.txt into a mapping (ccd,pdb) -> rmsd_float.
    This is used to build category-specific RMSD lists.
    """
    mapping = {}
    if not path.exists():
        return mapping
    text = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    for line in text:
        line = line.strip()
        if not line:
            continue
        if line.lower().startswith('ccd') or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        ccd = parts[0].strip()
        if ccd.lower().endswith("_ok"):
            ccd = ccd[:-3]
        ccd = ccd.upper()
        pdb = parts[1].strip().upper()
        v = to_float(parts[2])
        if v is None:
            continue
        mapping[(ccd, pdb)] = v
    return mapping

def parse_nuc_type_file(path: Path):
    """
    Parse nuc_type_count_per_pdb.txt into a mapping (ccd,pdb) -> category_label.
    Expected to be tab-separated with header containing 'CCD' and 'PDB' and
    a 'Category' column.
    """
    mapping = {}
    if not path.exists():
        print(f"[WARN] nuc_type_count_per_pdb.txt not found: {path} (category-specific curves will be skipped)")
        return mapping

    text = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    if not text:
        return mapping

    # find header line
    header_idx = 0
    header = text[0]
    # detect delimiter (tabs preferred)
    delim = '\t' if '\t' in header else None
    cols = header.split(delim) if delim else header.split()
    # find column indices
    try:
        idx_ccd = cols.index("CCD")
        idx_pdb = cols.index("PDB")
        idx_category = cols.index("Category")
    except ValueError:
        # fallback: try lowercase
        lower = [c.lower() for c in cols]
        try:
            idx_ccd = lower.index("ccd")
            idx_pdb = lower.index("pdb")
            idx_category = lower.index("category")
        except ValueError:
            # can't parse header: give up gracefully
            print(f"[WARN] nuc_type_count_per_pdb.txt header didn't contain CCD/PDB/Category: {path}")
            return mapping

    for line in text[1:]:
        if not line.strip():
            continue
        parts = line.split(delim) if delim else line.split()
        # guard against short lines
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

def parse_passing_pdbs(path: Path):
    """
    Return a set of normalized (CCD, PDB) pairs.

    Handles mixed naming styles:
        7Q7X_ok -> 7Q7X
        3CVU    -> 3CVU

    Output:
        {("1MA", "7Q7X"), ("64T", "3CVU"), ...}
    """

    pairs = set()

    if not path.exists():
        print(f"[WARN] passing_pdbs.txt not found: {path}")
        return pairs

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():

        line = line.strip()

        # Skip headers / separators
        if (
            not line
            or line.startswith("CCD_ID")
            or line.startswith("-")
            or line.startswith("=")
            or line.startswith("Total")
            or line.startswith("Pose")
            or line.startswith("====")
            or "PDB_IDs" in line
        ):
            continue

        parts = line.split(None, 1)

        if len(parts) < 2:
            continue

        # Normalize CCD
        ccd = parts[0].strip().upper()

        pdb_list = parts[1]

        for pdb in pdb_list.split(","):

            pdb = pdb.strip()

            if not pdb:
                continue

            # --- normalize PDB name ---
            pdb = pdb.upper()

            if pdb.endswith("_OK"):
                pdb = pdb[:-3]

            # keep only first 4 characters (true PDB ID)
            pdb = pdb[:4]

            pairs.add((ccd, pdb))

    print("Parsed passing pairs:", len(pairs))

    return pairs

def load_json(path: Path):
    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        try:
            txt = path.read_text(errors="ignore")
            return json.loads(txt)
        except Exception:
            return None

def make_and_save_success_rate_plot(rmsd_groups, success_rmsd_groups, pb_rates, out_path: Path, cutoffs=(3.0, 2.5, 2.0)):
    categories = [
        "DNA",
        "RNA",
        "DNA-protein",
        "RNA-protein",
        "DNA-RNA",
        "DNA-RNA-protein",
    ]

    t2 = 2.0
    t25 = 2.5
    t3 = 3.0

    x = np.arange(len(categories))
    width = 0.6

    plt.figure(figsize=(8.8, 4.8))

    labels = [
        f"RMSD ≤ {t2:g} Å",
        f"{t2:g} < RMSD ≤ {t25:g} Å",
        f"{t25:g} < RMSD ≤ {t3:g} Å",
    ]

    colors = [
        "#E6EAFB",
        "#D1E6D5",
        "#FFE0B8",
    ]

    for i, cat in enumerate(categories):
        vals = rmsd_groups.get(cat, [])       # ALL PDBs in this category
        total = len(vals)

        if total == 0:
            continue

        n_le2  = sum(1 for v in vals if v <= t2)
        n_le25 = sum(1 for v in vals if v <= t25)
        n_le3  = sum(1 for v in vals if v <= t3)

        seg1 = (n_le2 / total) * 100.0
        seg2 = max(n_le25 - n_le2,  0) / total * 100.0
        seg3 = max(n_le3  - n_le25, 0) / total * 100.0

        # stacked bar
        plt.bar(x[i], seg1, width=width, color=colors[0], label=labels[0] if i == 0 else None)
        plt.bar(x[i], seg2, width=width, bottom=seg1, color=colors[1], label=labels[1] if i == 0 else None)
        plt.bar(x[i], seg3, width=width, bottom=seg1 + seg2, color=colors[2], label=labels[2] if i == 0 else None)

        # count label: n_in_category / total_all
        total_all = len(rmsd_groups.get("All", []))
        plt.text(
            x[i], -6,
            f"{total}/{total_all}",
            ha="center", va="top", fontsize=9,
        )
        
        # PB-valid hatched overlay
        plt.bar(
            x[i], pb_rates.get(cat, 0),
            width=width,
            facecolor="none", edgecolor="black",
            hatch="///", linewidth=1.2,
            label="PB-valid" if i == 0 else None,
        )

    plt.xticks(x, categories, rotation=0, ha="center")
    plt.ylabel("Success rate (%)")
    plt.xlabel("PDB type")
    plt.title("RMSD success rate by PDB type")
    plt.ylim(0, 105)
    plt.grid(True, axis="y", linestyle=":", linewidth=0.5)
    plt.legend(title="Cutoffs / PB-valid")
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=200)
    plt.savefig(str(out_path.with_suffix(".pdf")), dpi=200)
    plt.close()
    
# --- write summary table ---
    total_all = len(rmsd_groups.get("All", []))
    txt_path = out_path.with_suffix(".txt")
    with txt_path.open("w", encoding="utf-8") as fh:
        fh.write(f"Molecule Type\tMODNAP (n = {total_all})\t\n")
        fh.write(f"\t% Success Rate at 2.5 Å\t% PB-Valid\n")
        for cat in categories:
            vals = rmsd_groups.get(cat, [])
            total = len(vals)
            if total == 0:
                fh.write(f"{cat} (n={total})\tN/A\tN/A\n")
                continue
            n_le25 = sum(1 for v in vals if v <= t25)
            sr = n_le25 / total * 100.0
            pb = pb_rates.get(cat, 0.0)
            fh.write(f"{cat} (n={total})\t{sr:.1f}\t{pb:.1f}\n")

    return out_path          # <-- this was missing

def ecdf_on_grid(values, grid):
    """
    Given a list of numeric values and a sorted grid, return ECDF(grid) values:
      ECDF(x) = fraction of values <= x
    This is done efficiently via np.searchsorted on sorted values.
    """
    if len(values) == 0:
        return np.zeros_like(grid, dtype=float)
    vals = np.sort(np.array(values))
    # for each x in grid count how many vals <= x (use 'right')
    counts = np.searchsorted(vals, grid, side='right')
    return counts / float(len(vals))

def make_and_save_cdf(x_vals, x_name, out_path: Path, grid_step=None):
    """
    Original single-line function — kept for compatibility (unchanged).
    """
    x_vals = [v for v in x_vals if v is not None]
    if not x_vals:
        print(f"[WARN] No numeric values for {x_name}; skipping plot.")
        return None

    # Determine grid
    min_val = float(min(x_vals))
    max_val = float(max(x_vals))

    if x_name.lower() == "rmsd":
        step = grid_step if grid_step is not None else 0.1
        x_min = 0.0
        x_max = min(8.0, max(math.ceil(max_val / step) * step, max_val))
        grid = np.arange(x_min, x_max + step/2.0, step)
    else:
        step = grid_step if grid_step is not None else 0.001
        x_min = 0.0
        x_max = 1.0
        grid = np.arange(x_min, x_max + step/2.0, step)

    cdf_vals = ecdf_on_grid(x_vals, grid)
    if x_name.lower() != "rmsd":
        cdf_vals = 1.0 - cdf_vals

    plt.figure(figsize=(6.5, 4.0))
    plt.plot(grid, cdf_vals, lw=2)
    plt.xlabel(x_name)
    plt.ylabel("Cumulative frequency")
    plt.title(f"Cumulative frequency — {x_name}")
    plt.ylim(-0.02, 1.02)
    plt.xlim(grid[0], grid[-1])
    plt.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=200)
    pdf_path = out_path.with_suffix('.pdf')
    plt.savefig(str(pdf_path), dpi=200)
    plt.close()
    return out_path

def make_and_save_cdf_multi(values_dict, x_name, out_path: Path, grid_step=None):
    """
    New: Plot multiple ECDF-derived CDF lines on the same axes.
    values_dict: dict label -> list of numeric values (floats)
      - include an 'All' label (or any label) to represent the overall dataset.
    x_name: name (e.g., "pTM", "ipTM", "TM score", "LDDT", "RMSD")
    """
    # filter empty lists but keep track for legend
    filtered = {lab: [v for v in vals if v is not None] for lab, vals in values_dict.items()}
    # if none have any data, skip
    if not any(filtered.values()):
        print(f"[WARN] No numeric values for {x_name} across all groups; skipping plot.")
        return None

    # Determine grid (use union-range from all data for RMSD/others)
    all_vals = []
    for vals in filtered.values():
        all_vals.extend(vals)
    if not all_vals:
        print(f"[WARN] No numeric values for {x_name}; skipping plot.")
        return None

    min_val = float(min(all_vals))
    max_val = float(max(all_vals))

    if x_name.lower() == "rmsd":
        step = grid_step if grid_step is not None else 0.1
        x_min = 0.0
        # cap at 8.0 for plotting as in original
        x_max = 8.0
        grid = np.arange(x_min, x_max + step/2.0, step)
    else:
        step = grid_step if grid_step is not None else 0.001
        x_min = 0.0
        x_max = 1.0
        grid = np.arange(x_min, x_max + step/2.0, step)

    plt.figure(figsize=(6.5, 4.0))

    # make sure 'All' is plotted first and maybe emphasized
    labels_order = []
    if "All" in filtered:
        labels_order.append("All")
    # other canonical categories in a stable order
    for lab in ("DNA", "RNA", "DNA-protein", "RNA-protein"):
        if lab in filtered and lab != "All":
            labels_order.append(lab)
    # add any remaining labels
    for lab in filtered:
        if lab not in labels_order:
            labels_order.append(lab)

    # color/linestyle cycle will be default; make 'All' line thicker
    for lab in labels_order:
        vals = filtered.get(lab, [])
        if not vals:
            continue
        cdf_vals = ecdf_on_grid(vals, grid)
        if x_name.lower() != "rmsd":
            cdf_vals = 1.0 - cdf_vals
        lw = 2.5 if lab == "All" else 1.5
        alpha = 0.9 if lab == "All" else 0.8
        plt.plot(grid, cdf_vals, lw=lw, alpha=alpha, label=lab)

    plt.xlabel(x_name)
    plt.ylabel("Cumulative frequency")
    plt.title(f"Cumulative frequency — {x_name}")
    plt.ylim(-0.02, 1.02)
    plt.xlim(grid[0], grid[-1])
    plt.grid(True, linestyle=':', linewidth=0.5)
    plt.legend(title="Group", fontsize=8)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=200)
    pdf_path = out_path.with_suffix('.pdf')
    plt.savefig(str(pdf_path), dpi=200)
    plt.close()
    return out_path

def strip_ok(s: str):
    """Helper: remove trailing '_ok' if present (for flexible matching)."""
    if not isinstance(s, str):
        return s
    return s[:-3] if s.endswith("_ok") else s

def find_category_for_record(rec, category_map):
    """
    Try several lookups to match the (ccd,pdb) key in category_map.
    This is a bit permissive to tolerate '_ok' suffix differences.
    """
    ccd = rec.get("ccd", "")
    if ccd.lower().endswith("_ok"):
        ccd = ccd[:-3]
    ccd = ccd.upper()
    pdb = rec.get("pdb", "").upper()
    candidates = [
        (ccd, pdb),
        (strip_ok(ccd), pdb),
        (ccd, strip_ok(pdb)),
        (strip_ok(ccd), strip_ok(pdb)),
    ]
    for k in candidates:
        if k in category_map:
            return category_map[k]
    return None

def collect_metric_values_by_category(records, metric_key, category_map, rmsd_map=None):
    """
    Given the openstructure records (list of dicts) collect metric values
    grouped by category label (as found in category_map).
    Returns: dict category_label -> list of floats
    """
    groups = {}
    for rec in records:
        ccd = rec.get("ccd", "")
        pdb = rec.get("pdb", "")
        cat = find_category_for_record(rec, category_map)
        if cat is None:
            # skip entries without a category in mapping
            continue
        if metric_key == "rmsd":
            if rmsd_map is None:
                continue
            # try exact and permissive lookups
            v = None
            for k in [(ccd, pdb), (strip_ok(ccd), pdb), (ccd, strip_ok(pdb)), (strip_ok(ccd), strip_ok(pdb))]:
                if k in rmsd_map:
                    v = rmsd_map[k]
                    break
            if v is None:
                continue
        else:
            v = to_float(rec.get(metric_key, None))
        if v is None:
            continue
        groups.setdefault(cat, []).append(v)
    return groups

def main():
    default_input = Path("../Benchmarking_Results_openstructure/openstructure_summary.txt")
    p = argparse.ArgumentParser(description="Plot smooth ECDF (CDF) for PTM, iPTM, RMSD, TM from openstructure_summary.txt")
    p.add_argument("--input", "-i", default=str(default_input), help="Path to openstructure_summary.txt")
    p.add_argument("--outdir", "-o", default=None, help="Output directory for PNG files (default: same dir as input)")
    args = p.parse_args()
    made = []   # <-- MUST be here, unconditional

    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        return

    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else input_path.parent
    records = parse_summary_file(input_path)
    inf_path = input_path.parent / "../Benchmarking_Results_inf/results_inf.csv"
    inf_records = parse_inf_csv(inf_path)

    tm_vals = [to_float(r["tm"]) for r in records]
    ptm_vals = [to_float(r["ptm"]) for r in records]
    iptm_vals = [to_float(r["iptm"]) for r in records]
    # --- RMSD: read from pymol_rmsd_summary.txt in the same directory as the input ---
    pymol_path = input_path.parent / "../Benchmarking_Results_rmsd/pymol_rmsd_summary.txt"
    pymol_rmsd_vals = parse_pymol_rmsd(pymol_path)
    if pymol_rmsd_vals is None:
        print("[ERROR] Could not read pymol_rmsd_summary.txt; aborting RMSD plotting.")
        return
    rmsd_vals = pymol_rmsd_vals

    lddt_vals = [to_float(r.get("lddt", None)) for r in records]

    # filter numeric only; RMSD ignore >8 as requested
    tm_num = [v for v in tm_vals if v is not None]
    ptm_num = [v for v in ptm_vals if v is not None]
    iptm_num = [v for v in iptm_vals if v is not None]
    rmsd_num_all = [v for v in rmsd_vals if v is not None]
    rmsd_num = rmsd_num_all  # use ALL values for frequency denominator

    lddt_num = [v for v in lddt_vals if v is not None]

    excluded_count = len(rmsd_num_all) - len(rmsd_num)
    if excluded_count > 0:
        print(f"[INFO] Excluded {excluded_count} RMSD values > 8.0 from RMSD plot.")

    print(f"Records total: {len(records)}")
    print(f"TM numeric: {len(tm_num)} | pTM numeric: {len(ptm_num)} | iPTM numeric: {len(iptm_num)} | RMSD numeric (<=8): {len(rmsd_num)}\n")

    out_ptm = outdir / "cdf_ptm.png"
    out_iptm = outdir / "cdf_iptm.png"
    out_tm = outdir / "cdf_tm.png"
    out_rmsd = outdir / "cdf_rmsd.png"
    out_lddt = outdir / "cdf_lddt.png"
    out_inf_all = outdir / "cdf_inf_all.png"
    out_inf_wc = outdir / "cdf_inf_wc.png"

    # ----------------- NEW: category-aware plotting -----------------
    # parse nuc_type_count_per_pdb.txt (if present) into mapping
    nuc_type_path = input_path.parent / "../../MODNAP/MODNAP_features/nuc_type_count_per_pdb.txt"
    category_map = parse_nuc_type_file(nuc_type_path)
    passing_path = input_path.parent / "../../posebusters/passing_pdbs.txt"
    passing_pairs = parse_passing_pdbs(passing_path)

    # parse pymol_rmsd into mapping for category-specific RMSD extraction
    pymol_rmsd_map = parse_pymol_rmsd_map(pymol_path)


    categories = [
        "DNA",
        "RNA",
        "DNA-protein",
        "RNA-protein",
        "DNA-RNA",
        "DNA-RNA-protein",
    ]
   
    # ============================================================
    # INF_all groups
    # ============================================================

    inf_all_groups = {
        "All": [
            r["inf_all"]
            for r in inf_records
            if r["inf_all"] is not None
        ]
    }

    cat_inf_all = {}

    for rec in inf_records:

        cat = find_category_for_record(
            {"ccd": rec["ccd"], "pdb": rec["pdb"]},
            category_map
        )

        if cat is None:
            continue

        if rec["inf_all"] is None:
            continue

        cat_inf_all.setdefault(cat, []).append(rec["inf_all"])

    for cat in categories:
        inf_all_groups[cat] = cat_inf_all.get(cat, [])


# ============================================================
# INF_WC groups
# ============================================================

    inf_wc_groups = {
        "All": [
            r["inf_wc"]
            for r in inf_records
            if r["inf_wc"] is not None
        ]
    }

    cat_inf_wc = {}

    for rec in inf_records:

        cat = find_category_for_record(
            {"ccd": rec["ccd"], "pdb": rec["pdb"]},
            category_map
        )

        if cat is None:
            continue

        if rec["inf_wc"] is None:
            continue
    
        cat_inf_wc.setdefault(cat, []).append(rec["inf_wc"])

    for cat in categories:
        inf_wc_groups[cat] = cat_inf_wc.get(cat, [])

    # PTM
    ptm_groups = {"All": ptm_num[:]}
    # collect category-specific ptm values from records
    cat_ptm = collect_metric_values_by_category(records, "ptm", category_map, rmsd_map=None)
    for cat in categories:
        ptm_groups[cat] = cat_ptm.get(cat, [])
    # iPTM
    iptm_groups = {"All": iptm_num[:]}
    cat_iptm = collect_metric_values_by_category(records, "iptm", category_map, rmsd_map=None)
    for cat in categories:
        iptm_groups[cat] = cat_iptm.get(cat, [])
    # TM
    tm_groups = {"All": tm_num[:]}
    cat_tm = collect_metric_values_by_category(records, "tm", category_map, rmsd_map=None)
    for cat in categories:
        tm_groups[cat] = cat_tm.get(cat, [])
    # LDDT
    lddt_groups = {"All": lddt_num[:]}
    cat_lddt = collect_metric_values_by_category(records, "lddt", category_map, rmsd_map=None)
    for cat in categories:
        lddt_groups[cat] = cat_lddt.get(cat, [])
    # RMSD: for overall we still use rmsd_num (already filtered to <=8),
    # but for categories we must extract per-(ccd,pdb) values from pymol_rmsd_map
    # RMSD groups (All + per category)
    rmsd_groups = {"All": rmsd_num_all[:]}

    cat_rmsd = {cat: [] for cat in categories}

    for (ccd, pdb), rmsd_val in pymol_rmsd_map.items():

        cat = find_category_for_record({"ccd": ccd, "pdb": pdb}, category_map)
        if cat is None:
            continue

        if cat in cat_rmsd:
            cat_rmsd[cat].append(rmsd_val)

    for cat in categories:
        rmsd_groups[cat] = cat_rmsd.get(cat, [])

    # build PB-valid set
    pb_set = set(passing_pairs)
    
    # initialize groups FIRST
    success_rmsd_groups = {cat: [] for cat in categories}
    
    # assign RMSD using category_map (NO JSON)
    for (ccd, pdb) in pb_set:

        rmsd_val = None
        for k in [
            (ccd, pdb),
            (strip_ok(ccd), pdb),
            (ccd, strip_ok(pdb)),
            (strip_ok(ccd), strip_ok(pdb)),
        ]:
            if k in pymol_rmsd_map:
                rmsd_val = pymol_rmsd_map[k]
                break

        if rmsd_val is None:
            continue

        cat = find_category_for_record({"ccd": ccd, "pdb": pdb}, category_map)
        if cat is None:
            continue

        if cat in success_rmsd_groups:
            success_rmsd_groups[cat].append(rmsd_val)
    
    # NOW compute totals AFTER groups exist
    category_totals = {
        cat: len(rmsd_groups.get(cat, []))
        for cat in categories
    }

    pb_rates = {
        cat: (
            len(success_rmsd_groups[cat]) / category_totals[cat] * 100
            if category_totals[cat] else 0
        )
        for cat in categories
    }

    out_success = outdir / "rmsd_success_rate.png"

    r = make_and_save_success_rate_plot(
        rmsd_groups,          # <-- ADD this (all-category RMSDs for stacked bars)
        success_rmsd_groups,  # <-- was the first arg before
        pb_rates,
        out_success,
        cutoffs=(3.0, 2.5, 2.0)
    )

    if r:
        made.append(r)

    # Now plot multi-line CDFs (All + each category)
    
    r = make_and_save_cdf_multi(ptm_groups, "pTM", out_ptm, grid_step=0.001)
    if r: made.append(r)
    r = make_and_save_cdf_multi(iptm_groups, "iPTM", out_iptm, grid_step=0.001)
    if r: made.append(r)
    r = make_and_save_cdf_multi(tm_groups, "TM score", out_tm, grid_step=0.001)
    if r: made.append(r)
    r = make_and_save_cdf_multi(rmsd_groups, "RMSD", out_rmsd, grid_step=0.1)
    if r: made.append(r)
    r = make_and_save_cdf_multi(lddt_groups, "LDDT", out_lddt, grid_step=0.001)
    if r: made.append(r)
    
    r = make_and_save_cdf_multi(
    inf_all_groups,
    "INF_all",
    out_inf_all,
    grid_step=0.001
    )
    if r:
        made.append(r)

    r = make_and_save_cdf_multi(
        inf_wc_groups,
        "INF_WC",
        out_inf_wc,
        grid_step=0.001
    )
    if r:
        made.append(r)

    # ---------------------------------------------------------------------

    if made:
        print("Saved plots:")
        for pth in made:
            print("  " + str(pth))
    else:
        print("[WARN] No plots were created (no numeric data found).")

if __name__ == "__main__":
    main()
