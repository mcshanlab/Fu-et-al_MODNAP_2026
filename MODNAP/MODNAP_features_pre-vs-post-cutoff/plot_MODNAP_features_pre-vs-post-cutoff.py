#!/usr/bin/env python3
"""
Usage:

python3 plot_MODNAP_features_pre-vs-post-cutoff.py --base ../MODNAP --excel ../../list-of-modified-nucleotides-in-the-PDB/all_nucleotide_CCDs_PDB.xlsx --outdir .

Pre/post cutoff text files are expected at ./pre_cutoff.txt and ./post_cutoff.txt
(same directory from which the script is run).  Each file must have at least two
whitespace-separated columns; the second column (index 1) is taken as the PDB ID.
Lines starting with '#' and blank lines are ignored.
"""
from pathlib import Path
import argparse
import json
import re
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import math
import csv

import matplotlib
matplotlib.use("Agg")

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.font_manager as fm

arial_exists = any("Arial" in f for f in fm.findSystemFonts())
matplotlib.rcParams['font.family'] = 'Arial' if arial_exists else 'DejaVu Sans'

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.use14corefonts'] = False

import matplotlib.pyplot as plt

LETTER_SUB_RE = re.compile(r"[A-Za-z]+")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

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


def collect_string_candidates(obj):
    candidates = []
    if isinstance(obj, dict):
        for v in obj.values():
            candidates.extend(collect_string_candidates(v))
    elif isinstance(obj, list) or isinstance(obj, tuple):
        for item in obj:
            candidates.extend(collect_string_candidates(item))
    elif isinstance(obj, str):
        candidates.append(obj)
    return candidates


def load_cutoff_set(path: Path) -> set:
    """
    Read a two-column (or more) whitespace-delimited file.
    Column index 1 (second column) is treated as the PDB ID.
    Returns a set of upper-cased PDB IDs.
    """
    path = Path(path)
    pdbs = set()
    if not path.exists():
        print(f"[WARN] Cutoff file not found: {path}")
        return pdbs
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                pdbs.add(parts[1].upper())
    return pdbs


# ---------------------------------------------------------------------------
# Length extraction
# ---------------------------------------------------------------------------

def extract_best_lengths_from_json(json_path: Path):
    data = load_json(json_path)
    if data is None:
        return (None, None)

    nuc_len = 0
    prot_len = 0
    found_nuc = False
    found_prot = False

    if "sequences" not in data:
        return (None, None)

    for entry in data["sequences"]:
        if not isinstance(entry, dict):
            continue
        if "dna" in entry:
            seq = entry["dna"].get("sequence", "")
            nuc_len += len(seq)
            found_nuc = True
        elif "rna" in entry:
            seq = entry["rna"].get("sequence", "")
            nuc_len += len(seq)
            found_nuc = True
        elif "protein" in entry:
            seq = entry["protein"].get("sequence", "")
            prot_len += len(seq)
            found_prot = True

    return (
        nuc_len if found_nuc else None,
        prot_len if found_prot else None,
    )


# ---------------------------------------------------------------------------
# DNA / RNA / protein classifier
# ---------------------------------------------------------------------------

def classify_json_nuc_prot(json_path):
    data = load_json(json_path)
    if data is None:
        return None, False

    found_dna = False
    found_rna = False
    found_protein = False

    def walk(obj):
        nonlocal found_dna, found_rna, found_protein
        if isinstance(obj, dict):
            for k, v in obj.items():
                key = str(k).strip().lower()
                if key == "dna":
                    found_dna = True
                elif key == "rna":
                    found_rna = True
                elif key == "protein":
                    found_protein = True
                walk(v)
        elif isinstance(obj, list):
            for x in obj:
                walk(x)

    walk(data)

    if found_dna and not found_rna:
        nuc_type = "DNA"
    elif found_rna and not found_dna:
        nuc_type = "RNA"
    elif found_dna and found_rna:
        nuc_type = "DNA+RNA"
    else:
        nuc_type = None

    return nuc_type, found_protein


# ---------------------------------------------------------------------------
# Modification base presence
# ---------------------------------------------------------------------------

def compute_mod_base_presence(json_path: Path):
    presence = {b: 0 for b in ("A", "T", "G", "C", "U")}
    data = load_json(json_path)
    if data is None:
        return presence

    found_structured = False
    if isinstance(data, dict):
        seqs = []
        if "sequences" in data and isinstance(data["sequences"], (list, tuple)):
            seqs = data["sequences"]
        else:
            for v in data.values():
                if isinstance(v, (list, tuple)):
                    for item in v:
                        if isinstance(item, dict) and ("dna" in item or "rna" in item):
                            seqs.append(item)

        for seq_entry in seqs:
            for key in ("dna", "rna"):
                if key in seq_entry and isinstance(seq_entry[key], dict):
                    nuc_block = seq_entry[key]
                    seq_str = nuc_block.get("sequence") or nuc_block.get("seq") or ""
                    if not isinstance(seq_str, str):
                        seq_str = str(seq_str)
                    seq_up = seq_str.upper()
                    mods = (nuc_block.get("modifications")
                            or nuc_block.get("mods")
                            or nuc_block.get("mod"))
                    if isinstance(mods, (list, tuple)):
                        found_structured = True
                        seen_positions = set()
                        for m in mods:
                            if not isinstance(m, dict):
                                continue
                            pos = None
                            for k in ("basePosition", "position", "pos",
                                      "index", "residueNumber"):
                                if k in m:
                                    try:
                                        pos = int(m[k])
                                        break
                                    except Exception:
                                        pos = None
                            if pos is None:
                                for v in m.values():
                                    if isinstance(v, int):
                                        pos = v
                                        break
                            if pos is None:
                                continue
                            idx = pos - 1 if pos > 0 else pos
                            if 0 <= idx < len(seq_up) and idx not in seen_positions:
                                base = seq_up[idx]
                                if base in presence:
                                    presence[base] = 1
                                    seen_positions.add(idx)

    if found_structured:
        return presence

    # Fallback heuristic
    strings = collect_string_candidates(data)
    keywords = ("mod", "modified", "modification", "mutation", "variant",
                 "residue", "position", "pos", "chem_mod", "mod_base")
    seen = set()
    for s in strings:
        s_low = s.lower()
        if any(k in s_low for k in keywords):
            for m in LETTER_SUB_RE.finditer(s):
                sub = m.group(0).upper()
                if len(sub) == 1 and sub in presence:
                    seen.add(sub)
                else:
                    for ch in sub:
                        if ch in presence:
                            seen.add(ch)
            for base in presence:
                if base in s or base.lower() in s:
                    seen.add(base)
    for b in presence:
        if b in seen:
            presence[b] = 1
    return presence


# ---------------------------------------------------------------------------
# Histogram helpers
# ---------------------------------------------------------------------------

def np_hist(values, bin_edges):
    vals = [int(v) for v in values]
    counts = [0] * (len(bin_edges) - 1)
    for v in vals:
        if v >= bin_edges[-1]:
            idx = len(bin_edges) - 2
        else:
            idx = None
            for i in range(len(bin_edges) - 1):
                if bin_edges[i] <= v < bin_edges[i + 1]:
                    idx = i
                    break
            if idx is None:
                continue
        counts[idx] += 1
    return counts, bin_edges


def make_grouped_histogram(values_pre, values_post, bin_width, xlabel, outpath,
                            bin_min=0):
    """
    Side-by-side bars for pre (blue) and post (red).
    Saves both PNG and PDF.
    """
    all_vals = list(values_pre) + list(values_post)
    if not all_vals:
        print(f"[WARN] No values for {xlabel}; skipping histogram")
        return None

    max_v = max(all_vals)
    bin_edges = list(range(bin_min, max_v + bin_width, bin_width))
    if bin_edges[-1] < max_v:
        bin_edges.append(bin_edges[-1] + bin_width)

    counts_pre,  _ = np_hist(values_pre,  bin_edges)
    counts_post, _ = np_hist(values_post, bin_edges)

    labels = [f"{bin_edges[i]}-{bin_edges[i+1]}"
              for i in range(len(bin_edges) - 1)]

    x = np.arange(len(labels))
    width = 0.4

    plt.figure(figsize=(12, 5))
    bars_pre  = plt.bar(x - width / 2, counts_pre,  width=width,
                        label="pre",  color="#9FC5E8", edgecolor="black",
                        linewidth=0.6)
    bars_post = plt.bar(x + width / 2, counts_post, width=width,
                        label="post", color="#FFB3B3", edgecolor="black",
                        linewidth=0.6)

    plt.xticks(x, labels, rotation=90, ha='right')
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(f"{xlabel}: pre vs post")
    plt.legend()

    for i, (a, b) in enumerate(zip(counts_pre, counts_post)):
        if a:
            plt.text(i - width / 2, a, str(a),
                     ha='center', va='bottom', fontsize=8)
        if b:
            plt.text(i + width / 2, b, str(b),
                     ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outpath), dpi=200)
    try:
        plt.savefig(str(outpath.with_suffix('.pdf')), format="pdf",
                    bbox_inches="tight")
    except Exception:
        pass
    plt.close()
    return outpath


# ---------------------------------------------------------------------------
# Protein vs NA ratio scatter (pre / post separately)
# ---------------------------------------------------------------------------

def plot_prot_na_ratio(per_pdb, outdir: Path):
    from matplotlib.patches import Patch

    bin_edges  = [0, 5, 10, 20, 50, 100, float("inf")]
    bin_labels = ["0–5", "5–10", "10–20", "20–50", "50–100", ">100"]
    colors_list = ["#FDDFDF", "#B6CDBD", "#c6dbef",
                   "#fdae6b", "#fb6a4a", "#a63603"]

    def collect(split):
        xs, ys, ratios = [], [], []
        for r in per_pdb:
            if r.get("note") != split:
                continue
            nuc  = r.get("nuc_len")
            prot = r.get("prot_len")
            if nuc is None or prot is None:
                continue
            try:
                nuc_v  = int(nuc)
                prot_v = int(prot)
            except Exception:
                continue
            if nuc_v <= 0 or prot_v <= 0:
                continue
            xs.append(prot_v)
            ys.append(nuc_v)
            ratios.append(prot_v / nuc_v)
        return xs, ys, ratios

    def make_plot(xs, ys, ratios, title, outpath):
        if not xs:
            print(f"[INFO] No data for {title}; skipping")
            return

        point_colors = []
        for rr in ratios:
            for i in range(len(bin_edges) - 1):
                if bin_edges[i] <= rr < bin_edges[i + 1]:
                    point_colors.append(colors_list[i])
                    break

        plt.figure(figsize=(8, 7))
        plt.scatter(xs, ys, s=40, c=point_colors, alpha=0.85,
                    edgecolors="k", linewidths=0.25)
        plt.gca().set_rasterized(False)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Number of Protein Residues")
        plt.ylabel("Number of Nucleic Acid Nucleotides")
        plt.title(title)

        x_vals = np.logspace(
            np.log10(min(xs)) - 0.2, np.log10(max(xs)) + 0.2, 400
        )
        for r in [1, 5, 10, 20, 50, 100]:
            y_line = x_vals / r
            plt.plot(x_vals, y_line, "--", color="#777777", linewidth=0.8)
            ix   = len(x_vals) // 3
            xpos = x_vals[ix]
            ypos = y_line[ix]
            plt.text(xpos, ypos * 1.12, f"{r}:1", fontsize=8,
                     bbox=dict(boxstyle="round,pad=0.2",
                               facecolor="white", alpha=0.6))

        legend_patches = [
            Patch(facecolor=colors_list[i], label=bin_labels[i])
            for i in range(len(bin_labels))
        ]
        plt.legend(handles=legend_patches, title="Protein:NA Ratio",
                   loc="upper left", frameon=True)
        plt.grid(which="both", linestyle=":", linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        plt.savefig(str(outpath.with_suffix(".png")), dpi=200)
        try:
            plt.savefig(str(outpath.with_suffix(".pdf")),
                        format="pdf", bbox_inches="tight")
        except Exception:
            pass
        plt.close()
        print(f"[INFO] Saved {outpath}")

    xs, ys, ratios = collect("pre")
    make_plot(xs, ys, ratios,
              title="Protein vs Nucleic Acid size – PRE",
              outpath=outdir / "prot_na_ratio_pre")

    xs, ys, ratios = collect("post")
    make_plot(xs, ys, ratios,
              title="Protein vs Nucleic Acid size – POST",
              outpath=outdir / "prot_na_ratio_post")


# ---------------------------------------------------------------------------
# Excel modification lookup
# ---------------------------------------------------------------------------

def load_mod_excel_lookup(path: Path) -> dict:
    if not path.exists():
        print(f"[WARN] Excel lookup not found: {path}")
        return {}

    df = pd.read_excel(path)
    df.columns = [c.strip() for c in df.columns]

    def is_yes(val):
        if isinstance(val, bool):
            return val
        if isinstance(val, (int, float)):
            return bool(val)
        return str(val).strip().lower() == "yes"

    lookup = {}
    for _, row in df.iterrows():
        ccd = str(row["CCD_ID"]).strip().upper()
        lookup[ccd] = {
            "base":      is_yes(row.get("Base Modified",      False)),
            "sugar":     is_yes(row.get("Sugar Modified",     False)),
            "phosphate": is_yes(row.get("Phosphate Modified", False)),
        }

    print(f"[INFO] Loaded {len(lookup)} CCD entries from Excel")
    return lookup


def get_ccd_from_json_path(json_path: Path) -> str:
    stem = json_path.stem
    if "_" in stem:
        return stem.split("_")[0].strip().upper()
    return json_path.parent.parent.name.strip().upper()


def get_modification_types_from_json(json_path: Path):
    ccd = get_ccd_from_json_path(json_path)
    return [ccd] if ccd else []


def classify_pdb_from_excel_lookup(mod_codes, excel_lookup):
    if not mod_codes:
        return "no modification"

    has_base = has_sugar = has_phos = False
    for code in mod_codes:
        entry = excel_lookup.get(code.upper().strip())
        if not entry:
            continue
        has_base  |= entry.get("base",      False)
        has_sugar |= entry.get("sugar",     False)
        has_phos  |= entry.get("phosphate", False)

    if not (has_base or has_sugar or has_phos):
        return "no modification"
    if has_base and has_sugar and has_phos:
        return "base+sugar+phosphate"
    if has_base and has_sugar:
        return "base+sugar"
    if has_base and has_phos:
        return "base+phosphate"
    if has_sugar and has_phos:
        return "sugar+phosphate"
    if has_base:
        return "base"
    if has_sugar:
        return "sugar"
    return "phosphate"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    default_base = (Path.home() / "Desktop" / "ML_aptamer_project"
                    / "step12_benchmark" / "benchmark")
    parser = argparse.ArgumentParser(
        description="Extract nucleotide/protein lengths from JSONs under "
                    "minibenchmark and plot histograms (pre vs post split)."
    )
    parser.add_argument("--base", "-b", default=str(default_base),
                        help=f"Base minibenchmark dir (default: {default_base})")
    parser.add_argument("--outdir", "-o", default=None,
                        help="Output directory for PNG/CSV/text (default: base dir)")
    parser.add_argument("--excel", "-e", default=None,
                        help="Path to all_nucleotide_CCDs_PDB.xlsx")
    parser.add_argument("--pre",  default="./pre_cutoff.txt",
                        help="Path to pre-cutoff PDB list  (default: ./pre_cutoff.txt)")
    parser.add_argument("--post", default="./post_cutoff.txt",
                        help="Path to post-cutoff PDB list (default: ./post_cutoff.txt)")
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    if not base.exists() or not base.is_dir():
        print(f"[ERROR] Base directory not found: {base}")
        return

    outdir = (Path(args.outdir).expanduser().resolve()
              if args.outdir else base)

    # ---- Load pre/post sets ----
    pre_set  = load_cutoff_set(args.pre)
    post_set = load_cutoff_set(args.post)
    print(f"[INFO] pre_set:  {len(pre_set)} PDBs")
    print(f"[INFO] post_set: {len(post_set)} PDBs")

    # ---- Load Excel modification lookup ----
    if args.excel:
        excel_path = Path(args.excel).expanduser().resolve()
    else:
        excel_path = base / "all_nucleotide_CCDs_PDB.xlsx"
    if not excel_path.exists():
        print(f"[ERROR] Excel lookup file not found: {excel_path}")
        print(f"[ERROR] Pass the correct path with --excel /path/to/file.xlsx")
        return
    MOD_EXCEL_LOOKUP = load_mod_excel_lookup(excel_path)
    if not MOD_EXCEL_LOOKUP:
        print(f"[ERROR] Excel lookup is empty: {excel_path}")
        return
    print(f"[INFO] Excel lookup loaded: {len(MOD_EXCEL_LOOKUP)} CCD entries")

    # -----------------------------------------------------------------------
    # Walk CCD -> PDB directories and collect lengths + split tag
    # -----------------------------------------------------------------------
    nuc_lengths  = []
    prot_lengths = []
    per_pdb      = []

    for ccd_dir in sorted(base.iterdir()):
        if not ccd_dir.is_dir():
            continue
        for pdb_dir in sorted(ccd_dir.iterdir()):
            if not pdb_dir.is_dir():
                continue
            json_files = sorted([
                p for p in pdb_dir.iterdir()
                if p.is_file()
                and p.suffix.lower() == ".json"
                and not p.name.startswith(".")
            ])

            pdb_id = pdb_dir.name.upper()

            if pdb_id in pre_set:
                split = "pre"
            elif pdb_id in post_set:
                split = "post"
            else:
                split = "unknown"

            if not json_files:
                per_pdb.append({
                    "ccd": ccd_dir.name, "pdb": pdb_id,
                    "nuc_len": None, "prot_len": None,
                    "json": None, "note": split,
                })
                continue

            jpath = json_files[0]
            nuc_len, prot_len = extract_best_lengths_from_json(jpath)

            if prot_len is None and nuc_len is not None:
                print(f"[INFO] {ccd_dir.name}/{pdb_id}: "
                      "nucleotide only; excluded from protein plot")

            if nuc_len  is not None:
                nuc_lengths.append(nuc_len)
            if prot_len is not None:
                prot_lengths.append(prot_len)

            per_pdb.append({
                "ccd": ccd_dir.name, "pdb": pdb_id,
                "nuc_len": nuc_len, "prot_len": prot_len,
                "json": str(jpath), "note": split,
            })

    # -----------------------------------------------------------------------
    # Missing-data report
    # -----------------------------------------------------------------------
    missing_json     = []
    json_parse_failed = []
    no_extractable   = []
    missing_nuc      = []
    missing_prot     = []

    for r in per_pdb:
        ccd  = r.get("ccd")
        pdb  = r.get("pdb")
        jstr = r.get("json")
        nuc  = r.get("nuc_len")
        prot = r.get("prot_len")

        if not jstr:
            missing_json.append((ccd, pdb, "no json"))
            continue

        json_path = Path(jstr)
        if not json_path.exists():
            missing_json.append((ccd, pdb, str(json_path)))
            continue

        data = load_json(json_path)
        if data is None:
            json_parse_failed.append((ccd, pdb, str(json_path)))
            continue

        if nuc is None and prot is None:
            no_extractable.append((ccd, pdb, str(json_path)))
        elif nuc is None:
            missing_nuc.append((ccd, pdb, str(json_path)))
        elif prot is None:
            missing_prot.append((ccd, pdb, str(json_path)))

    report_path = outdir / "lengths_missing_report.txt"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    groups = [
        ("PDBs with NO JSON file found",                     missing_json),
        ("PDBs with JSON parse failure",                      json_parse_failed),
        ("JSON parsed but NO extractable nucleotide/protein", no_extractable),
        ("NO extractable nucleotide (protein present)",       missing_nuc),
        ("NO extractable protein (nucleotide present)",       missing_prot),
    ]
    with report_path.open("w", encoding="utf-8") as fh:
        fh.write(f"Lengths missing / problematic report\n")
        fh.write(f"Generated for base: {base}\n\n")
        for title, lst in groups:
            fh.write(f"{title}  (count = {len(lst)})\n")
            fh.write("-" * 72 + "\n")
            if not lst:
                fh.write("(none)\n\n")
                continue
            for ccd, pdb, note in lst:
                fh.write(f"{ccd}/{pdb}\t{note}\n")
            fh.write("\n")

    print("\n[REPORT] Lengths missing / problematic summary:")
    print(f"  Report saved to: {report_path}")
    for title, lst in groups:
        print(f"  {title}: {len(lst)}")

    # -----------------------------------------------------------------------
    # CSV + text summary
    # -----------------------------------------------------------------------
    csv_path = outdir / "lengths_per_pdb.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["CCD", "PDB", "nuc_len", "prot_len", "split", "json"])
        for r in per_pdb:
            writer.writerow([
                r["ccd"], r["pdb"],
                r["nuc_len"] or "", r["prot_len"] or "",
                r["note"] or "", r["json"] or "",
            ])

    txt_path = outdir / "lengths_summary.txt"
    with txt_path.open("w", encoding="utf-8") as fh:
        fh.write(f"Summary generated for base: {base}\n")
        fh.write(f"Total PDB folders scanned: {len(per_pdb)}\n")
        fh.write(f"PDBs with nucleotide length: {len(nuc_lengths)}\n")
        fh.write(f"PDBs with protein length:    {len(prot_lengths)}\n\n")
        fh.write("Nucleotide length distribution:\n")
        fh.write(str(Counter(nuc_lengths)) + "\n\n")
        fh.write("Protein length distribution:\n")
        fh.write(str(Counter(prot_lengths)) + "\n\n")
        fh.write("CSV saved to: " + str(csv_path) + "\n")

    print(f"[INFO] Scanned {len(per_pdb)} PDB folders.")
    print(f"[INFO] Found {len(nuc_lengths)} nucleotide lengths "
          f"and {len(prot_lengths)} protein lengths.")
    print(f"[INFO] CSV saved to: {csv_path}")
    print(f"[INFO] Text summary saved to: {txt_path}")

    # -----------------------------------------------------------------------
    # Histograms: pre vs post side-by-side
    # -----------------------------------------------------------------------
    nuc_pre  = [r["nuc_len"]  for r in per_pdb
                if r["nuc_len"]  is not None and r["note"] == "pre"]
    nuc_post = [r["nuc_len"]  for r in per_pdb
                if r["nuc_len"]  is not None and r["note"] == "post"]

    prot_pre  = [r["prot_len"] for r in per_pdb
                 if r["prot_len"] is not None and r["note"] == "pre"]
    prot_post = [r["prot_len"] for r in per_pdb
                 if r["prot_len"] is not None and r["note"] == "post"]

    nuc_png  = outdir / "nuc_length_hist.png"
    prot_png = outdir / "prot_length_hist.png"

    nuc_plot  = make_grouped_histogram(
        nuc_pre, nuc_post, bin_width=10,
        xlabel="Nucleotide length", outpath=nuc_png)
    prot_plot = make_grouped_histogram(
        prot_pre, prot_post, bin_width=50,
        xlabel="Protein length", outpath=prot_png)

    if nuc_plot:
        print(f"[INFO] Nucleotide histogram saved to: {nuc_plot}")
    if prot_plot:
        print(f"[INFO] Protein histogram saved to: {prot_plot}")

    # -----------------------------------------------------------------------
    # Protein vs NA ratio scatter (separate plots for pre / post)
    # -----------------------------------------------------------------------
    plot_prot_na_ratio(per_pdb, outdir)

    # -----------------------------------------------------------------------
    # DNA / RNA / protein category bar chart (pre vs post)
    # -----------------------------------------------------------------------
    categories  = ["RNA", "DNA", "RNA-protein", "DNA-protein",
                   "DNA+RNA", "DNA+RNA-protein"]
    counts_pre  = {cat: 0 for cat in categories}
    counts_post = {cat: 0 for cat in categories}
    per_pdb_results = []

    for r in per_pdb:
        json_path_str = r.get("json")
        if not json_path_str:
            per_pdb_results.append(
                (r["ccd"], r["pdb"], None, False, "no_json", ""))
            continue

        json_path = Path(json_path_str)
        if not json_path.exists():
            per_pdb_results.append(
                (r["ccd"], r["pdb"], None, False, "missing_json",
                 json_path_str))
            continue

        nuc_type, has_prot = classify_json_nuc_prot(json_path)
        if nuc_type is None:
            per_pdb_results.append(
                (r["ccd"], r["pdb"], None, bool(has_prot),
                 "unknown", str(json_path)))
            continue

        target = counts_pre if r["note"] == "pre" else counts_post

        if has_prot:
            if nuc_type == "RNA":
                target["RNA-protein"] += 1;    cat_label = "RNA-protein"
            elif nuc_type == "DNA":
                target["DNA-protein"] += 1;    cat_label = "DNA-protein"
            else:
                target["DNA+RNA-protein"] += 1; cat_label = "DNA+RNA-protein"
                print(f"[DNA+RNA-protein]  CCD: {r['ccd']}  |  PDB: {r['pdb']}")
        else:
            if nuc_type == "RNA":
                target["RNA"] += 1;    cat_label = "RNA"
            elif nuc_type == "DNA":
                target["DNA"] += 1;    cat_label = "DNA"
            else:
                target["DNA+RNA"] += 1; cat_label = "DNA+RNA"
                print(f"[DNA+RNA]          CCD: {r['ccd']}  |  PDB: {r['pdb']}")

        per_pdb_results.append(
            (r["ccd"], r["pdb"], nuc_type, bool(has_prot),
             cat_label, str(json_path)))

    # Write per-PDB classification report
    report_path = outdir / "nuc_type_count_per_pdb.txt"
    with report_path.open("w", encoding="utf-8") as fh:
        fh.write("CCD\tPDB\tNucleotide\tHas_Protein\tCategory\tJSON\n")
        for ccd, pdb, nuc_type, has_prot, cat_label, jpath in per_pdb_results:
            fh.write(f"{ccd}\t{pdb}\t{nuc_type or ''}\t"
                     f"{'Yes' if has_prot else 'No'}\t{cat_label}\t{jpath}\n")
        fh.write("\nSummary counts:\n")
        for cat in categories:
            fh.write(f"{cat}: pre={counts_pre[cat]} "
                     f"post={counts_post[cat]}\n")
    print(f"[INFO] Per-PDB nucleotide/protein report saved to: {report_path}")

    # Bar chart: pre vs post side-by-side
    x     = np.arange(len(categories))
    width = 0.4

    plt.figure(figsize=(7, 4))
    bars_pre = plt.bar(
        x - width / 2,
        [counts_pre[c]  for c in categories],
        width=width, label="pre",
        color="#9FC5E8", edgecolor="black", linewidth=0.6,
    )
    bars_post = plt.bar(
        x + width / 2,
        [counts_post[c] for c in categories],
        width=width, label="post",
        color="#FFB3B3", edgecolor="black", linewidth=0.6,
    )
    plt.xticks(x, categories, rotation=0)
    plt.xlabel("Category")
    plt.ylabel("Count")
    plt.title("Count of Nucleotide/Protein Categories (pre vs post)")
    plt.legend()

    for bar in list(bars_pre) + list(bars_post):
        h = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, h,
                 str(int(h)), ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    cat_png = outdir / "nuc_type_count.png"
    cat_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(cat_png), dpi=200)
    try:
        plt.savefig(str(cat_png.with_suffix('.pdf')), format="pdf",
                    bbox_inches="tight")
    except Exception:
        pass
    plt.close()
    print(f"[INFO] Nucleotide/protein category count plot saved to: {cat_png}")

    # -----------------------------------------------------------------------
    # Modification base presence – separate pie charts for pre / post
    # -----------------------------------------------------------------------
    bases         = ("A", "T", "G", "C", "U")
    presence_pre  = {b: 0 for b in bases}
    presence_post = {b: 0 for b in bases}

    for r in per_pdb:
        json_path_str = r.get("json")
        split         = r.get("note")
        if split not in ("pre", "post") or not json_path_str:
            continue
        json_path = Path(json_path_str)
        if not json_path.exists():
            continue
        presence = compute_mod_base_presence(json_path)
        target   = presence_pre if split == "pre" else presence_post
        for b in bases:
            target[b] += presence.get(b, 0)

    total_votes = sum(presence_pre.values()) + sum(presence_post.values())
    if total_votes == 0:
        print("[INFO] No structured modification positions found; "
              "skipping mod_base pie charts")
    else:
        pie_colors = ['#C3DEDD', '#F0E2C3', '#F6C7B3', '#C8B6FF', '#F6B8D0']
        pie_labels = list(bases)

        for tag, pres in (("pre", presence_pre), ("post", presence_post)):
            plt.figure(figsize=(6, 6))
            plt.pie(
                [pres[l] for l in pie_labels],
                labels=pie_labels,
                colors=pie_colors,
                autopct=lambda pct: f"{pct:.1f}%" if pct > 0 else "",
                startangle=90,
            )
            plt.gca().set_rasterized(False)
            plt.title(f"Distribution of Modified Bases – {tag.upper()}")
            plt.axis("equal")

            png_path = outdir / f"mod_base_pie_{tag}.png"
            png_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(png_path), dpi=200)
            try:
                plt.savefig(str(png_path.with_suffix('.pdf')),
                            format="pdf", bbox_inches="tight")
            except Exception:
                pass
            plt.close()
            print(f"[INFO] Modification base pie chart saved to: {png_path}")

    # -----------------------------------------------------------------------
    # Modification type classification via Excel lookup (pre vs post)
    # -----------------------------------------------------------------------
    mod_type_categories = [
        "no modification", "base", "sugar", "phosphate",
        "base+sugar", "base+phosphate", "sugar+phosphate",
        "base+sugar+phosphate",
    ]
    mod_type_counts_pre  = {cat: 0 for cat in mod_type_categories}
    mod_type_counts_post = {cat: 0 for cat in mod_type_categories}
    per_pdb_mod_type     = []

    print("\n[INFO] Starting modification type classification (Excel mapping)...")

    for r in per_pdb:
        ccd           = r.get("ccd")
        pdb           = r.get("pdb")
        json_path_str = r.get("json")
        split         = r.get("note", "unknown")

        if not json_path_str:
            cat = "no modification"
        else:
            json_path_obj = Path(json_path_str)
            if not json_path_obj.exists():
                cat = "no modification"
            else:
                mod_codes = get_modification_types_from_json(json_path_obj)
                cat = classify_pdb_from_excel_lookup(mod_codes, MOD_EXCEL_LOOKUP)

        if split == "pre":
            mod_type_counts_pre[cat]  += 1
        elif split == "post":
            mod_type_counts_post[cat] += 1

        per_pdb_mod_type.append((ccd, pdb, split, cat))

    total_mod = (sum(mod_type_counts_pre.values())
                 + sum(mod_type_counts_post.values()))
    assert total_mod == len(per_pdb), (
        f"[ERROR] total {total_mod} != per_pdb {len(per_pdb)}"
    )

    # Write report
    mod_type_report_path = outdir / "mod_type_per_pdb.txt"
    with mod_type_report_path.open("w", encoding="utf-8") as fh:
        fh.write("CCD\tPDB\tSplit\tModification_Category\n")
        for ccd, pdb, split, cat in per_pdb_mod_type:
            fh.write(f"{ccd}\t{pdb}\t{split}\t{cat}\n")
        fh.write("\nSummary counts:\n")
        for cat in mod_type_categories:
            fh.write(f"  {cat}: "
                     f"pre={mod_type_counts_pre[cat]} "
                     f"post={mod_type_counts_post[cat]}\n")
        fh.write(f"\nTotal PDBs: {total_mod}\n")
    print(f"[INFO] Per-PDB modification type report saved to: "
          f"{mod_type_report_path}")

    # Bar chart: pre vs post side-by-side
    bar_xlabels = [
        "no\nmod", "base", "sugar", "phos",
        "base+\nsugar", "base+\nphos", "sugar+\nphos", "all three",
    ]
    pre_vals  = [mod_type_counts_pre[c]  for c in mod_type_categories]
    post_vals = [mod_type_counts_post[c] for c in mod_type_categories]

    x     = np.arange(len(mod_type_categories))
    width = 0.4

    plt.figure(figsize=(11, 5))
    bars_pre = plt.bar(
        x - width / 2, pre_vals,  width=width, label="pre",
        color="#9FC5E8", edgecolor='black', linewidth=0.6,
    )
    bars_post = plt.bar(
        x + width / 2, post_vals, width=width, label="post",
        color="#FFB3B3", edgecolor='black', linewidth=0.6,
    )
    plt.gca().set_rasterized(False)
    plt.xticks(x, bar_xlabels, fontsize=9)
    plt.ylabel("Number of PDBs")
    plt.title("Modification Type Distribution (pre vs post)")
    plt.legend()

    max_val = max(pre_vals + post_vals) if (pre_vals + post_vals) else 1
    y_off   = max(1, max_val * 0.01)
    for i, (a, b) in enumerate(zip(pre_vals, post_vals)):
        if a:
            plt.text(i - width / 2, a + y_off, f"n={a}",
                     ha="center", va="bottom", fontsize=8)
        if b:
            plt.text(i + width / 2, b + y_off, f"n={b}",
                     ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    mod_type_bar_path = outdir / "mod_type_bar.png"
    plt.savefig(str(mod_type_bar_path), dpi=200)
    try:
        plt.savefig(str(mod_type_bar_path.with_suffix(".pdf")),
                    format="pdf", bbox_inches="tight")
        print(f"[INFO] Modification type bar PDF saved to: "
              f"{mod_type_bar_path.with_suffix('.pdf')}")
    except Exception as e:
        print(f"[WARN] Could not save PDF for mod_type_bar: {e}")
    plt.close()
    print(f"[INFO] Modification type bar chart saved to: {mod_type_bar_path}")


if __name__ == "__main__":
    main()
