#!/usr/bin/env python3
"""
Usage:

python3 plot_MODNAP_features.py --base ../MODNAP --excel ../../list-of-modified-nucleotides-in-the-PDB/all_nucleotide_CCDs_PDB.xlsx --outdir .

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

# Already have these
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.font_manager as fm

# Automatically use Arial if available, otherwise fallback to DejaVu Sans
arial_exists = any("Arial" in f for f in fm.findSystemFonts())
matplotlib.rcParams['font.family'] = 'Arial' if arial_exists else 'DejaVu Sans'

# ADD THESE:
matplotlib.rcParams['text.usetex'] = False       # disable LaTeX rendering
matplotlib.rcParams['svg.fonttype'] = 'none'    # keep text as text in SVG
matplotlib.rcParams['pdf.use14corefonts'] = False

import matplotlib.pyplot as plt

# Regexes for detecting sequences


# tolerant regex to find letter-only substrings inside longer strings
LETTER_SUB_RE = re.compile(r"[A-Za-z]+")

def load_json(path: Path):
    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception as e:
        # fallback: try to read as text then parse
        try:
            txt = path.read_text(errors="ignore")
            return json.loads(txt)
        except Exception:
            return None

def collect_string_candidates(obj):
    """
    Recursively walk JSON-like obj and return a list of all string values found.
    """
    candidates = []
    if isinstance(obj, dict):
        for v in obj.values():
            candidates.extend(collect_string_candidates(v))
    elif isinstance(obj, list) or isinstance(obj, tuple):
        for item in obj:
            candidates.extend(collect_string_candidates(item))
    elif isinstance(obj, str):
        candidates.append(obj)
    else:
        # ignore numbers, booleans, None
        pass
    return candidates

def score_as_nucleotide(s: str):
    """
    Return fraction of characters in s that are nucleotide characters.
    """
    if not s:
        return 0.0
    letters = [ch for ch in s if ch.isalpha()]
    if not letters:
        return 0.0
    n_nuc = sum(1 for ch in letters if ch in NUC_CHARS)
    return n_nuc / len(letters)

def score_as_protein(s: str):
    """
    Return fraction of characters in s that are amino-acid letters.
    """
    if not s:
        return 0.0
    letters = [ch for ch in s if ch.isalpha()]
    if not letters:
        return 0.0
    n_aa = sum(1 for ch in letters if ch in AA_CHARS)
    return n_aa / len(letters)

def extract_best_lengths_from_json(json_path: Path):
    """
    Extract nucleotide and protein lengths using explicit AF3 JSON structure.
    No sequence-letter guessing.
    """
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
        prot_len if found_prot else None
    )

# ----------------- NEW: classify JSON for DNA/RNA and protein presence -----------------
def classify_json_nuc_prot(json_path):
    """
    Classify nucleotide type and protein presence using ONLY explicit JSON keys.
    No sequence-content (U/T) heuristics are used.
    """
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

    
# ------------------------------------------------------------------------------

# ----------------- NEW: compute modification base presence per PDB -----------------
def compute_mod_base_presence(json_path: Path):
    """
    Return a dict {'A':0/1,'T':0/1,'G':0/1,'C':0/1,'U':0/1} indicating whether
    that base is modified anywhere in this JSON (1 means modified in this PDB).
    Rules:
      - Prefer structured data under data['sequences'] -> items with 'dna' or 'rna'
        that have a 'sequence' string and a 'modifications' list with a position key.
      - For each modification, map basePosition -> base in sequence (1-based assumed).
      - Count each base position only once per PDB (seen_positions prevents duplicates).
      - If structured info not found, fall back to the older heuristic token-scan,
        but still count presence once per base.
    """
    presence = {b: 0 for b in ("A", "T", "G", "C", "U")}
    data = load_json(json_path)
    if data is None:
        return presence

    found_structured = False
    # Try structured route first
    if isinstance(data, dict):
        seqs = []
        if "sequences" in data and isinstance(data["sequences"], (list, tuple)):
            seqs = data["sequences"]
        else:
            # discover containers with dna/rna entries
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
                    mods = nuc_block.get("modifications") or nuc_block.get("mods") or nuc_block.get("mod")
                    if isinstance(mods, (list, tuple)):
                        found_structured = True
                        # record which indices seen for this PDB to avoid double counting same base position
                        seen_positions = set()
                        for m in mods:
                            if not isinstance(m, dict):
                                continue
                            # find a numeric pos in common keys
                            pos = None
                            for k in ("basePosition", "position", "pos", "index", "residueNumber"):
                                if k in m:
                                    try:
                                        pos = int(m[k])
                                        break
                                    except Exception:
                                        pos = None
                            if pos is None:
                                # fallback: try any integer-valued value in m
                                for v in m.values():
                                    if isinstance(v, int):
                                        pos = v
                                        break
                            if pos is None:
                                continue
                            # convert to 0-based index (assume 1-based if pos>0)
                            idx = pos - 1 if pos > 0 else pos
                            if 0 <= idx < len(seq_up) and idx not in seen_positions:
                                base = seq_up[idx]
                                if base in presence:
                                    presence[base] = 1
                                    seen_positions.add(idx)
    # If structured info found, return presence
    if found_structured:
        return presence

    # Fallback: heuristic scanning, but only mark presence (1) per base per JSON
    strings = collect_string_candidates(data)
    keywords = ("mod", "modified", "modification", "mutation", "variant", "residue", "position", "pos", "chem_mod", "mod_base")
    seen = set()
    for s in strings:
        s_low = s.lower()
        if any(k in s_low for k in keywords):
            for m in LETTER_SUB_RE.finditer(s):
                sub = m.group(0).upper()
                if len(sub) == 1 and sub in presence:
                    seen.add(sub)
                else:
                    # tokens like m6A or 5mC: check letters inside
                    for ch in sub:
                        if ch in presence:
                            seen.add(ch)
            # also check raw chars in the string
            for base in presence:
                if base in s or base.lower() in s:
                    seen.add(base)
    for b in presence:
        if b in seen:
            presence[b] = 1
    return presence
# ------------------------------------------------------------------------------

# --- inside make_histogram_and_save ---

def make_histogram_and_save(values, bin_width, xlabel, outpath, bin_min=0):
    counts, edges = np_hist(values, range(bin_min, max(values)+bin_width, bin_width))
    labels = [f"{int(edges[i])}-{int(edges[i+1])}" for i in range(len(edges)-1)]

    plt.figure(figsize=(10,5))
    x = range(len(counts))        # positions for bars
    cat_values = counts           # bar heights
    bars = plt.bar(
        x,
        cat_values,
        align='center',
        color=["#FFC9C8", "#FFD6A5", "#CDE5D9", "#D0F0FF"]
    )
    plt.gca().set_rasterized(False)
    plt.xticks(range(len(counts)), labels, rotation=90, ha='right')
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(f"Count of {xlabel}")

    # Add count labels on top of each bar
    try:
        max_count = max(counts) if counts else 0
        y_offset = max(1, max_count * 0.01)
        for i, c in enumerate(counts):
            plt.text(i, c + y_offset, str(c), ha='center', va='bottom', fontsize=8)
    except Exception:
        pass

    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outpath), dpi=200)

    # also save PDF version
    try:
        pdf_path = outpath.with_suffix('.pdf')
        plt.savefig(str(pdf_path), format="pdf", dpi=300, bbox_inches="tight")
    except Exception:
        pass
    plt.close()
    return outpath
    

def np_hist(values, bin_edges):
    """
    Lightweight histogram that returns counts and edges similar to numpy.histogram,
    but ensures we include last bin even if value == edges[-1]
    """
    # convert to ints
    vals = [int(v) for v in values]
    counts = [0] * (len(bin_edges) - 1)
    for v in vals:
        # find bin index
        # if v == last edge, place into last bin
        if v >= bin_edges[-1]:
            idx = len(bin_edges) - 2
        else:
            # find i s.t. bin_edges[i] <= v < bin_edges[i+1]
            # naive linear search (bins are small)
            idx = None
            for i in range(len(bin_edges)-1):
                if bin_edges[i] <= v < bin_edges[i+1]:
                    idx = i
                    break
            if idx is None:
                # if somehow not found, skip
                continue
        counts[idx] += 1
    return counts, bin_edges

def plot_prot_na_ratio(per_pdb, outdir: Path):
    """
    Scatter plot of Number of Nucleotide (y) vs Number of Protein residues (x)
    for PDBs that have both lengths. Points colored by Protein:NA ratio bins.
    Saves to outdir / "prot_na_ratio.png".
    """
    import numpy as np
    from matplotlib.patches import Patch

    # collect pairs with both lengths present and > 0
    xs = []
    ys = []
    ratios = []

    for r in per_pdb:
        try:
            nuc = r.get("nuc_len")
            prot = r.get("prot_len")
        except Exception:
            continue
        if nuc is None or prot is None:
            continue
        try:
            nuc_v = int(nuc)
            prot_v = int(prot)
        except Exception:
            continue
        if nuc_v <= 0 or prot_v <= 0:
            continue
        xs.append(prot_v)
        ys.append(nuc_v)
        ratios.append(prot_v / nuc_v)

    if not xs:
        print("[INFO] No PDBs with both nucleotide and protein lengths found; skipping prot_na_ratio.png")
        return None

    # ratio bins and colors (Protein:NA ratio)
    bin_edges = [0, 5, 10, 20, 50, 100, float("inf")]
    bin_labels = ["Ratio 0-5", "Ratio 5-10", "Ratio 10-20", "Ratio 20-50", "Ratio 50-100", "Ratio >100"]
    colors_list = ["#FDDFDF", "#B6CDBD", "#c6dbef", "#fdae6b", "#fb6a4a", "#a63603"]

    # assign colors per point
    point_colors = []
    for rr in ratios:
        for i in range(len(bin_edges)-1):
            if bin_edges[i] <= rr < bin_edges[i+1]:
                point_colors.append(colors_list[i])
                break

    # prepare figure
    plt.figure(figsize=(8,7))
    plt.scatter(
        xs,
        ys,
        s=40,
        c=point_colors,
        alpha=0.85,
        edgecolors="k",
        linewidths=0.25,
    )
    plt.gca().set_rasterized(False)   # <-- prevent rasterization for vector output

    # log scales
    plt.xscale("log")
    plt.yscale("log")

    # axis labels and title
    plt.xlabel("Number of Protein Residues")
    plt.ylabel("Number of Nucleic Acid Nucleotides")
    plt.title("Protein vs Nucleic Acid size (colored by Protein:NA ratio)")

    # diagonal ratio lines (Protein:NA = r  ->  y = x / r)
    x_min = min(xs)
    x_max = max(xs)
    x_vals = np.logspace(np.log10(x_min) - 0.2, np.log10(x_max) + 0.2, 400)
    ratio_line_vals = [1, 5, 10, 20, 50, 100]
    for r in ratio_line_vals:
        y_line = x_vals / r
        plt.plot(x_vals, y_line, linestyle="--", color="#777777", linewidth=0.8)
        ix = len(x_vals) // 3
        xpos = x_vals[ix]
        ypos = y_line[ix]
        plt.text(xpos, ypos * 1.12, f"{r}:1", fontsize=8, bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.6))

    # legend for ratio color bins
    legend_patches = [Patch(facecolor=colors_list[i], label=bin_labels[i]) for i in range(len(bin_labels))]
    plt.legend(handles=legend_patches, title="Protein:NA Ratio", loc="upper left", frameon=True)

    plt.grid(which="both", linestyle=":", linewidth=0.5, alpha=0.7)

    outpath = outdir / "prot_na_ratio.png"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(str(outpath), dpi=200)
    # also save PDF
    try:
        pdf_out = outpath.with_suffix('.pdf')
        plt.savefig(str(pdf_out))
    except Exception:
        pass
    plt.close()
    print(f"[INFO] Protein vs NA ratio plot saved to: {outpath}")
    return outpath


def load_mod_excel_lookup(path: Path) -> dict:
    if not path.exists():
        print(f"[WARN] Excel lookup not found: {path}")
        return {}

    df = pd.read_excel(path)
    df.columns = [c.strip() for c in df.columns]

    def is_yes(val):
        """Treat "Yes"/"yes"/"YES"/True/1 as True; everything else False."""
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
    """
    Extract the CCD code from the JSON filename.
    Convention: <CCD>_<PDBID>.json  e.g. 0C_1R3O.json -> CCD = "0C"
    Falls back to the parent directory name if no underscore is found.
    """
    stem = json_path.stem          # e.g. "0C_1R3O"
    if "_" in stem:
        return stem.split("_")[0].strip().upper()
    # fallback: try the CCD dir name (grandparent)
    return json_path.parent.parent.name.strip().upper()


def get_modification_types_from_json(json_path: Path):
    """
    Return the CCD code for this JSON as a one-element list, derived purely
    from the filename (e.g. 0C_1R3O.json -> ["0C"]).
    No JSON parsing needed for the modification lookup.
    """
    ccd = get_ccd_from_json_path(json_path)
    return [ccd] if ccd else []


def classify_pdb_from_excel_lookup(mod_codes, excel_lookup):
    """
    Map CCD codes -> (base/sugar/phosphate) using Excel lookup.
    """
    if not mod_codes:
        return "no modification"

    has_base = has_sugar = has_phos = False

    for code in mod_codes:
        entry = excel_lookup.get(code.upper().strip())
        if not entry:
            continue

        has_base |= entry.get("base", False)
        has_sugar |= entry.get("sugar", False)
        has_phos |= entry.get("phosphate", False)

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


def main():
    default_base = Path.home() / "Desktop" / "ML_aptamer_project" / "step12_benchmark" / "benchmark"
    parser = argparse.ArgumentParser(description="Extract nucleotide/protein lengths from JSONs under minibenchmark and plot histograms.")
    parser.add_argument("--base", "-b", default=str(default_base), help=f"Base minibenchmark dir (default: {default_base})")
    parser.add_argument("--outdir", "-o", default=None, help="Output directory for PNG/CSV/text (default: base dir)")
    parser.add_argument("--excel", "-e", default=None, help="Path to all_nucleotide_CCDs_PDB.xlsx (default: <base>/all_nucleotide_CCDs_PDB.xlsx)")
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    if not base.exists() or not base.is_dir():
        print(f"[ERROR] Base directory not found: {base}")
        return

    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else base

    # Load Excel modification lookup (scoped inside main, passed explicitly where needed)
    if args.excel:
        excel_path = Path(args.excel).expanduser().resolve()
    else:
        excel_path = base / "all_nucleotide_CCDs_PDB.xlsx"
    if not excel_path.exists():
        print(f"[ERROR] Excel lookup file not found: {excel_path}")
        print(f"[ERROR] Pass the correct path with --excel /path/to/all_nucleotide_CCDs_PDB.xlsx")
        return
    MOD_EXCEL_LOOKUP = load_mod_excel_lookup(excel_path)
    if not MOD_EXCEL_LOOKUP:
        print(f"[ERROR] Excel lookup loaded but is empty — check the file: {excel_path}")
        return
    print(f"[INFO] Excel lookup loaded: {len(MOD_EXCEL_LOOKUP)} CCD entries from {excel_path}")

    # Walk CCD -> PDB -> find first .json in each PDB folder
    nuc_lengths = []
    prot_lengths = []
    per_pdb = []  # list of dicts

    for ccd_dir in sorted(base.iterdir()):
        if not ccd_dir.is_dir():
            continue
        for pdb_dir in sorted(ccd_dir.iterdir()):
            if not pdb_dir.is_dir():
                continue
            # find JSON files in this PDB folder
            json_files = sorted([p for p in pdb_dir.iterdir() if p.is_file() and p.suffix.lower() == ".json" and not p.name.startswith(".")])
            if not json_files:
                per_pdb.append({
                    "ccd": ccd_dir.name,
                    "pdb": pdb_dir.name,
                    "nuc_len": None,
                    "prot_len": None,
                    "json": None,
                    "note": "no json"
                })
                continue
            # choose first json (if there are multiple, you can change strategy)
            jpath = json_files[0]
            nuc_len, prot_len = extract_best_lengths_from_json(jpath)

            # ADD: log nucleotide-only PDBs
            if prot_len is None and nuc_len is not None:
                print(f"[INFO] {ccd_dir.name}/{pdb_dir.name}: nucleotide only; excluded from protein plot")

            if nuc_len is not None:
                nuc_lengths.append(nuc_len)
            if prot_len is not None:
                prot_lengths.append(prot_len)

            per_pdb.append({
                "ccd": ccd_dir.name,
                "pdb": pdb_dir.name,
                "nuc_len": nuc_len,
                "prot_len": prot_len,
                "json": str(jpath),
                "note": None
            })

    # ----------------- report which PDBs did not contribute to the nuc/prot histograms -------------
    missing_json = []
    json_parse_failed = []
    no_extractable = []
    missing_nuc = []
    missing_prot = []

    for r in per_pdb:
        ccd = r.get("ccd")
        pdb = r.get("pdb")
        jstr = r.get("json")
        nuc = r.get("nuc_len")
        prot = r.get("prot_len")

        if not jstr:
            missing_json.append((ccd, pdb, "no json"))
            continue

        json_path = Path(jstr)
        if not json_path.exists():
            alt = Path(jstr) if Path(jstr).is_absolute() else (outdir / Path(jstr))
            if alt.exists():
                json_path = alt

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
    with report_path.open("w", encoding="utf-8") as fh:
        fh.write(f"Lengths missing / problematic report\n")
        fh.write(f"Generated for base: {base}\n\n")

        groups = [
            ("PDBs with NO JSON file found", missing_json),
            ("PDBs with JSON parse failure", json_parse_failed),
            ("JSON parsed but NO extractable nucleotide/protein", no_extractable),
            ("NO extractable nucleotide (protein present)", missing_nuc),
            ("NO extractable protein (nucleotide present)", missing_prot),
        ]

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

    # Save CSV
    csv_path = outdir / "lengths_per_pdb.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["CCD", "PDB", "nuc_len", "prot_len", "json", "note"])
        for r in per_pdb:
            writer.writerow([r["ccd"], r["pdb"], r["nuc_len"] or "", r["prot_len"] or "", r["json"] or "", r["note"] or ""])

    # Save text summary
    txt_path = outdir / "lengths_summary.txt"
    with txt_path.open("w", encoding="utf-8") as fh:
        fh.write(f"Summary generated for base: {base}\n")
        fh.write(f"Total PDB folders scanned: {len(per_pdb)}\n")
        fh.write(f"PDBs with nucleotide length: {len(nuc_lengths)}\n")
        fh.write(f"PDBs with protein length: {len(prot_lengths)}\n\n")
        fh.write("Nucleotide length distribution (counts) sample:\n")
        fh.write(str(Counter(nuc_lengths)) + "\n\n")
        fh.write("Protein length distribution (counts) sample:\n")
        fh.write(str(Counter(prot_lengths)) + "\n\n")
        fh.write("CSV saved to: " + str(csv_path) + "\n")

    print(f"[INFO] Scanned {len(per_pdb)} PDB folders.")
    print(f"[INFO] Found {len(nuc_lengths)} nucleotide lengths and {len(prot_lengths)} protein lengths.")
    print(f"[INFO] CSV saved to: {csv_path}")
    print(f"[INFO] Text summary saved to: {txt_path}")

    # Make histograms: nuc bins width=10, prot bins width=50
    nuc_png = outdir / "nuc_length_hist.png"
    prot_png = outdir / "prot_length_hist.png"

    nuc_plot = make_histogram_and_save(nuc_lengths, bin_width=10, xlabel="Nucleotide length", outpath=nuc_png, bin_min=0)
    prot_plot = make_histogram_and_save(prot_lengths, bin_width=50, xlabel="Protein length", outpath=prot_png, bin_min=0)

    if nuc_plot:
        print(f"[INFO] Nucleotide histogram saved to: {nuc_plot}")
    if prot_plot:
        print(f"[INFO] Protein histogram saved to: {prot_plot}")

    # generate protein vs nucleotide ratio plot (only PDBs that have both lengths)
    plot_prot_na_ratio(per_pdb, outdir)

    # ----------------- tally DNA/RNA and protein combinations and plot -------------
    categories = ["RNA", "DNA", "RNA-protein", "DNA-protein", "DNA+RNA", "DNA+RNA-protein"]
    counts = dict((cat, 0) for cat in categories)

    per_pdb_results = []  # tuples: (ccd, pdb, nuc_type_or_None, has_protein_bool, category_label, json_path_str)

    for r in per_pdb:
        json_path_str = r.get("json")
        if not json_path_str:
            per_pdb_results.append((r["ccd"], r["pdb"], None, False, "no_json", ""))
            continue

        json_path = Path(json_path_str)
        if not json_path.exists():
            json_path = (outdir if json_path_str.startswith("/") == False else Path(json_path_str))
            if not json_path.exists():
                per_pdb_results.append((r["ccd"], r["pdb"], None, False, "missing_json", json_path_str))
                continue

        nuc_type, has_prot = classify_json_nuc_prot(json_path)
        if nuc_type is None:
            per_pdb_results.append((r["ccd"], r["pdb"], None, bool(has_prot), "unknown", str(json_path)))
            continue

        if has_prot:
            if nuc_type == "RNA":
                counts["RNA-protein"] += 1
                cat_label = "RNA-protein"
            elif nuc_type == "DNA":
                counts["DNA-protein"] += 1
                cat_label = "DNA-protein"
            else:  # DNA+RNA
                counts["DNA+RNA-protein"] += 1
                cat_label = "DNA+RNA-protein"
                print(f"[DNA+RNA-protein]  CCD: {r['ccd']}  |  PDB: {r['pdb']}")
        else:
            if nuc_type == "RNA":
                counts["RNA"] += 1
                cat_label = "RNA"
            elif nuc_type == "DNA":
                counts["DNA"] += 1
                cat_label = "DNA"
            else:  # DNA+RNA
                counts["DNA+RNA"] += 1
                cat_label = "DNA+RNA"
                print(f"[DNA+RNA]          CCD: {r['ccd']}  |  PDB: {r['pdb']}")

        per_pdb_results.append((r["ccd"], r["pdb"], nuc_type, bool(has_prot), cat_label, str(json_path)))

    # Write per-PDB classification report
    report_path = outdir / "nuc_type_count_per_pdb.txt"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as fh:
        fh.write("CCD\tPDB\tNucleotide\tHas_Protein\tCategory\tJSON\n")
        for ccd, pdb, nuc_type, has_prot, cat_label, jpath in per_pdb_results:
            fh.write(f"{ccd}\t{pdb}\t{nuc_type or ''}\t{'Yes' if has_prot else 'No'}\t{cat_label}\t{jpath}\n")
        fh.write("\nSummary counts:\n")
        for cat in categories:
            fh.write(f"{cat}: {counts.get(cat, 0)}\n")

    print(f"[INFO] Per-PDB nucleotide/protein report saved to: {report_path}")

    # Categorical bar plot
    cat_labels = categories
    cat_values = [counts[c] for c in cat_labels]

    plt.figure(figsize=(6,4))
    x = range(len(cat_labels))
    bars = plt.bar(
        x,
        cat_values,
        align='center',
        color=["#FFC9C8", "#FFD6A5", "#CDE5D9", "#D0F0FF", "#E8C5F0", "#B5D5C5"]
    )
    plt.xticks(x, cat_labels, rotation=0)
    plt.xlabel("Category")
    plt.ylabel("Count")
    plt.title("Count of Nucleotide/Protein Categories")

    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            str(int(height)),
            ha='center',
            va='bottom',
            fontsize=9
        )

    plt.tight_layout()
    cat_png = outdir / "nuc_type_count.png"
    cat_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(cat_png), dpi=200)
    try:
        cat_pdf = cat_png.with_suffix('.pdf')
        plt.savefig(str(cat_pdf))
    except Exception:
        pass
    plt.close()
    print(f"[INFO] Nucleotide/protein category count plot saved to: {cat_png}")

    # ----------------- modification base presence pie chart -------------
    total_presence = {b: 0 for b in ("A", "T", "G", "C", "U")}
    for r in per_pdb:
        json_path_str = r.get("json")
        if not json_path_str:
            continue
        json_path = Path(json_path_str)
        if not json_path.exists():
            json_path = (outdir if json_path_str.startswith("/") == False else Path(json_path_str))
            if not json_path.exists():
                continue
        presence = compute_mod_base_presence(json_path)
        for b in total_presence:
            total_presence[b] += presence.get(b, 0)

    total_votes = sum(total_presence.values())
    mod_pie_path = outdir / "mod_base_pie.png"
    if total_votes == 0:
        print("[INFO] No structured modification positions found in JSONs; skipping mod_base_pie.png")
    else:
        labels = ["A", "T", "G", "C", "U"]
        sizes = [total_presence[l] for l in labels]
        plt.figure(figsize=(6,6))
        colors = ['#C3DEDD', '#F0E2C3', '#F6C7B3', '#C8B6FF', '#F6B8D0']

        plt.pie(
            sizes,
            labels=labels,
            colors=colors,
            autopct=lambda pct: f"{pct:.1f}%" if pct > 0 else "",
            startangle=90
        )
        plt.gca().set_rasterized(False)
        plt.title("Distribution of Modified Bases (A,T,G,C,U) — per-PDB presence")
        plt.axis("equal")
        mod_pie_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(mod_pie_path), dpi=200)
        try:
            mod_pdf = mod_pie_path.with_suffix('.pdf')
            plt.savefig(str(mod_pdf))
        except Exception:
            pass
        plt.close()
        print(f"[INFO] Modification base pie chart saved to: {mod_pie_path}")

    # ================= MODIFICATION TYPE CLASSIFICATION =================
    mod_type_categories = [
        "no modification", "base", "sugar", "phosphate",
        "base+sugar", "base+phosphate", "sugar+phosphate",
        "base+sugar+phosphate",
    ]
    mod_type_counts = {cat: 0 for cat in mod_type_categories}
    per_pdb_mod_type = []

    print("\n[INFO] Starting modification type classification (Excel mapping)...")

    for r in per_pdb:
        ccd = r.get("ccd")
        pdb = r.get("pdb")
        json_path_str = r.get("json")

        if not json_path_str:
            cat = "no modification"
        else:
            json_path_obj = Path(json_path_str)
            if not json_path_obj.exists():
                cat = "no modification"
            else:
                mod_codes = get_modification_types_from_json(json_path_obj)
                cat = classify_pdb_from_excel_lookup(mod_codes, MOD_EXCEL_LOOKUP)

        mod_type_counts[cat] += 1
        per_pdb_mod_type.append((ccd, pdb, cat))

    assert sum(mod_type_counts.values()) == len(per_pdb)

    # Write report
    mod_type_report_path = outdir / "mod_type_per_pdb.txt"
    with mod_type_report_path.open("w", encoding="utf-8") as fh:
        fh.write("CCD\tPDB\tModification_Category\n")
        for ccd, pdb, cat in per_pdb_mod_type:
            fh.write(f"{ccd}\t{pdb}\t{cat}\n")
        fh.write("\nSummary counts:\n")
        for cat in mod_type_categories:
            fh.write(f"{cat}: {mod_type_counts[cat]}\n")

    print(f"[INFO] Per-PDB modification type report saved to: {mod_type_report_path}")

    # Bar plot
    bar_order = mod_type_categories
    bar_values = [mod_type_counts[c] for c in bar_order]

    plt.figure(figsize=(11, 5))
    x = range(len(bar_order))

    plt.bar(x, bar_values, color=[
        "#D3D3D3", "#FFB3B3", "#FFE599", "#9FC5E8",
        "#FFCBA4", "#C9B1D9", "#B6D7A8", "#EA9999",
    ])

    plt.xticks(x, [
        "no\nmod", "base", "sugar", "phos",
        "base+\nsugar", "base+\nphos", "sugar+\nphos",
        "all three"
    ], fontsize=9)

    plt.ylabel("Number of PDBs")
    plt.title("Modification Type Distribution")

    # Add n = x labels on top of each bar
    max_val = max(bar_values) if bar_values else 1
    y_offset = max(1, max_val * 0.01)
    for i, v in enumerate(bar_values):
            plt.text(i, v + y_offset, f"n={v}", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    mod_type_bar_path = outdir / "mod_type_bar.png"
    plt.savefig(mod_type_bar_path, dpi=200)
    try:
        mod_type_bar_pdf = mod_type_bar_path.with_suffix('.pdf')
        plt.savefig(str(mod_type_bar_pdf), format="pdf", bbox_inches="tight")
        print(f"[INFO] Modification type bar chart PDF saved to: {mod_type_bar_pdf}")
    except Exception as e:
        print(f"[WARN] Could not save PDF for mod_type_bar: {e}")
    plt.close()

    print(f"[INFO] Modification type bar chart saved to: {mod_type_bar_path}")
    # ================= END MODIFICATION TYPE BLOCK =================


if __name__ == "__main__":
    main()
