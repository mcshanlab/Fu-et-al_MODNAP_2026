#!/usr/bin/env python3
"""
usage: python3 convert-CCD-to-images.py <input.csv|xlsx> <output.html>

python3 convert-CCD-to-images.py all_nucleotide_CCDs_PDB.xlsx list_of_modified_nucleotides_CCD_IDs_final.html

Generate an HTML document with CCD ligands using SVG images from RCSB.
Each CCD ID is a clickable link to its RCSB ligand page.
Modifications (Base, Sugar, Phosphate) are read directly from the input CSV/XLSX file.
"""

import os
import sys
import pandas as pd

# ---------------- Helpers ----------------

def load_dataframe(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xls", ".xlsx"]:
        df = pd.read_excel(path, dtype=str)
    elif ext == ".csv":
        df = pd.read_csv(path, dtype=str)
    else:
        raise ValueError(f"Unsupported file type: {ext}")
    return df.fillna('')

def get_svg_url(ccd_id):
    first_char = ccd_id[0].upper()
    return f"https://cdn.rcsb.org/images/ccd/labeled/{first_char}/{ccd_id.upper()}.svg"

def get_rcsb_link(ccd_id):
    return f"https://www.rcsb.org/ligand/{ccd_id.upper()}"

# ---------------- HTML Builder ----------------

def build_html(input_file, output_file="ligands.html"):
    df = load_dataframe(input_file)

    # Count modifications
    base_count = df['Base Modified'].str.strip().str.lower().eq('yes').sum()
    sugar_count = df['Sugar Modified'].str.strip().str.lower().eq('yes').sum()
    phosphate_count = df['Phosphate Modified'].str.strip().str.lower().eq('yes').sum()

    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>CCD Ligands Report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; }}
h1 {{ text-align: center; }}
h2 {{ margin-top: 40px; }}
img {{ width: 500px; height: auto; margin-top: 10px; }}  /* adjust size as needed */
.ligand {{ margin-bottom: 50px; }}
a {{ text-decoration: none; color: #1a0dab; }}
a:hover {{ text-decoration: underline; }}
.mod {{ font-weight: bold; }}
.summary {{ font-size: 1.1em; margin-bottom: 30px; }}
</style>
</head>
<body>
<h1>CCD Ligands Report</h1>
<p class="summary">Summary: Base modified = {base_count}, Sugar modified = {sugar_count}, Phosphate modified = {phosphate_count}
<br>
<p class="page">Click the CCD ID name to go to its page on the PDB</p>
"""

    for idx, row in df.iterrows():
        ccd_id = str(row.get("CCD_ID","")).strip()
        name = str(row.get("Name","")).strip()
        base_mod = str(row.get("Base Modified","")).strip()
        sugar_mod = str(row.get("Sugar Modified","")).strip()
        phosphate_mod = str(row.get("Phosphate Modified","")).strip()

        svg_url = get_svg_url(ccd_id)
        rcsb_url = get_rcsb_link(ccd_id)

        mod_text = f"Base: {base_mod}, Sugar: {sugar_mod}, Phosphate: {phosphate_mod}"

        html_content += f"""
<div class="ligand">
    <h2><a href="{rcsb_url}" target="_blank">{ccd_id}</a></h2>
    <p>{name}</p>
    <p class="mod">Modifications: {mod_text}</p>
    <img src="{svg_url}" alt="{ccd_id}" loading="eager" decoding="sync">
</div>
"""
        print(f"Added {ccd_id} - {name} | {mod_text}")

    html_content += """
</body>
</html>
"""

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"✅ HTML report generated: {output_file}")

# ---------------- Main ----------------

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python ccd_html_from_csv.py <input.csv|xlsx> <output.html>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    build_html(input_file, output_file)
