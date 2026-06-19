import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import Descriptors

# -----------------------------
# 1. Load Excel (first sheet)
# -----------------------------
file = "all_nucleotide_CCDs_PDB.xlsx"
df = pd.read_excel(file, sheet_name=0)

# Expected columns:
# CCD_ID | Name | SMILES
required_cols = ["CCD_ID", "Name", "SMILES"]

for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing required column: {col}")

# -----------------------------
# 2. Define canonical unmodified nucleotides
# -----------------------------
unmodified = set([
    "ADENOSINE",
    "ADENOSINE-5'-MONOPHOSPHATE",
    "ADENOSINE-5'-DIPHOSPHATE",
    "ADENOSINE-5'-TRIPHOSPHATE",

    "CYTIDINE",
    "CYTIDINE-5'-MONOPHOSPHATE",
    "CYTIDINE-5'-DIPHOSPHATE",
    "CYTIDINE-5'-TRIPHOSPHATE",

    "THYMIDINE",
    "THYMIDINE-5'-MONOPHOSPHATE",
    "THYMIDINE-5'-DIPHOSPHATE",
    "THYMIDINE-5'-TRIPHOSPHATE",

    "GUANOSINE",
    "GUANOSINE-5'-MONOPHOSPHATE",
    "GUANOSINE-5'-DIPHOSPHATE",
    "GUANOSINE-5'-TRIPHOSPHATE",

    "URIDINE",
    "URIDINE-5'-MONOPHOSPHATE",
    "URIDINE-5'-DIPHOSPHATE",
    "URIDINE-5'-TRIPHOSPHATE"
])

# -----------------------------
# 3. Normalize names
# -----------------------------
def normalize(name):
    if pd.isna(name):
        return ""

    name = str(name).upper().strip()
    name = re.sub(r"\s+", " ", name)
    name = name.replace("5' -", "5'-")
    name = name.replace("5' ", "5'-")
    name = name.replace("’", "'")

    return name

df["Name_norm"] = df["Name"].apply(normalize)

# -----------------------------
# 4. Set modification columns
# -----------------------------
mod_cols = [
    "Base Modified",
    "Sugar Modified",
    "Phosphate Modified"
]

# create columns if missing
for col in mod_cols:
    if col not in df.columns:
        df[col] = "No"

# clean values
for col in mod_cols:
    df[col] = (
        df[col]
        .fillna("No")
        .astype(str)
        .str.strip()
    )

# -----------------------------
# 4b. Classify modified/unmodified
# -----------------------------
mask_modified = (
    (df["Base Modified"] == "Yes") |
    (df["Sugar Modified"] == "Yes") |
    (df["Phosphate Modified"] == "Yes")
)

df["class"] = mask_modified.map({
    True: "modified",
    False: "unmodified"
})

# -----------------------------
# 5. Calculate molecular weight
# -----------------------------
def calc_mw(smiles):
    if pd.isna(smiles):
        return None

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None

    return Descriptors.MolWt(mol)

df["MW"] = df["SMILES"].apply(calc_mw)

# remove failed calculations
valid_mw = df["MW"].dropna()

mw_min = valid_mw.min()
mw_max = valid_mw.max()
mw_mean = valid_mw.mean()

# -----------------------------
# 6. Counts
# -----------------------------
total = len(df)
n_unmod = (df["class"] == "unmodified").sum()
n_mod = (df["class"] == "modified").sum()

# count Polymer == Yes
n_polymer = (
    df["Polymer"]
    .fillna("")
    .astype(str)
    .str.strip()
    .eq("Yes")
    .sum()
)

# count Standalone == Yes
n_standalone = (
    df["Standalone"]
    .fillna("")
    .astype(str)
    .str.strip()
    .eq("Yes")
    .sum()
)

n_polymer_modified = (
    (df["Polymer"].fillna("").astype(str).str.strip() == "Yes") &
    (df["class"] == "modified")
).sum()

n_standalone_modified = (
    (df["Standalone"].fillna("").astype(str).str.strip() == "Yes") &
    (df["class"] == "modified")
).sum()


print("\n=== SUMMARY ===")
print(f"Total nucleotides: {total}")
print(f"Unmodified:        {n_unmod}")
print(f"Modified:          {n_mod}")
print(f"Polymer = Yes:     {n_polymer}")
print(f"Polymer = Yes & Modified: {n_polymer_modified}")
print(f"Standalone = Yes & Modified: {n_standalone_modified}")
print(f"Standalone = Yes:  {n_standalone}")

print("\n=== MOLECULAR WEIGHT ===")
print(f"Lowest MW:         {mw_min:.2f} Da")
print(f"Highest MW:        {mw_max:.2f} Da")
print(f"Mean MW:           {mw_mean:.2f} Da")

# identify actual compounds
low_row = df.loc[df["MW"].idxmin()]
high_row = df.loc[df["MW"].idxmax()]

print("\nLowest MW compound:")
print(f"{low_row['CCD_ID']} : {low_row['Name']} ({low_row['MW']:.2f} Da)")

print("\nHighest MW compound:")
print(f"{high_row['CCD_ID']} : {high_row['Name']} ({high_row['MW']:.2f} Da)")

# -----------------------------
# 7. Save output
# -----------------------------
df[[
    "CCD_ID",
    "Name",
    "SMILES",
    "class",
    "MW",
    "Base Modified",
    "Sugar Modified",
    "Phosphate Modified"
]].to_csv(
    "CCD_nucleotide_classification_with_MW.csv",
    index=False
)

print("\nSaved: CCD_nucleotide_classification_with_MW.csv")
