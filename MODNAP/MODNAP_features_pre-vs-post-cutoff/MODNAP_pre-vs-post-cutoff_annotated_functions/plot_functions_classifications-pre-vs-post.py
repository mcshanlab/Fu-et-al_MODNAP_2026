import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

# --------------------------------------
# INPUT FILES
# --------------------------------------
entry_file = "../../MODNAP_entry_list.txt"
pre_file = "../pre_cutoff.txt"
post_file = "../post_cutoff.txt"

# --------------------------------------
# STEP 1: Load PDB → function mapping
# --------------------------------------
pdb_to_function = {}

with open(entry_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        full_id = parts[0]
        pdb = full_id.split("-")[1]
        func = " ".join(parts[1:])

        pdb_to_function[pdb] = func

# --------------------------------------
# STEP 2: Load pre/post PDB sets
# --------------------------------------
def load_pdbs(filepath):
    pdbs = set()
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("CCD"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                pdbs.add(parts[1])
    return pdbs

pre_pdbs = load_pdbs(pre_file)
post_pdbs = load_pdbs(post_file)

print("Pre entries:", len(pre_pdbs))
print("Post entries:", len(post_pdbs))

# --------------------------------------
# STEP 3: Map PDB → functions
# --------------------------------------
pre_functions = []
post_functions = []

for pdb, func in pdb_to_function.items():
    if pdb in pre_pdbs:
        pre_functions.append(func)
    elif pdb in post_pdbs:
        post_functions.append(func)

print("Pre functions:", len(pre_functions))
print("Post functions:", len(post_functions))

pre_counts = Counter(pre_functions)
post_counts = Counter(post_functions)

# --------------------------------------
# STEP 4: SORT BY PRE COUNTS
# --------------------------------------
labels = sorted(pre_counts.keys(), key=lambda x: pre_counts[x], reverse=True)

for k in post_counts:
    if k not in labels:
        labels.append(k)

pre_values = [pre_counts.get(x, 0) for x in labels]
post_values = [post_counts.get(x, 0) for x in labels]

# --------------------------------------
# OUTPUT TABLE
# --------------------------------------
with open("pre_vs_post_function_counts.txt", "w") as f:
    f.write("Function\tPre\tPost\n")
    for l, p1, p2 in zip(labels, pre_values, post_values):
        f.write(f"{l}\t{p1}\t{p2}\n")

print("Saved table")

# --------------------------------------
# PLOT (CORRECT: vertical grouped bars)
# --------------------------------------
x = np.arange(len(labels))
width = 0.4

plt.figure(figsize=(max(10, 0.5 * len(labels)), 6))

plt.bar(x - width/2, pre_values, width=width, label="Pre", color="#A7C7E7")
plt.bar(x + width/2, post_values, width=width, label="Post", color="#F7C6C7")

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Count")
plt.xlabel("Function")

plt.legend()

plt.tight_layout()

plt.savefig("pre_vs_post_function_histogram.pdf", dpi=300)
plt.savefig("pre_vs_post_function_histogram.png", dpi=300)
plt.close()

print("Saved correct plot")
