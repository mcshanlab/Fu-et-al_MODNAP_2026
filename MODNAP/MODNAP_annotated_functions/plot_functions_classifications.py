import matplotlib

# Ensure editable text in PDF
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

import matplotlib.pyplot as plt
from collections import Counter

# --------------------------------------
# Input / Output
# --------------------------------------
input_file = "../MODNAP_entry_list.txt"
output_pdf = "pdb_function_histogram.pdf"
output_png = "pdb_function_histogram.png"
output_txt = "pdb_function_counts.txt"
output_pdf_vertical = "pdb_function_histogram_vertical.pdf"
output_png_vertical = "pdb_function_histogram_vertical.png"

# --------------------------------------
# Parse file and extract functions
# --------------------------------------
functions = []

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        # Function is everything AFTER first column
        # Handles multi-word labels like "protein-DNA Transcription"
        func = " ".join(parts[1:])
        functions.append(func)

# --------------------------------------
# Count occurrences
# --------------------------------------
counts = Counter(functions)

# Sort by count (descending)
sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

labels = [x[0] for x in sorted_counts]
values = [x[1] for x in sorted_counts]

# --------------------------------------
# Output stats
# --------------------------------------
total_entries = len(functions)

print(f"Total entries: {total_entries}")
print("\nCounts per function:")

with open(output_txt, "w") as f:
    f.write(f"Total entries: {total_entries}\n\n")
    f.write("Function\tCount\tPercentage\n")

    for func, count in sorted_counts:
        pct = count / total_entries * 100
        print(f"{func}: {count} ({pct:.2f}%)")
        f.write(f"{func}\t{count}\t{pct:.2f}%\n")

print(f"\nSaved counts to: {output_txt}")

# --------------------------------------
# Plot (horizontal bar = histogram style)
# --------------------------------------
plt.figure(figsize=(6, 4 + 0.3 * len(labels)))

light_colors = [
    "#A7C7E7", "#B8E0D2", "#F9D5E5", "#FFF1B6", "#D6CDEA",
    "#F7C6C7", "#CDEAC0", "#FFD6A5", "#B5EAD7", "#E2F0CB"
]

# Repeat colors if there are more bars than colors
colors = [light_colors[i % len(light_colors)] for i in range(len(labels))]

bars = plt.barh(labels, values, color=colors)

plt.xlabel("Count")
plt.ylabel("Classification")

# Invert y-axis so highest count is on top
plt.gca().invert_yaxis()

# Add a bit of room on the right for labels
plt.xlim(0, max(values) * 1.15)

# --------------------------------------
# Add percentage labels to the right
# --------------------------------------
for i, val in enumerate(values):
    pct = val / total_entries * 100
    plt.text(
        val + max(values) * 0.01,   # small offset
        i,
        f"{pct:.1f}%",
        va='center',
        fontsize=9
    )

# --------------------------------------
# Clean style
# --------------------------------------
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

plt.tight_layout()

# --------------------------------------
# Save PDF
# --------------------------------------
plt.savefig(output_pdf, dpi=300)
plt.savefig(output_png, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf}")
print(f"Saved figure to: {output_png}")
# --------------------------------------
# SECOND FIGURE: vertical bars (axes reversed)
# --------------------------------------
plt.figure(figsize=(6 + 0.3 * len(labels), 5))

# Use same colors
bars = plt.bar(labels, values, color=colors)

plt.xlabel("Classification")
plt.ylabel("Count")

# Rotate x labels so they don't overlap
plt.xticks(rotation=45, ha="right")

# Add percentage labels above bars
for i, val in enumerate(values):
    pct = val / total_entries * 100
    plt.text(
        i,
        val + max(values) * 0.01,
        f"{pct:.1f}%",
        ha='center',
        fontsize=9
    )

# Clean style
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

plt.tight_layout()

# Save vertical versions
plt.savefig(output_pdf_vertical, dpi=300)
plt.savefig(output_png_vertical, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf_vertical}")
print(f"Saved figure to: {output_png_vertical}")
