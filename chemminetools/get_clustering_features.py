import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics import pairwise_distances

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

# =============================
# GLOBAL STYLE (Arial + PDF vector)
# =============================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# -----------------------------
# 1. Load dataset
# -----------------------------
df = pd.read_csv("HierarchicalClustering_canvas.tsv", sep="\t")

if not np.issubdtype(df.iloc[:, 0].dtype, np.number):
    ids = df.iloc[:, 0]
    X = df.iloc[:, 1:].copy()
else:
    ids = None
    X = df.copy()

X = X.apply(pd.to_numeric, errors="coerce")

# -----------------------------
# 2. Drop fully empty columns
# -----------------------------
empty_cols = X.columns[X.isna().all()].tolist()

if len(empty_cols) > 0:
    print("\nDropping fully empty features:")
    for c in empty_cols:
        print(" -", c)
    X = X.drop(columns=empty_cols)

# -----------------------------
# 3. Impute missing values
# -----------------------------
imputer = SimpleImputer(strategy="median")
X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=X.columns)

# -----------------------------
# 4. Standardize
# -----------------------------
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X_imputed), columns=X.columns)

# -----------------------------
# 5. PCA
# -----------------------------
n_components = min(10, X_scaled.shape[1])
pca = PCA(n_components=n_components)
pca.fit(X_scaled)

loadings = pd.DataFrame(
    pca.components_.T,
    index=X_scaled.columns,
    columns=[f"PC{i+1}" for i in range(n_components)]
)

explained = pca.explained_variance_ratio_

pca_importance = (loadings.abs() * explained).sum(axis=1)

# -----------------------------
# 6. Correlation redundancy
# -----------------------------
corr = X_scaled.corr().abs()
redundancy = corr.sum() - 1

# -----------------------------
# 7. Mutual information
# -----------------------------
pc1_scores = pca.transform(X_scaled)[:, 0]
mi = mutual_info_regression(X_scaled, pc1_scores)
mi_series = pd.Series(mi, index=X_scaled.columns)

# -----------------------------
# 8. Combine importance
# -----------------------------
importance = pd.DataFrame({
    "PCA_importance": pca_importance,
    "Redundancy_penalty": redundancy,
    "Mutual_information": mi_series
})

importance_norm = (
    importance - importance.min()
) / (importance.max() - importance.min() + 1e-12)

importance["Composite_score"] = (
    0.5 * importance_norm["PCA_importance"]
    + 0.2 * importance_norm["Mutual_information"]
    - 0.3 * importance_norm["Redundancy_penalty"]
)

importance_sorted = importance.sort_values("Composite_score", ascending=False)

print("\nTop physicochemical features driving clustering:\n")
print(importance_sorted["Composite_score"].head(20))

importance_sorted.to_csv("feature_importance_clustering.csv")

# -----------------------------
# 9. Hierarchical clustering
# -----------------------------
dist = pairwise_distances(X_scaled, metric="euclidean")
Z = linkage(dist, method="average")

cluster_labels = fcluster(Z, t=5, criterion="maxclust")

pc = pca.transform(X_scaled)

pc1_var = explained[0] * 100
pc2_var = explained[1] * 100

# =============================
# FIGURE
# =============================
fig = plt.figure(figsize=(18, 12))

# -----------------------------
# A. Feature importance
# -----------------------------
ax1 = plt.subplot(2, 2, 1)

top_feats = importance_sorted.head(15)

ax1.barh(
    top_feats.index[::-1],
    top_feats["Composite_score"][::-1]
)
ax1.set_title("A. Key physicochemical drivers of clustering")
ax1.set_xlabel("Composite score")

# -----------------------------
# B. PCA clustering space (UNCHANGED)
# -----------------------------
ax2 = plt.subplot(2, 2, 2)

unique_clusters = np.unique(cluster_labels)
cluster_sizes = {c: np.sum(cluster_labels == c) for c in unique_clusters}
sorted_clusters = sorted(unique_clusters, key=lambda c: cluster_sizes[c], reverse=True)

cmap = plt.get_cmap("tab10")

for i, c in enumerate(sorted_clusters):
    mask = cluster_labels == c
    ax2.scatter(
        pc[mask, 0],
        pc[mask, 1],
        label=f"Cluster {c}",
        s=40,
        color=cmap(i % 10),
        zorder=2 if i else 1
    )

ax2.set_title("B. PCA chemical space (clustered)")
ax2.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
ax2.set_ylabel(f"PC2 ({pc2_var:.1f}%)")
ax2.legend(frameon=False, fontsize=8)

# -----------------------------
# C. Correlation heatmap
# -----------------------------
ax3 = plt.subplot(2, 2, 3)

top_corr_feats = importance_sorted.head(12).index
corr_top = X_scaled[top_corr_feats].corr()

sns.heatmap(
    corr_top,
    cmap="coolwarm",
    center=0,
    square=True,
    ax=ax3
)

ax3.set_title("C. Descriptor redundancy structure")
ax3.text(
    0.5, 1.08,
    "Pearson correlation coefficient matrix (red = +, blue = -)",
    transform=ax3.transAxes,
    ha="center",
    fontsize=9
)

# -----------------------------
# D. DENDROGRAM (FIXED COLORING)
# -----------------------------
ax4 = plt.subplot(2, 2, 4)

# IMPORTANT FIX: use a real cutoff so colors appear
max_d = Z[:, 2].max()
color_threshold = 0.7 * max_d

dendrogram(
    Z,
    no_labels=True,
    color_threshold=color_threshold,
    ax=ax4
)

ax4.set_title("D. Hierarchical clustering dendrogram")
ax4.set_ylabel("Distance")

ax4.text(
    0.5, 1.08,
    "Branch colors indicate clusters at cutoff distance",
    transform=ax4.transAxes,
    ha="center",
    fontsize=9
)

plt.tight_layout()

# -----------------------------
# SAVE PDF (Illustrator-ready)
# -----------------------------
plt.savefig("Figure_clustering_analysis.pdf", bbox_inches="tight")
plt.show()

print("\nSaved: Figure_clustering_analysis.pdf")
