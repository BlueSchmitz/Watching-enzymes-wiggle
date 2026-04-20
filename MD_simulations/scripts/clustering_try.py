#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import pandas as pd

def compute_medoids(pc1, pc2, labels):
    medoids = {}

    for c in np.unique(labels):
        if c == -1:
            continue

        mask = labels == c
        idx = np.where(mask)[0]

        x = pc1[mask]
        y = pc2[mask]

        cx = x.mean()
        cy = y.mean()

        dists = (x - cx)**2 + (y - cy)**2
        medoid_local = np.argmin(dists)
        medoid_global = idx[medoid_local]

        medoids[c] = medoid_global

    return medoids

# Load data
df_pca = pd.read_csv("pca_projection.dat", delim_whitespace=True, comment="#", names=["time", "PC1", "PC2"])
df_clust = pd.read_csv("cluster_assignments_try.csv")

# Extract arrays
pc1_all = df_pca["PC1"].values
pc2_all = df_pca["PC2"].values
time_all = df_pca["time"].values

frames = df_clust["frame"].values
labels = df_clust["cluster"].values

# Select only clustered frames
pc1 = pc1_all[frames]
pc2 = pc2_all[frames]
times = time_all[frames]

# variance
df2 = pd.read_csv("pc_variance.csv")
pc1_var = df2[df2["PC"] == 1]["variance_percent"].values[0]
pc2_var = df2[df2["PC"] == 2]["variance_percent"].values[0]

# map clusters
unique, counts = np.unique(labels, return_counts=True)
cluster_sizes = dict(zip(unique, counts))
# keep only top 10 clusters by size
selected_clusters = sorted(cluster_sizes, key=cluster_sizes.get, reverse=True)[:10]
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(selected_clusters)))

# new mapping for consistent cluster numbering
cluster_map = {old: new for new, old in enumerate(selected_clusters, start=1)}
# Apply mapping
mapped_labels = np.array([cluster_map.get(l, -1) for l in labels])

# Plot
plt.figure(figsize=(5, 4))
# Plot noise points in light grey
noise_mask = mapped_labels == -1
plt.scatter(pc1[noise_mask], pc2[noise_mask],
            color="lightgrey",
            s=5,
            label="-1",
            alpha=0.5)
# Plot clusters 1–10
for i in range(1, len(selected_clusters) + 1):
    mask = mapped_labels == i
    plt.scatter(pc1[mask], pc2[mask],
                color=colors[i-1],
                s=5,
                label=f"{i}",
                alpha=0.7)
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.legend(frameon=True, title="Clusters", fontsize=6)
plt.tight_layout()
plt.savefig("pc1_vs_pc2_clusters.pdf")
plt.close()
