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
occupancy = counts / counts.sum()
cluster_occ = dict(zip(unique, occupancy))
# keep only clusters >10%
selected_clusters = [c for c, occ in cluster_occ.items() if occ > 0.10]
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, len(selected_clusters)))

# Plot
plt.figure(figsize=(5, 4))
for i, c in enumerate(selected_clusters):
    mask = labels == c
    plt.scatter(pc1[mask], pc2[mask],
                color=colors[i],
                s=5,
                label=f"{c}",
                alpha=0.7)
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.legend(frameon=True, title="Clusters", fontsize=6)
plt.tight_layout()
plt.savefig("pc1_vs_pc2_clusters.pdf")
plt.close()
