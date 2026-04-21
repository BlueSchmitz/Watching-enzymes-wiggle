#!/usr/bin/env python3
# usage: python plot_clustering.py
'''Visualisation of clustering results.'''

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # non-interactive
import seaborn as sns
import pandas as pd

### Plotting
df_clust = pd.read_csv("cluster_assignments.csv")

# Load PCA data
df_pca = pd.read_csv("pca_projection.dat", delim_whitespace=True,
                     comment="#", names=["time", "PC1", "PC2"])

pc1_all = df_pca["PC1"].values
pc2_all = df_pca["PC2"].values
time_all = df_pca["time"].values

# Load variance data
df2 = pd.read_csv("pc_variance.csv")
pc1_var = df2[df2["PC"] == 1]["variance_percent"].values[0]
pc2_var = df2[df2["PC"] == 2]["variance_percent"].values[0]

for (method, selection, cutoff), subdf in df_clust.groupby(["method", "selection", "cutoff"], dropna=False):

    frames = subdf["frame"].values
    labels = subdf["cluster"].values
    times = subdf["time_ps"].values

    # Map PCA values
    pc1 = pc1_all[frames]
    pc2 = pc2_all[frames]

    # Cluster sizes
    unique, counts = np.unique(labels, return_counts=True)
    cluster_sizes = dict(zip(unique, counts))
    cluster_sizes.pop(-1, None) # remove noise from ranking

    selected_clusters = sorted(cluster_sizes,
                               key=cluster_sizes.get,
                               reverse=True)[:10]
    
    # colors
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0, 1, len(selected_clusters)))

    # mapping
    cluster_map = {old: new for new, old in enumerate(selected_clusters, start=1)}
    mapped_labels = np.array([cluster_map.get(l, -1) for l in labels])

    # Plotting
    plt.figure(figsize=(5, 4))

    # noise
    noise_mask = mapped_labels == -1
    plt.scatter(pc1[noise_mask], pc2[noise_mask],
                color="lightgrey", s=5, alpha=0.5, label="-1")

    # clusters 1–10
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

    # Filename
    cutoff_str = "none" if pd.isna(cutoff) else str(cutoff).replace(".", "_")

    fname = f"pc1_pc2_{method}_{selection}_cutoff_{cutoff_str}.pdf"
    plt.savefig(fname)
    plt.close()

    print(f"Saved {fname}")