#!/usr/bin/env python3
# usage: python clustering.py topol.tpr traj.xtc
'''Visualisation of clustering results.'''

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # non-interactive
import seaborn as sns
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import hdbscan

# Load closed trajectory
tpr = sys.argv[1]
trj = sys.argv[2]
u = mda.Universe(tpr, trj)

# Define tail and barrel selections
tail = u.select_atoms("resid 213:221")
barrel = u.select_atoms("protein and not resid 213:221")
n_tail = len(tail.residues)
n_barrel = len(barrel.residues)

# Align trajectory to core (TIM barrel-like region)
align.AlignTraj(u, u, select="backbone and not resid 213:221", in_memory=False).run()

# Downsample trajectory for clustering
stride = 10
max_frames = 5000
sampled_frames = []
times = []

for i, ts in enumerate(u.trajectory[::stride]):
    if i >= max_frames:
        break
    sampled_frames.append(ts.frame)
    times.append(ts.time)

n_frames = len(sampled_frames)
print(f"Using {n_frames} frames for clustering")

# Build RMSD distance matrices
def compute_rmsd_matrix(universe, selection, frame_indices):
    sel = universe.select_atoms(selection)

    coords = []

    for idx in frame_indices:
        universe.trajectory[idx]
        coords.append(sel.positions.copy())

    coords = np.array(coords)
    n = len(coords)

    print(f"Building RMSD matrix of size {n} x {n}")

    dist_matrix = np.zeros((n, n), dtype=np.float64)

    for i in range(n):
        for j in range(i + 1, n):
            val = rmsd(coords[i], coords[j],
                       center=True, superposition=True)
            dist_matrix[i, j] = dist_matrix[j, i] = val

    return dist_matrix

# Whole protein
dist_protein = compute_rmsd_matrix(u, "protein", sampled_frames)
np.save("dist_protein.npy", dist_protein)

# Tail only
dist_tail = compute_rmsd_matrix(u, "resid 213:221", sampled_frames)
np.save("dist_tail.npy", dist_tail)

# Hierarchical clustering
def hierarchical_clustering(dist_matrix, cutoff=2.5):
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')
    labels = fcluster(Z, cutoff, criterion='distance')
    return labels, Z

labels_protein_1, Zp = hierarchical_clustering(dist_protein, cutoff=1)
labels_tail_1, Zt = hierarchical_clustering(dist_tail, cutoff=1)
labels_protein_2_5, Zp = hierarchical_clustering(dist_protein, cutoff=2.5)
labels_tail_2_5, Zt = hierarchical_clustering(dist_tail, cutoff=2.5)
labels_protein_2, Zp = hierarchical_clustering(dist_protein, cutoff=2)
labels_tail_2, Zt = hierarchical_clustering(dist_tail, cutoff=2)
labels_protein_3, Zp = hierarchical_clustering(dist_protein, cutoff=3)
labels_tail_3, Zt = hierarchical_clustering(dist_tail, cutoff=3)

# HDBSCAN
def run_hdbscan(dist_matrix):
    clusterer = hdbscan.HDBSCAN(
        metric='precomputed',
        min_cluster_size=50,
        min_samples=10
    )
    labels = clusterer.fit_predict(dist_matrix)
    return labels, clusterer

hdb_labels_protein, hdb_prot = run_hdbscan(dist_protein)
hdb_labels_tail, hdb_tail = run_hdbscan(dist_tail)

records = []

def store(labels, method, selection, cutoff):
    for i, lab in enumerate(labels):
        records.append({
            "index": i,
            "frame": sampled_frames[i],
            "time_ps": times[i],
            "method": method,
            "selection": selection,
            "cutoff": cutoff,
            "cluster": int(lab)
        })

# Store hierarchical results
store(labels_protein_1, "hierarchical", "protein", 1.0)
store(labels_tail_1, "hierarchical", "tail", 1.0)

store(labels_protein_2, "hierarchical", "protein", 2.0)
store(labels_tail_2, "hierarchical", "tail", 2.0)

store(labels_protein_2_5, "hierarchical", "protein", 2.5)
store(labels_tail_2_5, "hierarchical", "tail", 2.5)

store(labels_protein_3, "hierarchical", "protein", 3.0)
store(labels_tail_3, "hierarchical", "tail", 3.0)

# Store HDBSCAN
store(hdb_labels_protein, "hdbscan", "protein", None)
store(hdb_labels_tail, "hdbscan", "tail", None)

df = pd.DataFrame(records)
df.to_csv("cluster_assignments.csv", index=False)

# compute medoids for all sets of labels
all_label_sets = [
    ("hierarchical", "protein", 1.0, labels_protein_1, dist_protein),
    ("hierarchical", "tail", 1.0, labels_tail_1, dist_tail),

    ("hierarchical", "protein", 2.0, labels_protein_2, dist_protein),
    ("hierarchical", "tail", 2.0, labels_tail_2, dist_tail),

    ("hierarchical", "protein", 2.5, labels_protein_2_5, dist_protein),
    ("hierarchical", "tail", 2.5, labels_tail_2_5, dist_tail),

    ("hierarchical", "protein", 3.0, labels_protein_3, dist_protein),
    ("hierarchical", "tail", 3.0, labels_tail_3, dist_tail),

    ("hdbscan", "protein", None, hdb_labels_protein, dist_protein),
    ("hdbscan", "tail", None, hdb_labels_tail, dist_tail),
]

def compute_medoids(dist_matrix, labels):
    medoids = {}
    for c in np.unique(labels):
        if c == -1:
            continue
        idx = np.where(labels == c)[0]
        sub = dist_matrix[np.ix_(idx, idx)]
        medoids[c] = idx[np.argmin(sub.sum(axis=1))]
    return medoids

medoid_records = []
def store_medoids(method, selection, cutoff, labels, dist):
    medoids = compute_medoids(dist, labels)

    for c, idx in medoids.items():
        medoid_records.append({
            "method": method,
            "selection": selection,
            "cutoff": cutoff,
            "cluster": int(c),
            "frame_index": int(idx),
            "frame": int(sampled_frames[idx]),
            "time_ps": float(times[idx])
        })


for method, selection, cutoff, labels, dist in all_label_sets:
    store_medoids(method, selection, cutoff, labels, dist)

df_medoids = pd.DataFrame(medoid_records)
df_medoids.to_csv("medoids.csv", index=False)

### Plotting
df_clust = pd.read_csv("cluster_assignments.csv")

# Load PCA data
df_pca = pd.read_csv("pca_projection.dat", delim_whitespace=True,
                     comment="#", names=["time", "PC1", "PC2"])

pc1_all = df_pca["PC1"].values
pc2_all = df_pca["PC2"].values
time_all = df_pca["time"].values

# variance
df2 = pd.read_csv("pc_variance.csv")
pc1_var = df2[df2["PC"] == 1]["variance_percent"].values[0]
pc2_var = df2[df2["PC"] == 2]["variance_percent"].values[0]

# medoids
df_medoids = pd.read_csv("medoids.csv")

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
        
    # overlay medoids
    sub_medoids = df_medoids[
        (df_medoids["method"] == method) &
        (df_medoids["selection"] == selection) &
        (df_medoids["cutoff"].fillna(-999) == (cutoff if cutoff == cutoff else -999))
    ]

    for _, row in sub_medoids.iterrows():
        x = pc1_all[frames]
        y = pc2_all[frames]

        i = cluster_map[row["cluster"]] - 1

        plt.scatter(x, y,
                    color=colors[i],
                    s=80,
                    marker="*",
                    edgecolor="black",
                    linewidth=0.7,
                    zorder=10)

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