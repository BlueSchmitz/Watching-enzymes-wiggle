#!/usr/bin/env python3
# usage: python clustering.py topol.tpr traj.xtc
'''Visualisation of hbonds, hydrophobic contacts and salt bridges.'''

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
tail = u.select_atoms("resid 249:259")
barrel = u.select_atoms("protein and not resid 249:259")
n_tail = len(tail.residues)
n_barrel = len(barrel.residues)

# Align trajectory to core (TIM barrel-like region)
align.AlignTraj(u, u, select="backbone and not resid 249:259", in_memory=False).run()

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

# Tail only
dist_tail = compute_rmsd_matrix(u, "resid 150:200", sampled_frames)

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

times = [ts.time for ts in u.trajectory]

def store(labels, method, selection, cutoff):
    for i, lab in enumerate(labels):
        records.append({
            "frame": i,
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

# Plotting
df = pd.read_csv("pca_projection.dat", delim_whitespace=True)
pc1 = df["PC1"].values
pc2 = df["PC2"].values

df2 = pd.read_csv("pc_variance.csv")
pc1_var = df2[df2["PC"] == 1]["variance_percent"].values[0]
pc2_var = df2[df2["PC"] == 2]["variance_percent"].values[0]

# plot pc1 vs pc2 colored by cluster
def plot_clustering(labels, method):
    assert len(labels) == len(pc1)
    plt.figure(figsize=(6, 4))
    plt.scatter(pc1, pc2, s=5, c=labels,cmap="viridis", alpha=0.7)
    plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
    plt.tight_layout()
    plt.savefig(f"pc1_vs_pc2_clusters_{method}.pdf")
    plt.close()

plots = [
    (labels_protein_1, "protein_hierarchical_1"),
    (labels_tail_1, "tail_hierarchical_1"),
    (labels_protein_2, "protein_hierarchical_2"),
    (labels_tail_2, "tail_hierarchical_2"),
    (labels_protein_2_5, "protein_hierarchical_2_5"),
    (labels_tail_2_5, "tail_hierarchical_2_5"),
    (labels_protein_3, "protein_hierarchical_3"),
    (labels_tail_3, "tail_hierarchical_3"),
    (hdb_labels_protein, "protein_hdbscan"),
    (hdb_labels_tail, "tail_hdbscan")
]

for labels, name in plots:
    plot_clustering(labels, name)