#!/usr/bin/env python3
# usage: python run_clustering_Bb.py topol.tpr traj.xtc
'''Run clustering.'''

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
tail = u.select_atoms("resid 220:229")
barrel = u.select_atoms("protein and not resid 220:229")
n_tail = len(tail.residues)
n_barrel = len(barrel.residues)

# Align trajectory to core (TIM barrel-like region)
align.AlignTraj(u, u, select="backbone and not resid 220:229", in_memory=False).run()

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
dist_tail = compute_rmsd_matrix(u, "resid 220:229", sampled_frames)
np.save("dist_tail.npy", dist_tail)

# Hierarchical clustering
def hierarchical_clustering(dist_matrix, cutoff=2.5):
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')
    labels = fcluster(Z, cutoff, criterion='distance')
    return labels, Z

labels_protein_2, Zp = hierarchical_clustering(dist_protein, cutoff=2)
labels_tail_2, Zt = hierarchical_clustering(dist_tail, cutoff=2)
labels_protein_2_5, Zp = hierarchical_clustering(dist_protein, cutoff=2.5)
labels_tail_2_5, Zt = hierarchical_clustering(dist_tail, cutoff=2.5)
labels_protein_3, Zp = hierarchical_clustering(dist_protein, cutoff=3)
labels_tail_3, Zt = hierarchical_clustering(dist_tail, cutoff=3)
labels_protein_3_5, Zp = hierarchical_clustering(dist_protein, cutoff=3.5)
labels_tail_3_5, Zt = hierarchical_clustering(dist_tail, cutoff=3.5)
labels_protein_4, Zp = hierarchical_clustering(dist_protein, cutoff=4)
labels_tail_4, Zt = hierarchical_clustering(dist_tail, cutoff=4)
labels_tail_5, Zt = hierarchical_clustering(dist_tail, cutoff=5)
labels_tail_6, Zt = hierarchical_clustering(dist_tail, cutoff=6)

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
store(labels_protein_2, "hierarchical", "protein", 2.0)
store(labels_tail_2, "hierarchical", "tail", 2.0)

store(labels_protein_2_5, "hierarchical", "protein", 2.5)
store(labels_tail_2_5, "hierarchical", "tail", 2.5)

store(labels_protein_3, "hierarchical", "protein", 3.0)
store(labels_tail_3, "hierarchical", "tail", 3.0)

store(labels_protein_3_5, "hierarchical", "protein", 3.5)
store(labels_tail_3_5, "hierarchical", "tail", 3.5)

store(labels_protein_4, "hierarchical", "protein", 4.0)
store(labels_tail_4, "hierarchical", "tail", 4.0)

store(labels_tail_5, "hierarchical", "tail", 5.0)
store(labels_tail_6, "hierarchical", "tail", 6.0)

# Store HDBSCAN
store(hdb_labels_protein, "hdbscan", "protein", None)
store(hdb_labels_tail, "hdbscan", "tail", None)

df = pd.DataFrame(records)
df.to_csv("cluster_assignments.csv", index=False)

# compute medoids for all sets of labels
all_label_sets = [
    ("hierarchical", "protein", 2.0, labels_protein_2, dist_protein),
    ("hierarchical", "tail", 2.0, labels_tail_2, dist_tail),

    ("hierarchical", "protein", 2.5, labels_protein_2_5, dist_protein),
    ("hierarchical", "tail", 2.5, labels_tail_2_5, dist_tail),

    ("hierarchical", "protein", 3.0, labels_protein_3, dist_protein),
    ("hierarchical", "tail", 3.0, labels_tail_3, dist_tail),

    ("hierarchical", "protein", 3.5, labels_protein_3_5, dist_protein),
    ("hierarchical", "tail", 3.5, labels_tail_3_5, dist_tail),

    ("hierarchical", "protein", 4.0, labels_protein_4, dist_protein),
    ("hierarchical", "tail", 4.0, labels_tail_4, dist_tail),

    ("hierarchical", "tail", 5.0, labels_tail_5, dist_tail),
    ("hierarchical", "tail", 6.0, labels_tail_6, dist_tail),

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