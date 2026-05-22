#!/usr/bin/env python3
# usage:
# python plot_distance_distributions.py topol.tpr traj.xtc

"""
Generate distance distribution plots for:
- hydrogen bonds
- hydrophobic contacts

Uses previously saved:
- hbonds.csv
- hydrophobic_contacts_matrix.csv
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# ---------------------------------------------------------
# Input
# ---------------------------------------------------------

tpr = sys.argv[1]
trj = sys.argv[2]

u = mda.Universe(tpr, trj)

# ---------------------------------------------------------
# Load saved data
# ---------------------------------------------------------

hbonds_df = pd.read_csv("hbonds.csv")

hydrophobic_data = pd.read_csv(
    "hydrophobic_contacts_matrix.csv",
    index_col=0
)

# ---------------------------------------------------------
# Residue labels
# ---------------------------------------------------------

res_lookup = {
    r.resid: f"{r.resname}{r.resid}"
    for r in u.residues
}

# ---------------------------------------------------------
# Output directories
# ---------------------------------------------------------

os.makedirs("hbond_kdes", exist_ok=True)
os.makedirs("hydrophobic_kdes", exist_ok=True)

HBOND_CUTOFF = 3.5
HYDROPHOBIC_CUTOFF = 4.5

# ---------------------------------------------------------
# H-bond distance distributions
# ---------------------------------------------------------

print("Generating H-bond distance distributions...")

hbond_pairs = (
    hbonds_df[[
        "tail_resid",
        "barrel_resid",
        "donor_idx",
        "acceptor_idx"
    ]]
    .drop_duplicates()
)

for _, row in hbond_pairs.iterrows():

    tail_resid = int(row["tail_resid"])
    barrel_resid = int(row["barrel_resid"])

    donor_idx = int(row["donor_idx"])
    acceptor_idx = int(row["acceptor_idx"])

    distances_over_time = []

    for ts in u.trajectory:

        d = distances.distance_array(
            u.atoms[[donor_idx]].positions,
            u.atoms[[acceptor_idx]].positions,
            box=u.dimensions
        )[0, 0]

        distances_over_time.append(d)

    distances_over_time = np.array(distances_over_time)

    plt.figure(figsize=(5, 5))

    sns.histplot(
        distances_over_time,
        bins=50,
        stat="density",
        kde=True
    )

    plt.axvline(
        HBOND_CUTOFF,
        color="grey",
        linestyle=":",
        linewidth=2
    )

    plt.xlim(0, 10)

    plt.xlabel("Donor-acceptor distance (Å)")
    plt.ylabel("Density")

    plt.title(
        f"{res_lookup[tail_resid]} ↔ "
        f"{res_lookup[barrel_resid]}"
    )

    plt.tight_layout()

    outfile = (
        f"hbond_kdes/"
        f"{res_lookup[tail_resid]}_"
        f"{res_lookup[barrel_resid]}.pdf"
    )

    plt.savefig(outfile)
    plt.close()

# ---------------------------------------------------------
# Hydrophobic selections
# ---------------------------------------------------------

tail = u.select_atoms("resid 249:259")
barrel = u.select_atoms("protein and not resid 249:259")

hydrophobic_sel = (
    "resname ALA VAL LEU ILE MET PHE TRP TYR PRO "
    "and not name H*"
)

tail_hphob = tail.select_atoms(hydrophobic_sel)
barrel_hphob = barrel.select_atoms(hydrophobic_sel)

# ---------------------------------------------------------
# Hydrophobic distance distributions
# ---------------------------------------------------------

print("Generating hydrophobic distance distributions...")

# rows = barrel residues
# columns = tail residues
barrel_resids = hydrophobic_data.index.astype(int)
tail_resids = hydrophobic_data.columns.astype(int)

hydrophobic_pairs = []

for barrel_resid in barrel_resids:

    for tail_resid in tail_resids:

        occupancy = hydrophobic_data.loc[
            str(barrel_resid),
            str(tail_resid)
        ]

        if occupancy > 0:
            hydrophobic_pairs.append(
                (tail_resid, barrel_resid)
            )

for tail_resid, barrel_resid in hydrophobic_pairs:

    tail_atoms_pair = tail_hphob.select_atoms(
        f"resid {tail_resid}"
    )

    barrel_atoms_pair = barrel_hphob.select_atoms(
        f"resid {barrel_resid}"
    )

    min_distances = []

    for ts in u.trajectory:

        d = distances.distance_array(
            tail_atoms_pair.positions,
            barrel_atoms_pair.positions,
            box=u.dimensions
        )

        min_distances.append(np.min(d))

    min_distances = np.array(min_distances)

    plt.figure(figsize=(5, 5))

    sns.histplot(
        min_distances,
        bins=50,
        stat="density",
        kde=True
    )

    plt.axvline(
        HYDROPHOBIC_CUTOFF,
        color="grey",
        linestyle=":",
        linewidth=2
    )

    plt.xlim(0, 10)

    plt.xlabel("Minimum heavy-atom distance (Å)")
    plt.ylabel("Density")

    plt.title(
        f"{res_lookup[tail_resid]} ↔ "
        f"{res_lookup[barrel_resid]}"
    )

    plt.tight_layout()

    outfile = (
        f"hydrophobic_kdes/"
        f"{res_lookup[tail_resid]}_"
        f"{res_lookup[barrel_resid]}.pdf"
    )

    plt.savefig(outfile)
    plt.close()

# ---------------------------------------------------------
# Per-frame interaction counts
# ---------------------------------------------------------

print("Calculating per-frame interaction counts...")

# =========================================================
# H-bonds per frame
# =========================================================

# Count H-bonds detected in each frame
hbond_counts = (
    hbonds_df.groupby("frame")
    .size()
    .reindex(
        range(len(u.trajectory)),
        fill_value=0
    )
)

hbond_counts_df = pd.DataFrame({
    "frame": hbond_counts.index,
    "n_hbonds": hbond_counts.values
})

hbond_counts_df.to_csv(
    "hbonds_per_frame.csv",
    index=False
)

# Plot distribution
plt.figure(figsize=(5,4))

sns.histplot(
    hbond_counts_df["n_hbonds"],
    bins=np.arange(
        -0.5,
        hbond_counts_df["n_hbonds"].max() + 1.5,
        1
    ),
    stat="probability"
)

plt.xlabel("Number of hydrogen bonds")
plt.ylabel("Probability")

plt.tight_layout()
plt.savefig("hbonds_count_distribution.pdf")
plt.close()


# =========================================================
# Hydrophobic contacts per frame
# =========================================================

hydrophobic_counts = []

for ts in u.trajectory:

    d = distances.distance_array(
        tail_hphob.positions,
        barrel_hphob.positions,
        box=u.dimensions
    )

    contacts = np.where(d < HYDROPHOBIC_CUTOFF)

    unique_pairs = set()

    for i, j in zip(*contacts):

        tail_resid = tail_hphob[i].resid
        barrel_resid = barrel_hphob[j].resid

        unique_pairs.add(
            (tail_resid, barrel_resid)
        )

    hydrophobic_counts.append(len(unique_pairs))

hydrophobic_counts_df = pd.DataFrame({
    "frame": np.arange(len(hydrophobic_counts)),
    "n_hydrophobic_contacts": hydrophobic_counts
})

hydrophobic_counts_df.to_csv(
    "hydrophobic_contacts_per_frame.csv",
    index=False
)

# Plot distribution
plt.figure(figsize=(5,4))

sns.histplot(
    hydrophobic_counts_df["n_hydrophobic_contacts"],
    bins=np.arange(
        -0.5,
        hydrophobic_counts_df[
            "n_hydrophobic_contacts"
        ].max() + 1.5,
        1
    ),
    stat="probability"
)

plt.xlabel("Number of hydrophobic contacts")
plt.ylabel("Probability")

plt.tight_layout()
plt.savefig(
    "hydrophobic_contacts_count_distribution.pdf"
)
plt.close()

print("Done.")