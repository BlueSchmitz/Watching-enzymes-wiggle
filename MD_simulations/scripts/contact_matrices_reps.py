#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

import seaborn as sns
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances

try:
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
except ImportError:
    from MDAnalysis.analysis.hbonds import HydrogenBondAnalysis as HBA


# =========================
# INPUT
# =========================
tpr = sys.argv[1]
xtc_files = sys.argv[2:]

n_reps = len(xtc_files)
print(f"Running analysis on {n_reps} replicas")


# =========================
# SYSTEM SETUP (from first replica)
# =========================
u0 = mda.Universe(tpr, xtc_files[0])

res_lookup = {r.resid: f"{r.resname}{r.resid}" for r in u0.residues}

tail_sel = "resid 249:259"
barrel_sel = "protein and not resid 249:259"

tail = u0.select_atoms(tail_sel)
barrel = u0.select_atoms(barrel_sel)

tail_res = sorted(tail.residues.resids)
barrel_res = sorted(barrel.residues.resids)

hydrophobic_sel = "resname ALA VAL LEU ILE MET PHE TRP TYR PRO and not name H*"


# =========================
# STORAGE FOR DISTRIBUTIONS
# =========================
hb_dist_all = []
hydro_dist_all = []

# STORAGE FOR HEATMAPS
hb_maps = []
hydro_maps = []


# =========================
# FUNCTIONS
# =========================
def compute_hbond_matrix(u):

    h = HBA(
        u,
        between=[[tail_sel, barrel_sel]],
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=120
    )
    h.run()

    hbonds = h.results.hbonds
    donor_ix = hbonds[:,1].astype(int)
    acceptor_ix = hbonds[:,3].astype(int)

    df = pd.DataFrame({
        "frame": hbonds[:,0].astype(int),
        "donor_resid": u.atoms[donor_ix].resids,
        "acceptor_resid": u.atoms[acceptor_ix].resids,
    })

    df["barrel_resid"] = np.where(
        df["donor_resid"].isin(tail_res),
        df["acceptor_resid"],
        df["donor_resid"]
    )

    df["tail_resid"] = np.where(
        df["donor_resid"].isin(tail_res),
        df["donor_resid"],
        df["acceptor_resid"]
    )

    pair_counts = (
        df.groupby(["tail_resid", "barrel_resid"])
          .agg(n_frames=("frame", "nunique"))
          .reset_index()
    )

    n_frames = len(u.trajectory)
    pair_counts["occupancy"] = pair_counts["n_frames"] / n_frames

    mat = pair_counts.pivot(
        index="barrel_resid",
        columns="tail_resid",
        values="occupancy"
    ).reindex(index=barrel_res, columns=tail_res).fillna(0)

    return mat


def compute_hydrophobic_matrix(u):

    tail_atoms = u.select_atoms(tail_sel).select_atoms(hydrophobic_sel)
    barrel_atoms = u.select_atoms(barrel_sel).select_atoms(hydrophobic_sel)

    contact_counts = {}
    hydro_per_frame = []

    for ts in u.trajectory:

        d = distances.distance_array(
            tail_atoms.positions,
            barrel_atoms.positions,
            box=u.dimensions
        )

        hydro_per_frame.append(np.sum(d < 4.5))

        contacts = np.where(d < 4.5)
        seen = set()

        for i, j in zip(*contacts):

            pair = (tail_atoms[i].resid, barrel_atoms[j].resid)
            contact_counts.setdefault(pair, set()).add(ts.frame)

            if pair not in seen:
                seen.add(pair)

    rows = []
    n_frames = len(u.trajectory)

    for (t, b), frames in contact_counts.items():
        rows.append({
            "tail_resid": t,
            "barrel_resid": b,
            "occupancy": len(frames) / n_frames
        })

    df = pd.DataFrame(rows)

    mat = df.pivot(
        index="barrel_resid",
        columns="tail_resid",
        values="occupancy"
    ).reindex(index=barrel_res, columns=tail_res).fillna(0)

    return mat, np.array(hydro_per_frame)


# =========================
# LOOP OVER REPLICAS
# =========================
for i, trj in enumerate(xtc_files):

    print(f"Processing replica {i+1}: {trj}")

    u = mda.Universe(tpr, trj)

    # ---- H-bonds ----
    h = HBA(
        u,
        between=[[tail_sel, barrel_sel]],
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=120
    )
    h.run()

    hb_dist_all.append(h.count_by_time())

    # ---- Hydrophobic ----
    hydro_mat, hydro_counts = compute_hydrophobic_matrix(u)

    hydro_dist_all.append(hydro_counts)

    # ---- Heatmaps ----
    hb_maps.append(compute_hbond_matrix(u))
    hydro_maps.append(hydro_mat)


# =========================
# DISTRIBUTIONS (OVERLAID HISTOGRAMS)
# =========================

# ---- H-bonds ----
plt.figure(figsize=(5,4))

for i, data in enumerate(hb_dist_all):
    plt.hist(data, bins=30, density=True, alpha=0.4, label=f"rep{i+1}")

plt.xlabel("Number of H-bonds")
plt.ylabel("Probability density")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("hbonds_distribution.pdf")
plt.close()


# ---- Hydrophobic ----
plt.figure(figsize=(5,4))

for i, data in enumerate(hydro_dist_all):
    plt.hist(data, bins=30, density=True, alpha=0.4, label=f"rep{i+1}")

plt.xlabel("Number of hydrophobic contacts")
plt.ylabel("Probability density")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("hydrophobic_distribution.pdf")
plt.close()


# =========================
# MEAN HEATMAPS
# =========================

hb_mean = sum(hb_maps) / len(hb_maps)
hydro_mean = sum(hydro_maps) / len(hydro_maps)


plt.figure(figsize=(8,6))
sns.heatmap(hb_mean, cmap="viridis", vmin=0, vmax=1, linewidths=0.2)
plt.title("Mean H-bond occupancy (replicas)")
plt.tight_layout()
plt.savefig("hbonds_mean_heatmap.pdf")
plt.close()


plt.figure(figsize=(8,6))
sns.heatmap(hydro_mean, cmap="viridis", vmin=0, vmax=1, linewidths=0.2)
plt.title("Mean hydrophobic occupancy (replicas)")
plt.tight_layout()
plt.savefig("hydrophobic_mean_heatmap.pdf")
plt.close()


print("\nDONE")
print("Outputs:")
print("- hbonds_distribution.pdf")
print("- hydrophobic_distribution.pdf")
print("- hbonds_mean_heatmap.pdf")
print("- hydrophobic_mean_heatmap.pdf")