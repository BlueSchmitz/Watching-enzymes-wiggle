#!/usr/bin/env python3
# usage: python contact_matrices.py topol.tpr traj.xtc
'''Visualisation of hbonds, hydrophobic contacts and salt bridges.'''

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # non-interactive
import seaborn as sns
import pandas as pd
import MDAnalysis as mda
try: # MDAnalysis v2.0+ has moved the hbond analysis module
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
except ImportError:
    from MDAnalysis.analysis.hbonds import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis import distances

# Load closed trajectory
tpr = sys.argv[1]
trj = sys.argv[2]
u = mda.Universe(tpr, trj)

res_lookup = {r.resid: "{}{}".format(r.resname, r.resid) for r in u.residues}

# Define tail and barrel selections
tail = u.select_atoms("resid 249:259")
barrel = u.select_atoms("protein and not resid 249:259")
n_tail = len(tail.residues)
n_barrel = len(barrel.residues)

### hbond analysis ###
h = HBA(
    u,
    between=[["resid 249:259", "protein and not resid 249:259"]],  # expects strings
    #update_selections=False, # speed up by not re-evaluating selections every frame
    d_a_cutoff=3.5,
    d_h_a_angle_cutoff=120
)

h.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
)

hbonds = h.results.hbonds
donor_ix = hbonds[:,1].astype(int)
hydrogen_ix = hbonds[:,2].astype(int)
acceptor_ix = hbonds[:,3].astype(int)

# Make dataframe for storing results 
df = pd.DataFrame({
    "frame": hbonds[:,0].astype(int),
    "donor": u.atoms[donor_ix].resnames + u.atoms[donor_ix].resids.astype(str),
    "donor_resid": u.atoms[donor_ix].resids,
    "donor_atom": u.atoms[donor_ix].names,
    "donor_idx": donor_ix,
    "hydrogen_idx": hydrogen_ix,
    "acceptor": u.atoms[acceptor_ix].resnames + u.atoms[acceptor_ix].resids.astype(str),
    "acceptor_resid": u.atoms[acceptor_ix].resids,
    "acceptor_atom": u.atoms[acceptor_ix].names,
    "acceptor_idx": acceptor_ix,
    "DA_distance": hbonds[:,4],
    "DHA_angle": hbonds[:,5],
})
df["pair"] = df["donor"] + "-" + df["acceptor"]
df.to_csv("hbonds.csv", index=False)


### Prepare heatmap ###
# Define tail and barrel residue lists for mapping
tail_res = sorted(tail.residues.resids)
barrel_res = sorted(barrel.residues.resids)

# force tail to be columns and barrel to be rows
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

# count frames per hbond pair
pair_counts = (
    df.groupby(["tail_resid", "barrel_resid"])
      .agg(n_frames=("frame", "nunique"))
      .reset_index()
)
n_total_frames = len(u.trajectory)
pair_counts["occupancy"] = pair_counts["n_frames"] / n_total_frames

# Save occupancy data to CSV
pair_counts.to_csv("hbonds_occupancy.csv", index=False)

# Huge matrix of all tail residues and barrel residues with hbonds
heatmap_data = pd.DataFrame(
    0.0,
    index=tail_res,
    columns=barrel_res
)

# keep all tail residues
tail_res_to_plot = sorted(tail.residues.resids)

# only barrel residues that were involved in hbonds
barrel_res_to_plot = sorted(pair_counts['barrel_resid'].unique())

# pivot your data for the heatmap
heatmap_data = pair_counts.pivot(
    index="barrel_resid",  # barrel on rows → y-axis
    columns="tail_resid",  # tail on columns → x-axis
    values="occupancy"
).reindex(
    index=barrel_res_to_plot,
    columns=tail_res_to_plot
).fillna(0)

# Map resid to resname+resid for better labels
heatmap_data.index = [res_lookup[i] for i in heatmap_data.index]
heatmap_data.columns = [res_lookup[i] for i in heatmap_data.columns]

# Plot heatmap of hbond occupancy
plt.figure(figsize=(8,3))

#mask = heatmap_data == 0  # True where occupancy = 0
sns.heatmap(
    heatmap_data,
    cmap="viridis",
    vmin=0,
    vmax=1,
    #mask=mask,
    linewidths=0.2
)

plt.xlabel("C-terminal tail residues")
plt.ylabel("TIM barrel residues")
plt.tight_layout()
plt.savefig("hbonds_heatmap.png", dpi=300)
plt.close()

# Plot hydrogen bonds over time 
plt.figure(figsize=(8, 4))
plt.plot(h.times, h.count_by_time(), lw=2)
plt.xlabel("Time (ps)")
plt.ylabel("Number of hydrogen bonds")
plt.tight_layout()
plt.savefig("hbonds_over_time.png", dpi=300)
plt.close()

### Hydrophobic contacts and salt bridges ###
# hydrophobic residues: ALA VAL LEU ILE MET PHE TRP TYR PRO
hydrophobic_sel = "resname ALA VAL LEU ILE MET PHE TRP TYR PRO and not name H*"
tail_hphob = tail.select_atoms(hydrophobic_sel)
barrel_hphob = barrel.select_atoms(hydrophobic_sel)

def hydrophobic_matrix(u, tail_atoms, barrel_atoms, tail_res, barrel_res, cutoff=4.5):

    contact_counts = {}

    n_frames = len(u.trajectory)

    for ts in u.trajectory:

        d = distances.distance_array(
            tail_atoms.positions,
            barrel_atoms.positions,
            box=u.dimensions
        )

        contacts = np.where(d < cutoff)

        for i, j in zip(*contacts):

            tail_resid = tail_atoms[i].resid
            barrel_resid = barrel_atoms[j].resid

            pair = (tail_resid, barrel_resid)

            contact_counts.setdefault(pair, set()).add(ts.frame)

    rows = []

    for (tail_resid, barrel_resid), frames in contact_counts.items():

        rows.append({
            "tail_resid": tail_resid,
            "barrel_resid": barrel_resid,
            "occupancy": len(frames) / n_frames
        })

    df = pd.DataFrame(rows)

    matrix = df.pivot(
        index="barrel_resid",   # rows → barrel residues
        columns="tail_resid",   # columns → tail residues
        values="occupancy"
    ).reindex(
        index=barrel_res,  # barrel rows
        columns=tail_res   # tail columns
    ).fillna(0)

    return matrix

hydrophobic_data = hydrophobic_matrix(
    u,
    tail_hphob,
    barrel_hphob,
    tail_res,
    barrel_res
)
hydrophobic_data.to_csv("hydrophobic_contacts_matrix.csv")

# Only barrel residues that actually have contacts
barrel_res_with_contacts = hydrophobic_data.index[hydrophobic_data.sum(axis=1) > 0]
barrel_res_to_plot = sorted(barrel_res_with_contacts)
tail_res_to_plot = sorted(hydrophobic_data.columns)  # tail residues

# Slice the matrix for plotting
hydrophobic_data_plot = hydrophobic_data.loc[barrel_res_to_plot, tail_res_to_plot]

# Map resid → resname+resid
hydrophobic_data_plot.index = [res_lookup[i] for i in hydrophobic_data_plot.index]
hydrophobic_data_plot.columns = [res_lookup[i] for i in hydrophobic_data_plot.columns]

# Plot heatmap of hydrophobic contact occupancy
plt.figure(figsize=(8,3))

#mask = hydrophobic_data_plot == 0
sns.heatmap(
    hydrophobic_data_plot,
    cmap="viridis",
    vmin=0,
    vmax=1,
    #mask=mask,
    linewidths=0.2
)

plt.xlabel("C-terminal tail residues")
plt.ylabel("TIM barrel residues")
plt.tight_layout()
plt.savefig("hydrophobic_contacts.png", dpi=300)
plt.close()
