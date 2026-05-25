#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda

# Load topology only (fast)
u = mda.Universe("../rep1.00/topol.tpr")

# Residue label lookup
res_lookup = {r.resid: f"{r.resname}{r.resid}" for r in u.residues}

### HBOND HEATMAP ###
pair_counts = pd.read_csv("hbonds_occupancy.csv")

tail_res_to_plot = sorted(pair_counts["tail_resid"].unique())
barrel_res_to_plot = sorted(pair_counts["barrel_resid"].unique())

heatmap_data = pair_counts.pivot(
    index="barrel_resid",
    columns="tail_resid",
    values="occupancy"
).reindex(
    index=barrel_res_to_plot,
    columns=tail_res_to_plot
).fillna(0)

# Rename labels
heatmap_data.index = [res_lookup[i] for i in heatmap_data.index]
heatmap_data.columns = [res_lookup[i] for i in heatmap_data.columns]

plt.figure(figsize=(8,6))

mask = heatmap_data == 0

sns.heatmap(
    heatmap_data,
    cmap="viridis",
    vmin=0,
    vmax=1,
    mask=mask,
    linewidths=0.2
)

plt.xlabel("C-terminal tail residues")
plt.ylabel("TIM barrel residues")
plt.tight_layout()
plt.savefig("hbonds_heatmap_replot.pdf")
plt.close()


### HYDROPHOBIC CONTACT HEATMAP ###
hydrophobic_data = pd.read_csv(
    "hydrophobic_contacts_matrix.csv",
    index_col=0
)

# Convert indices back to integers
hydrophobic_data.index = hydrophobic_data.index.astype(int)
hydrophobic_data.columns = hydrophobic_data.columns.astype(int)

# Keep only rows with contacts
barrel_res_with_contacts = hydrophobic_data.index[
    hydrophobic_data.sum(axis=1) > 0
]

hydrophobic_data_plot = hydrophobic_data.loc[
    sorted(barrel_res_with_contacts),
    sorted(hydrophobic_data.columns)
]

# Rename labels
hydrophobic_data_plot.index = [
    res_lookup[i] for i in hydrophobic_data_plot.index
]

hydrophobic_data_plot.columns = [
    res_lookup[i] for i in hydrophobic_data_plot.columns
]

plt.figure(figsize=(8,6))

mask = hydrophobic_data_plot == 0

sns.heatmap(
    hydrophobic_data_plot,
    cmap="viridis",
    vmin=0,
    vmax=1,
    mask=mask,
    linewidths=0.2
)

plt.xlabel("C-terminal tail residues")
plt.ylabel("TIM barrel residues")
plt.tight_layout()
plt.savefig("hydrophobic_contacts_replot.pdf")
plt.close()