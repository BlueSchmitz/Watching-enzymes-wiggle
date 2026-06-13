#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as colors

def truncate_cmap(cmap_name, minval=0.2, maxval=1.0, n=256):
    cmap = plt.get_cmap(cmap_name)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        f"trunc({cmap_name},{minval},{maxval})",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

blues_dark = truncate_cmap("Blues", 0.3, 1.0)
oranges_dark = truncate_cmap("Oranges", 0.3, 1.0)

# =========================
# LOAD H-BOND DATA
# =========================
pair_counts = pd.read_csv("hbonds_occupancy.csv")

hbonds = pair_counts.pivot(
    index="barrel_resid",
    columns="tail_resid",
    values="occupancy"
).fillna(0)

bg = np.ones_like(hbonds)  # same shape, all ones

# =========================
# LOAD HYDROPHOBIC DATA
# =========================
hydrophobic = pd.read_csv(
    "hydrophobic_contacts_matrix.csv",
    index_col=0
)

hydrophobic.index = hydrophobic.index.astype(int)
hydrophobic.columns = hydrophobic.columns.astype(int)

# =========================
# ALIGN MATRICES (FIXED)
# =========================

# keep only barrel residues that have ANY interaction in either dataset
active_barrel = sorted(
    set(hbonds.index[hbonds.sum(axis=1) > 0]) |
    set(hydrophobic.index[hydrophobic.sum(axis=1) > 0])
)

# same logic for tail residues
active_tail = sorted(
    set(hbonds.columns[hbonds.sum(axis=0) > 0]) |
    set(hydrophobic.columns[hydrophobic.sum(axis=0) > 0])
)

hbonds = hbonds.reindex(index=active_barrel, columns=active_tail).fillna(0)
hydrophobic = hydrophobic.reindex(index=active_barrel, columns=active_tail).fillna(0)

# =========================
# PLOT
# =========================
fig, ax = plt.subplots(figsize=(9, 7))

# colorbars
cbar_ax_hphob = fig.add_axes([0.92, 0.55, 0.02, 0.30])
cbar_ax_hbond = fig.add_axes([0.92, 0.15, 0.02, 0.30])
hb_mappable = cm.ScalarMappable(norm=colors.Normalize(0,1), cmap=blues_dark)
hp_mappable = cm.ScalarMappable(norm=colors.Normalize(0,1), cmap=oranges_dark)

# --- GRID (background) ---
sns.heatmap(
    np.ones_like(hbonds),
    cmap=["#f7f7f7"],
    cbar=False,
    linewidths=0.3,
    linecolor="#d0d0d0",
    ax=ax
)

# --- HYDROPHOBIC (orange) ---
sns.heatmap(
    hydrophobic,
    cmap=oranges_dark,
    vmin=0,
    vmax=1,
    mask=(hydrophobic == 0),
    linewidths=0,
    cbar=False,
    ax=ax
)

# --- H-BONDS (blue) ---
sns.heatmap(
    hbonds,
    cmap=blues_dark,
    vmin=0,
    vmax=1,
    mask=(hbonds == 0),
    linewidths=0,
    cbar=False,
    ax=ax
)

cbar_ax_hbond.set_title("H-bonds", fontsize=9)
cbar_ax_hphob.set_title("Hydrophobic", fontsize=9)
fig.colorbar(hb_mappable, cax=cbar_ax_hbond)
fig.colorbar(hp_mappable, cax=cbar_ax_hphob)

ax.set_xlabel("C-terminal tail residues")
ax.set_ylabel("TIM barrel residues")

plt.tight_layout(rect=[0, 0, 0.9, 1])

plt.savefig(
    "combined_interaction_heatmap.pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.close()