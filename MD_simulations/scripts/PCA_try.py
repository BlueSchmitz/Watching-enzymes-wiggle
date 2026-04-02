#!/usr/bin/env python3
# usage: python PCA.py proj.xvg eigenvalues.xvg distance.xvg proj_20_pcs.xvg
'''Visualisation of the PCA analysis.'''

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from scipy.stats import gaussian_kde

cmap = plt.get_cmap("viridis")

### 1 Plot PC1 vs. PC2 ###
# Load projection file (skip comments)
proj_file = sys.argv[1]
eigen_file = sys.argv[2]

pc1 = []
pc2 = []
time = []

mode = 0

with open(proj_file) as f:
    for line in f:
        if line.startswith("&"):
            mode = 1
            continue
        if line.startswith("@") or line.startswith("#"):
            continue

        t, val = map(float, line.split())

        if mode == 0:
            time.append(t)
            pc1.append(val)
        else:
            pc2.append(val)

time = np.array(time)
pc1 = np.array(pc1)
pc2 = np.array(pc2)
data = np.column_stack((time.astype(int), pc1, pc2))

eig = np.loadtxt(eigen_file, comments=["@", "#"])
variance = eig[:,1]

variance_percent = variance / np.sum(variance) * 100
variance_percent = variance_percent[:10]  # show only first 10 PCs
pc_index = np.arange(1, len(variance_percent) + 1)

### 2 Free energy landscape ###
kB = 0.008314  # kJ/mol/K
T = 298  # K

data_kde = np.vstack([pc1, pc2])
kde = gaussian_kde(data_kde, bw_method=0.2)

### Outline plot
fig = plt.figure(figsize=(6, 6))

# Use KDE-based FEL for smoother contours
# Prepare grid
pc1_min, pc1_max = pc1.min(), pc1.max()
pc2_min, pc2_max = pc2.min(), pc2.max()
# Add padding
pad1 = 0.2 * (pc1_max - pc1_min)
pad2 = 0.2 * (pc2_max - pc2_min)
x_min, x_max = pc1_min - pad1, pc1_max + pad1
y_min, y_max = pc2_min - pad2, pc2_max + pad2
# Build grid 
X, Y = np.meshgrid(
    np.linspace(x_min, x_max, 100),
    np.linspace(y_min, y_max, 100)
)
grid_points = np.vstack([X.ravel(), Y.ravel()])
# Evaluate KDE
P = kde(grid_points).reshape(100, 100) # probability density
# mask unsampled regions (keep top 95% density)
threshold = np.percentile(P, 5)
mask = np.log(P) < np.log(P.max()) - 6
# Free energy
F = -kB * T * np.log(P + 1e-12)
F = F - np.nanmin(F) # normalize to min=0
F_masked = np.ma.array(F, mask=mask)

cmap = plt.cm.viridis.copy()
cmap.set_bad(color='white')  # unsampled regions = white background
cs = plt.contourf(X, Y, F_masked, levels=15, cmap=cmap)
plt.contour(X, Y, F_masked, levels=15, colors="black", linewidths=0.3)
plt.colorbar(cs, label="Free energy (kJ/mol)")

pc1_var = variance_percent[0]
pc2_var = variance_percent[1]
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.tick_params(direction='in')

plt.savefig("fel_try.png", dpi=300, bbox_inches="tight")
plt.close()

### 2 Free energy landscape ###
kB = 0.008314  # kJ/mol/K
T = 298  # K

H, xedges, yedges = np.histogram2d(pc1, pc2, bins=50)
P = H / np.max(H) # normalize to max bin count --> most populated bin F=0, others F>0
F = -kB * T * np.log(P + 1e-12)  # avoid log(0)
F[P == 0] = np.nan               # hide unsampled regions

# Bin centers
xc = 0.5 * (xedges[:-1] + xedges[1:])
yc = 0.5 * (yedges[:-1] + yedges[1:])

# Contour plot
fig, ax = plt.subplots(figsize=(7, 5))

levels = np.linspace(np.nanmin(F), np.nanpercentile(F, 95), 15)

cf = ax.contourf(xc, yc, F.T, levels=levels, cmap="viridis")
c = ax.contour(xc, yc, F.T, levels=levels, colors="k", linewidths=0.4, alpha=0.5)

fig, ax = plt.subplots(figsize=(7, 5))
levels = np.linspace(np.nanmin(F), np.nanpercentile(F, 95), 15)
cf = ax.contourf(xc, yc, F.T, levels=levels, cmap="viridis")
c = ax.contour(xc, yc, F.T, levels=levels, colors="k", linewidths=0.4, alpha=0.5)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label("Free energy (kJ/mol)")
plt.tight_layout()
plt.savefig("free_energy_try.png", dpi=300)
plt.close()