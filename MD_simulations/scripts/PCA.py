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
from scipy.ndimage import gaussian_filter
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from sklearn.cluster import AgglomerativeClustering
import hdbscan

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
np.savetxt(
    "pca_projection.dat",
    data,
    header="time PC1 PC2",
    fmt=["%d", "%.6f", "%.6f"]
)

plt.figure(figsize=(6, 4))
plt.scatter(pc1, pc2, s=5, color=cmap(0.5), alpha=0.7)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("pc1_vs_pc2.pdf")
plt.close()

### Plot PC1 and PC2 vs distance ###
dist_file = sys.argv[3]
dist_time = []
distance = []

with open(dist_file) as f:
    for line in f:
        if line.startswith("@") or line.startswith("#"):
            continue

        t, d = map(float, line.split())
        dist_time.append(t * 1000)  # convert ns --> ps
        distance.append(d)

dist_time = np.array(dist_time)
distance = np.array(distance)

# match frames by time
common_times = np.intersect1d(time, dist_time)

pc1_match = []
pc2_match = []
dist_match = []

for t in common_times:
    i = np.where(time == t)[0][0]
    j = np.where(dist_time == t)[0][0]

    pc1_match.append(pc1[i])
    pc2_match.append(pc2[i])
    dist_match.append(distance[j])

pc1_match = np.array(pc1_match)
pc2_match = np.array(pc2_match)
dist_match = np.array(dist_match)

corr_pc1 = np.corrcoef(pc1_match, dist_match)[0,1]
corr_pc2 = np.corrcoef(pc2_match, dist_match)[0,1]

### Plot distance vs PC1 ###
plt.figure(figsize=(6,4))
plt.scatter(pc1_match, dist_match, s=5, color=cmap(0.5), alpha=0.7)
plt.xlabel("PC1")
plt.ylabel("Lys167-Tyr259 distance (nm)")
plt.text(
    0.05, 0.95,
    f"r = {corr_pc1:.2f}",
    transform=plt.gca().transAxes,
    verticalalignment='top'
)
plt.tight_layout()
plt.savefig("distance_vs_PC1.pdf")
plt.close()

### Plot distance vs PC2 ###
plt.figure(figsize=(6,4))
plt.scatter(pc2_match, dist_match, s=5, color=cmap(0.5), alpha=0.7)
plt.xlabel("PC2")
plt.ylabel("Lys167-Tyr259 distance (nm)")
plt.text(
    0.05, 0.95,
    f"r = {corr_pc2:.2f}",
    transform=plt.gca().transAxes,
    verticalalignment='top'
)
plt.tight_layout()
plt.savefig("distance_vs_PC2.pdf")
plt.close()

### Plot PC1 vs PC2 colored by distance ###
plt.figure(figsize=(5,4))
plt.scatter(pc1, pc2, c=distance, cmap='viridis', s=5)
plt.colorbar(label="Lys167-Tyr259 distance (nm)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("pc1_pc2_colored_by_distance.pdf")
plt.close()

### Calculate correlation between distance and PCs
proj_20_file = sys.argv[4]

time = []
pcs = []  
current_pc = []

mode = 0  

with open(proj_20_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("@"):
            continue
        if line.startswith("&"):
            if current_pc:
                pcs.append(current_pc)
            current_pc = []
            mode += 1
            continue

        parts = line.split()
        t = float(parts[0])
        val = float(parts[1])

        if mode == 0:
            time.append(t)

        current_pc.append(val)

# Append the last PC block
if current_pc:
    pcs.append(current_pc)

# Convert to numpy arrays
time = np.array(time)
pcs = [np.array(pc) for pc in pcs]

# Stack into a single array: first column = time, then PCs
data = np.column_stack([time] + pcs)
pcs_array = np.column_stack(pcs)  # shape: (n_timepoints, 20)

# Read distance data
dist_time = []
distance = []

with open(dist_file) as f:
    for line in f:
        if line.startswith("@") or line.startswith("#"):
            continue

        t, d = map(float, line.split())

        dist_time.append(t * 1000)  # ns --> ps
        distance.append(d)

dist_time = np.array(dist_time)
distance = np.array(distance)


# Match time points between PCs and distance
common_times = np.intersect1d(time, dist_time)

pc_match = np.array([pcs_array[np.where(time == t)[0][0]] for t in common_times])
dist_match = np.array([distance[np.where(dist_time == t)[0][0]] for t in common_times])

correlations = [np.corrcoef(pc_match[:, i], dist_match)[0,1] for i in range(pcs_array.shape[1])]
correlations = np.array(correlations)


# Save correlations to file
n_pcs = pcs_array.shape[1]
out = np.column_stack((np.arange(1, n_pcs+1), correlations))

np.savetxt(
    "pc_distance_correlations.dat",
    out,
    header="PC  Pearson_r",
    fmt=["%d", "%.5f"]
)

print("Saved correlations to pc_distance_correlations.dat")

### 4 Plot variance explained by PCs ###
eig = np.loadtxt(eigen_file, comments=["@", "#"])
variance = eig[:,1]

variance_percent = variance / np.sum(variance) * 100
variance_percent = variance_percent[:10]  # show only first 10 PCs
pc_index = np.arange(1, len(variance_percent) + 1)

df = pd.DataFrame({
    "PC": pc_index,
    "variance_percent": variance_percent
})
df.to_csv("pc_variance.csv", index=False)

plt.figure(figsize=(8, 4))
plt.bar(pc_index, variance_percent, color=cmap(0.5), edgecolor = "black")
plt.xlabel("PC index")
plt.ylabel("Variance explained (%)")
plt.xticks(pc_index)
plt.tight_layout()
plt.savefig("variance_explained.pdf")
plt.close()

pc1_var = variance_percent[0]
pc2_var = variance_percent[1]

### 2 Free energy landscape ###
kB = 0.008314  # kJ/mol/K
T = 298  # K

H, xedges, yedges = np.histogram2d(pc1, pc2, bins=50)
P = H / np.max(H) # normalize to max bin count --> most populated bin F=0, others F>0
F = -kB * T * np.log(P + 1e-12)  # avoid log(0)
F[P == 0] = np.nan               # hide unsampled regions

plt.figure(figsize=(5, 4))
plt.imshow(F.T, origin='lower',
           extent=[xedges[0], xedges[-1],
                   yedges[0], yedges[-1]],
           aspect='auto')

plt.colorbar(label="Free Energy (kJ/mol)")
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.tight_layout()
plt.savefig("free_energy.pdf")
plt.close()

### Optional: Probability density 
# 2D histogram with density=True
H, xedges, yedges = np.histogram2d(pc1, pc2, bins=50, density=True)
# Apply Gaussian smoothing to the histogram
sigma = 1  # standard deviation in bins, adjust for more/less smoothing
H_smooth = gaussian_filter(H, sigma=sigma)
# probability density -> free energy
F = -kB * T * np.log(H_smooth + 1e-12)  # avoid log(0)
# shift minimum free energy to zero
F = F - np.nanmin(F)
F[P == 0] = np.nan # hide unsampled regions

plt.figure(figsize=(5, 5))
plt.imshow(F.T, origin='lower',
           extent=[xedges[0], xedges[-1],
                   yedges[0], yedges[-1]],
           aspect='auto')

plt.colorbar(label="Free Energy (kJ/mol)")
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.tight_layout()
plt.savefig("free_energy_density.pdf")
plt.close()

### Optional: KDE-based free energy landscape ###
data_kde = np.vstack([pc1, pc2])
kde = gaussian_kde(data_kde, bw_method=0.2)
# Evaluate on a regular grid
x_low, x_high = np.percentile(pc1, [1, 99])
y_low, y_high = np.percentile(pc2, [1, 99])
xgrid = np.linspace(x_low, x_high, 100)
ygrid = np.linspace(y_low, y_high, 100)
X, Y = np.meshgrid(xgrid, ygrid)
grid_points = np.vstack([X.ravel(), Y.ravel()])
P = kde(grid_points).reshape(100, 100) # probability density
F = -kB * T * np.log(P + 1e-12)
F = F - np.nanmin(F) # normalize to min=0
threshold = np.percentile(P, 10)  # keep top 95% density
F[P < threshold] = np.nan

# Contour plot
plt.figure(figsize=(6, 5))
cs = plt.contourf(X, Y, F, levels=20)
plt.colorbar(cs, label="Free energy (kJ/mol)")
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.tight_layout()
plt.savefig("fel_kde.pdf")
plt.close()

# marginals
dx = xgrid[1] - xgrid[0]
dy = ygrid[1] - ygrid[0]
P_x = np.sum(P, axis=0) * dy
P_y = np.sum(P, axis=1) * dx
# Normalize
P_x /= np.max(P_x)
P_y /= np.max(P_y)


# Contour plot (standard for FEL)
fig = plt.figure(figsize=(8, 7))
# Grid layout
gs = fig.add_gridspec(
    2, 2,
    width_ratios=[4, 1],
    height_ratios=[1, 4],
    hspace=0.05,
    wspace=0.05
)
ax_fel = fig.add_subplot(gs[1, 0])
ax_top = fig.add_subplot(gs[0, 0], sharex=ax_fel)
ax_right = fig.add_subplot(gs[1, 1], sharey=ax_fel)

# FEL 
cs = ax_fel.contourf(X, Y, F, levels=20)
cbar = fig.colorbar(cs, ax=ax_fel)
cbar.set_label("Free energy (kJ/mol)")

pc1_var = variance_percent[0]
pc2_var = variance_percent[1]
ax_fel.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
ax_fel.set_ylabel(f"PC2 ({pc2_var:.1f}%)")
ax_fel.set_xlim(-4, 4)
ax_fel.set_ylim(-4, 4)
ax_fel.set_aspect('equal')
ax_fel.tick_params(direction='in')

# Top marginal (PC1)
ax_top.plot(xgrid, P_x)
ax_top.set_ylabel("P(PC1)")
ax_top.tick_params(labelbottom=False)

# Right marginal (PC2)
ax_right.plot(P_y, ygrid)
ax_right.set_xlabel("P(PC2)")
ax_right.tick_params(labelleft=False)

# Clean up spines
ax_top.spines["right"].set_visible(False)
ax_top.spines["top"].set_visible(False)
ax_right.spines["right"].set_visible(False)
ax_right.spines["top"].set_visible(False)

plt.savefig("fel_with_marginals.pdf", bbox_inches="tight")
plt.close()

### Outline plot
fig = plt.figure(figsize=(8, 7))
# Grid layout
gs = fig.add_gridspec(
    2, 2,
    width_ratios=[4, 1],
    height_ratios=[1, 4],
    hspace=0.05,
    wspace=0.05
)
ax_fel = fig.add_subplot(gs[1, 0])
ax_top = fig.add_subplot(gs[0, 0], sharex=ax_fel)
ax_right = fig.add_subplot(gs[1, 1], sharey=ax_fel)

# FEL 
# Use KDE-based FEL for smoother contours
# Define domain from data percentiles to avoid outliers dominating the plot
x_low, x_high = np.percentile(pc1, [1, 99])
y_low, y_high = np.percentile(pc2, [1, 99])
# Padding to avoid cutting off contours at edges
pad_x = 0.1 * (x_high - x_low)
pad_y = 0.1 * (y_high - y_low)
x_low -= pad_x
x_high += pad_x
y_low -= pad_y
y_high += pad_y
# Build grid 
xgrid = np.linspace(x_low, x_high, 100)
ygrid = np.linspace(y_low, y_high, 100)
X, Y = np.meshgrid(xgrid, ygrid)
# Evaluate KDE
grid_points = np.vstack([X.ravel(), Y.ravel()])
P = kde(grid_points).reshape(100, 100) # probability density
# Free energy
F = -kB * T * np.log(P + 1e-12)
F = F - np.nanmin(F) # normalize to min=0
# Masking
threshold = np.percentile(P, 10)  # keep top 95% density
F[P < threshold] = np.nan

cs = ax_fel.contourf(X, Y, F, levels=15, cmap="viridis")
ax_fel.contour(X, Y, F, levels=15, colors="black", linewidths=0.3)
cbar = fig.colorbar(cs, ax=ax_fel, fraction=0.05, pad=0.03)
fig.subplots_adjust(right=0.88)
cbar.set_label("Free energy (kJ/mol)")

pc1_var = variance_percent[0]
pc2_var = variance_percent[1]
ax_fel.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
ax_fel.set_ylabel(f"PC2 ({pc2_var:.1f}%)")
ax_fel.set_xlim(-4, 4)
ax_fel.set_ylim(-4, 4)
ax_fel.set_aspect('equal')
ax_fel.tick_params(direction='in')

# Top marginal (PC1)
ax_top.plot(xgrid, P_x, color="black", alpha=0.1, linewidth=1)
ax_top.fill_between(xgrid, P_x, color="black", alpha=0.1)
ax_top.set_ylabel("P(PC1)")
ax_top.tick_params(labelbottom=False)
ax_top.set_yticks([])
ax_top.set_xticks([])
ax_top.set_xlim(ax_fel.get_xlim())

# Right marginal (PC2)
ax_right.plot(P_y, ygrid, color="black", alpha=0.1, linewidth=1)
ax_right.fill_betweenx(ygrid, P_y, color="black", alpha=0.1)
ax_right.set_xlabel("P(PC2)")
ax_right.tick_params(labelleft=False)
ax_right.set_xticks([])
ax_right.set_yticks([])
ax_right.set_ylim(ax_fel.get_ylim())

# Clean up spines
ax_top.spines["right"].set_visible(False)
ax_top.spines["top"].set_visible(False)
ax_right.spines["right"].set_visible(False)
ax_right.spines["top"].set_visible(False)

plt.savefig("fel_with_marginals_lines.pdf", bbox_inches="tight")
plt.close()

print(pc1.min(), pc1.max())
print(pc2.min(), pc2.max())
print(x_low, x_high)
print(y_low, y_high)

# Cluster 
X = np.vstack([pc1, pc2]).T

for i in range(2, 10):
    clustering = AgglomerativeClustering(n_clusters=5)
    labels = clustering.fit_predict(X)
    # plot pc1 vs pc2 colored by cluster
    plt.figure(figsize=(6, 4))
    plt.scatter(pc1, pc2, s=5, color=cmap(labels), alpha=0.7)
    plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
    plt.tight_layout()
    plt.savefig(f"pc1_vs_pc2_clusters{i}.pdf")
    plt.close()

#hdbscan
clusterer = hdbscan.HDBSCAN(min_cluster_size=100)
labels = clusterer.fit_predict(X)
# plot pc1 vs pc2 colored by cluster
plt.figure(figsize=(6, 4))
plt.scatter(pc1, pc2, s=5, color=cmap(labels), alpha=0.7)
plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
plt.ylabel(f"PC2 ({pc2_var:.1f}%)")
plt.tight_layout()
plt.savefig(f"pc1_vs_pc2_clusters_hdbscan.pdf")
plt.close()