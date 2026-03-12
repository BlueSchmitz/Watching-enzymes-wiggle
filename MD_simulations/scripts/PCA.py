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
plt.scatter(pc1, pc2, s=5, color = "#287c8e", alpha=0.7)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("pc1_vs_pc2.png", dpi=300)
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
plt.scatter(pc1_match, dist_match, s=5, color = "#287c8e", alpha=0.7)
plt.xlabel("PC1")
plt.ylabel("Lys167-Tyr259 distance (nm)")
plt.text(
    0.05, 0.95,
    f"r = {corr_pc1:.2f}",
    transform=plt.gca().transAxes,
    verticalalignment='top'
)
plt.tight_layout()
plt.savefig("distance_vs_PC1.png", dpi=300)
plt.close()

### Plot distance vs PC2 ###
plt.figure(figsize=(6,4))
plt.scatter(pc2_match, dist_match, s=5, color = "#287c8e", alpha=0.7)
plt.xlabel("PC2")
plt.ylabel("Lys167-Tyr259 distance (nm)")
plt.text(
    0.05, 0.95,
    f"r = {corr_pc2:.2f}",
    transform=plt.gca().transAxes,
    verticalalignment='top'
)
plt.tight_layout()
plt.savefig("distance_vs_PC2.png", dpi=300)
plt.close()

### Plot PC1 vs PC2 colored by distance ###
plt.figure(figsize=(6,4))
plt.scatter(pc1, pc2, c=distance, cmap='viridis', s=5)
plt.colorbar(label="Lys167-Tyr259 distance (nm)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("pc1_pc2_colored_by_distance.png", dpi=300)
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


### 2 Free energy landscape ###
kB = 0.008314  # kJ/mol/K
T = 298  # K

H, xedges, yedges = np.histogram2d(pc1, pc2, bins=50)
P = H / np.max(H) # normalize to max bin count --> most populated bin F=0, others F>0
F = -kB * T * np.log(P + 1e-12)  # avoid log(0)
#F[P == 0] = np.nan               # hide unsampled regions

plt.figure(figsize=(8, 4))
plt.imshow(F.T, origin='lower',
           extent=[xedges[0], xedges[-1],
                   yedges[0], yedges[-1]],
           aspect='auto')

plt.colorbar(label="Free Energy (kJ/mol)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("free_energy.png", dpi=300)
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
plt.figure(figsize=(8, 4))
plt.imshow(F.T, origin='lower',
           extent=[xedges[0], xedges[-1],
                   yedges[0], yedges[-1]],
           aspect='auto')

plt.colorbar(label="Free Energy (kJ/mol)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("free_energy_density.png", dpi=300)
plt.close()

### Optional: KDE-based free energy landscape ###
data_kde = np.vstack([pc1, pc2])
kde = gaussian_kde(data_kde)
# Evaluate on a regular grid
xgrid = np.linspace(min(pc1), max(pc1), 100)
ygrid = np.linspace(min(pc2), max(pc2), 100)
X, Y = np.meshgrid(xgrid, ygrid)
grid_points = np.vstack([X.ravel(), Y.ravel()])

P = kde(grid_points).reshape(100, 100)  # probability density
F = -kB * T * np.log(P + 1e-12)
F = F - np.nanmin(F)  # normalize to min=0
F[P < 1e-10] = np.nan  # mask very low probabilities

# Contour plot (standard for FEL)
plt.figure(figsize=(6, 5))
cs = plt.contourf(X, Y, F, levels=20)
plt.colorbar(cs, label="Free energy (kJ/mol)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("fel_kde.png", dpi=300)
plt.close()

### 3 Extract timepoints of representative structures ###
min_pc1 = np.argmin(pc1)
max_pc1 = np.argmax(pc1)
min_pc2 = np.argmin(pc2)
max_pc2 = np.argmax(pc2)

# corresponding times
min_pc1_time = int(time[min_pc1])
max_pc1_time = int(time[max_pc1])
min_pc2_time = int(time[min_pc2])
max_pc2_time = int(time[max_pc2])

with open("pc_extreme_frames.dat", "w") as f:
    f.write(f"min_pc1 {min_pc1_time}\n")
    f.write(f"max_pc1 {max_pc1_time}\n")
    f.write(f"min_pc2 {min_pc2_time}\n")
    f.write(f"max_pc2 {max_pc2_time}\n")

### 4 Plot variance explained by PCs ###
eig = np.loadtxt(eigen_file, comments=["@", "#"])
variance = eig[:,1]

variance_percent = variance / np.sum(variance) * 100
variance_percent = variance_percent[:10]  # show only first 10 PCs
pc_index = np.arange(1, len(variance_percent) + 1)

plt.figure(figsize=(8, 4))
plt.bar(pc_index, variance_percent, color = "#287c8e", edgecolor = "black")
plt.xlabel("PC index")
plt.ylabel("Variance explained (%)")
plt.xticks(pc_index)
plt.tight_layout()
plt.savefig("variance_explained.png", dpi=300)
plt.close()