#!/usr/bin/env python3
# usage: python PCA.py proj.xvg eigenvalues.xvg
'''Visualisation of the PCA analysis.'''

import sys
import numpy as np
import matplotlib.pyplot as plt
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

plt.figure(figsize=(8, 4))
plt.scatter(pc1, pc2, s=5, color = "#287c8e", alpha=0.5)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("pc1_vs_pc2.png", dpi=300)
plt.close()

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