#!/usr/bin/env python3
# usage: python PCA.py proj.xvg eigenvalues.xvg
'''Visualisation of the PCA analysis.'''

import sys
import numpy as np
import matplotlib.pyplot as plt

### 1 Plot PC1 vs. PC2 ###
# Load projection file (skip comments)
proj_file = sys.argv[1]
eigen_file = sys.argv[2]
data = np.loadtxt(proj_file, comments=["@", "#"])

pc1 = data[:, 1]
pc2 = data[:, 2]

plt.figure(figsize=(8, 4))
plt.scatter(pc1, pc2, s=5, alpha=0.5)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("pc1_vs_pc2.png", dpi=300)
plt.close()

### 2 Free energy landscape ###
kB = 0.008314  # kJ/mol/K
T = 300

H, xedges, yedges = np.histogram2d(pc1, pc2, bins=100)
P = H / np.sum(H)
F = -kB * T * np.log(P + 1e-12)  # avoid log(0)

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

### 3 Extract representative structures ###
min_pc1 = np.argmin(pc1)
max_pc1 = np.argmax(pc1)
min_pc2 = np.argmin(pc2)
max_pc2 = np.argmax(pc2)

print("Frame with minimum PC1:", min_pc1)
print("Frame with maximum PC1:", max_pc1)
print("Frame with minimum PC2:", min_pc2)
print("Frame with maximum PC2:", max_pc2)

### 4 Plot variance explained by PCs ###
eig = np.loadtxt(eigen_file, comments=["@", "#"])
variance = eig[:,1]

plt.figure(figsize=(8, 4))
plt.plot(variance / np.sum(variance) * 100)
plt.xlabel("PC index")
plt.ylabel("Variance explained (%)")
plt.tight_layout()
plt.savefig("variance_explained.png", dpi=300)
plt.close()