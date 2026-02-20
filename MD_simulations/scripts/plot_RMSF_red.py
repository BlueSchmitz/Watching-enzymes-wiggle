#!/usr/bin/env python3
# usage: RMSF_plot.py path/to/rmsf.xvg
'''Plot the RMSF of a trajectory.'''

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Handle command-line input
if len(sys.argv) < 2:
    print("Usage: python plot_RMSF.py <rmsf.xvg>.")
    sys.exit(1)

# Input / output
rmsf_file = Path(sys.argv[1])
out_png = "rmsf.png"

# Residues to color red (e.g., DERA tail/active site regions)
residues_to_scale = set().union(
    range(19,26),
    range(75,83),
    range(168,179),
    range(202,209),
    range(249,260)
)

# Load RMSF data (skip GROMACS comments)
resid = []
rmsf = []

with open(rmsf_file) as f:
    for line in f:
        if line.startswith(('#', '@')):
            continue
        cols = line.split()
        resid.append(int(cols[0]))
        rmsf.append(float(cols[1]))

resid = np.array(resid)
rmsf = np.array(rmsf)

# Plot with continuous blue backbone + red highlights
plt.figure(figsize=(8, 4))

# First: plot entire line in blue (provides connections)
plt.plot(resid, rmsf, color='tab:blue', linewidth=2)

# Then: overplot ONLY the specified red residues
red_mask = np.isin(resid, list(residues_to_scale))
red_indices = np.where(red_mask)[0]

if len(red_indices) > 0:
    # Find continuous red segments
    start = red_indices[0]

    for i in range(1, len(red_indices)):
        if red_indices[i] != red_indices[i - 1] + 1:
            end = red_indices[i - 1] + 1
            plt.plot(resid[start:end], rmsf[start:end],
                     color='red', linewidth=2)
            start = red_indices[i]

    # Final segment
    end = red_indices[-1] + 1
    plt.plot(resid[start:end], rmsf[start:end],
             color='red', linewidth=2)

plt.xlabel("Residue number")
plt.ylabel("RMSF (nm)")
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(out_png, dpi=300)
plt.close()

print(f"Saved RMSF plot to {out_png}")
