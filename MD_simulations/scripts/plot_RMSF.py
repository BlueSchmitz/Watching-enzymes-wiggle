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

# Plot
plt.figure(figsize=(8, 4))
plt.plot(resid, rmsf, linewidth=2)
plt.xlabel("Residue number")
plt.ylabel("RMSF (nm)")
plt.title("Cα RMSF per residue")
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(out_png, dpi=300)
plt.close()

print(f"Saved RMSF plot to {out_png}")