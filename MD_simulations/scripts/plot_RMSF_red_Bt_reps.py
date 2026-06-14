#!/usr/bin/env python3

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap("viridis")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42

files = sys.argv[1:]

if len(files) < 2 or "-h" in files or "--help" in files:
    print("Usage: python RMSF_plot.py rep1.xvg rep2.xvg rep3.xvg")
    sys.exit(1)

# Residues to color red (e.g., DERA tail/active site regions)
residues_to_scale = set().union(
    range(1,11),
    range(26,33),
    range(77,85),
    range(169,180),
    range(198,205),
    range(222,225),
    range(233,239)
)

all_resid = None
all_rmsf = []

# Load all replicas
for f in files:
    resid = []
    rmsf = []

    with open(f) as fh:
        for line in fh:
            if line.startswith(('#', '@')):
                continue
            cols = line.split()
            resid.append(int(cols[0]))
            rmsf.append(float(cols[1]))

    resid = np.array(resid)
    rmsf = np.array(rmsf)

    if all_resid is None:
        all_resid = resid
    else:
        if not np.array_equal(all_resid, resid):
            raise ValueError("Residue numbering mismatch between replicas")

    all_rmsf.append(rmsf)

all_rmsf = np.array(all_rmsf)

# Compute statistics
mean_rmsf = np.mean(all_rmsf, axis=0)
std_rmsf = np.std(all_rmsf, axis=0)

# Plot
plt.figure(figsize=(8, 4))

plt.plot(all_resid, mean_rmsf, linewidth=2, label="mean RMSF")

plt.fill_between(
    all_resid,
    mean_rmsf - std_rmsf,
    mean_rmsf + std_rmsf,
    alpha=0.3,
    label="±1 SD"
)

# Highlight regions in red 
highlight_ranges = [
    (1, 10),
    (26, 32),
    (77, 84),
    (169, 179),
    (198, 204),
    (222, 224),
    (233, 238),
]

for start_res, end_res in highlight_ranges:
    mask = (all_resid >= start_res) & (all_resid <= end_res)
    plt.plot(
        all_resid[mask],
        mean_rmsf[mask],
        color='red',
        linewidth=2
    )

plt.xlabel("Residue number")
plt.ylabel("RMSF (nm)")
plt.legend()

plt.tight_layout()
plt.savefig("rmsf_mean_sd.pdf")
plt.close()

print("Saved: rmsf_mean_sd.pdf")