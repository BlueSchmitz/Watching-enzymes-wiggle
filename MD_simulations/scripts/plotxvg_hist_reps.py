### This script plots histograms from multiple .xvg files on one plot ###

import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import re

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42

cmap = plt.get_cmap("viridis")

# Get input files
xvg_files = sys.argv[1:]

# Help
if "-h" in xvg_files or "--help" in xvg_files:
    print("""
    Use:
    python plot_histogram_reps.py file1.xvg file2.xvg file3.xvg
    """)
    exit()

# Color map for replicas
colors = [cmap(i) for i in np.linspace(0.2, 0.8, len(xvg_files))]

plt.figure(figsize=(5,3))

xlabel_set = False

for i, xvg_filename in enumerate(xvg_files):

    # Load data
    x, y = np.loadtxt(xvg_filename, comments=["@", "#"], unpack=True)

    # Bin width assumption (GROMACS-style histogram)
    bin_width = x[1] - x[0]
    print(f"{xvg_filename} summed probability density:",
          np.sum(y * bin_width))

    x_centers = x + bin_width / 2

    plt.bar(
        x_centers,
        y,
        width=bin_width,
        align="center",
        edgecolor="black",
        linewidth=0.6,
        alpha=0.4,
        color=colors[i],
        label=f"rep{i+1}"
    )

    # Extract xlabel once
    if not xlabel_set:
        with open(xvg_filename, "r") as f:
            for line in f:
                if "xaxis  label" in line:
                    res = re.search('xaxis  label "(.*)"', line)
                    xlabel = res.group(1)
                    xlabel_set = True

# Reference line
plt.axvline(x=0.6, color='grey', linestyle='--', linewidth=1)

# Labels (fixed y-label for histograms)
plt.xlabel(xlabel if xlabel_set else "Distance (nm)")
plt.ylabel(r"Probability density (nm$^{-1}$)")

# Axis formatting
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(0.4))
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))

plt.xlim(0,3.2)
plt.ylim(0, 6)

plt.legend(frameon=False)

plt.tight_layout()

outname = xvg_files[0].split(".")[0] + "_reps_histogram.pdf"
plt.savefig(outname)

print(f"Saved: {outname}")
plt.close()