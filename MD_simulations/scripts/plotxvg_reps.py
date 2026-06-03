### This script plots multiple xvg files on one plot ###

import sys
import matplotlib.pyplot as plt
import numpy as np
import re

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42

# Get input files
xvg_files = sys.argv[1:]

# Help
if "-h" in xvg_files or "--help" in xvg_files:
    print("""
    Use:
    python plotxvg_reps.py file1.xvg file2.xvg file3.xvg
    """)
    exit()

# Set color map for replicas
cmap = plt.get_cmap("viridis")
colors = [cmap(i) for i in np.linspace(0.2, 0.8, len(xvg_files))]

plt.figure(figsize=(5,4))

xlabel = None
ylabel = None

for i, xvg_filename in enumerate(xvg_files):

    # Load data
    if xvg_filename == 'gyrate.xvg':
        r = np.loadtxt(xvg_filename, comments=["@", "#"], unpack=True)
        x = r[0]
        y = r[1]
    else:
        x, y = np.loadtxt(xvg_filename, comments=["@", "#"], unpack=True)

    # Read labels only once (from first file)
    if xlabel is None or ylabel is None:
        with open(xvg_filename, "r") as f:
            for line in f:
                if "xaxis  label" in line:
                    xlabel = re.search('xaxis  label "(.*)"', line).group(1)
                if "yaxis  label" in line:
                    ylabel = re.search('yaxis  label "(.*)"', line).group(1)

    # Plot each replica
    rep_label = f"rep{i+1}"
    plt.scatter(x, y, s=3, alpha=0.5, linewidths=0, rasterized=True, color=colors[i], label=rep_label)

# Reference line (keep your original)
#plt.ylim(0, 1)

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(frameon=False)

plt.tight_layout()

# Output name based on first file
outname = xvg_files[0].split(".")[0] + "_reps_plot.pdf"
plt.savefig(outname)
plt.close()