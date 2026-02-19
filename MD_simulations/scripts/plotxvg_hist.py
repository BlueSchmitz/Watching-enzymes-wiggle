### This script plots histograms from .xvg files ###

# Import libraries
import sys
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# Get xvg file
xvg_filename = sys.argv[1]

# Print help command
if xvg_filename in ("-h", "--help"):
    print('''
          Use the script plot_histogram.py using the syntax below:
          
          python3 plot_histogram.py {filename}.xvg
          ''')
    exit()

# Load data from .xvg file (ignore @ and # lines)
x, y = np.loadtxt(xvg_filename, comments = ["@", "#"], unpack = True)

# Plot histogram of distances
plt.figure(figsize=(6,4))
bin_width = x[1] - x[0]
print("Summed probability density:", np.sum(y * bin_width))
x_centers = x + bin_width / 2
plt.bar(
    x_centers,
    y,
    width=bin_width,
    align="center",
    edgecolor="black",
    linewidth=0.8
)
plt.xlabel("Distance (nm)")
plt.ylabel(r"Probability density (nm$^{-1}$)")
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(0.4)) # Set x-axis ticks every 0.4 nm
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))  # remove unnecessary zeros
plt.xlim(left=0)

# Save figure
plt.tight_layout()
plt.savefig(xvg_filename.split(".")[0] + "_histogram.png")
print(f"Histogram saved as {xvg_filename.split('.')[0]}_histogram.png")
