### This script plots simple x,y coordinates from .xvg files ###

#Import libraries.
import sys
import matplotlib.pyplot as plt
import numpy as np
import re

# Get xvg file.
xvg_filename = sys.argv[1]
cmap = plt.get_cmap("viridis")

# Print help command.
if xvg_filename == "-h" or xvg_filename == "--help":
    print('''
          Use the script plotxvg.py using the syntax below:
          
          python3 plotxvg.py {filename}.xvg"
          ''')
    exit()


# Get the x and y coordinates from the xvg file.

if xvg_filename == 'gyrate.xvg':
    r = np.loadtxt(xvg_filename, comments = ["@", "#"], unpack = True)
    x = r[0]
    y = r[1]
else :
    x, y = np.loadtxt(xvg_filename, comments = ["@", "#"], unpack = True)


# Open xvg file and loop through through lines.
f = open(xvg_filename, "r")
for line in f:
    # Get x axis label.
    if "xaxis  label " in line:
        res_x = re.search('xaxis  label "(.*)"', line)
        print(res_x.group(1))
    # Get y axis label.
    if "yaxis  label " in line:
        res_y = re.search('yaxis  label "(.*)"', line)
        print(res_y.group(1))
# Close file.
f.close()

# Plot x and y with labels.
plt.scatter(x, y, s=3, alpha=0.5, linewidths=0)
plt.xlabel(res_x.group(1))
plt.ylabel(res_y.group(1))
plt.tight_layout()
plt.savefig(xvg_filename.split(".")[0] + "_" + 'plot.png')
plt.close()