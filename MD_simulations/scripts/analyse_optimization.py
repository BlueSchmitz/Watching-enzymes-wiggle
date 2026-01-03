#!/usr/bin/env python3
# usage: python analyse_optimization.py path/to/tune_summary.csv

import sys
import os
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from tabulate import tabulate
matplotlib.use("Agg")  # Use non-GUI backend


# Load the CSV
if len(sys.argv) != 2:
    print("Usage: python analyse_optimization.py path/to/tune_summary.csv")
    sys.exit(1)
file_path = sys.argv[1]
df = pd.read_csv(file_path)

# Sort by Performance (ns/day) descending
df_sorted = df.sort_values("Performance (ns/day)", ascending=False)

# Round numeric columns for readability
df_sorted = df_sorted.round({
    "Performance (ns/day)": 2,
    "Performance (h/ns)": 2,
    "Speed per core (ns/day/core)": 2
})

# Print as a nice table
print(tabulate(df_sorted, headers='keys', tablefmt='fancy_grid', showindex=False))

# Save sorted table to CSV
output_csv = "tune_summary_sorted.csv"
df_sorted.to_csv(output_csv, index=False)
print(f"\nSorted table saved to: {output_csv}")

# Scatter plot
plt.figure(figsize=(10, 6))
sc = plt.scatter(
    df_sorted["Total cores"],
    df_sorted["Performance (ns/day)"],
    c=df_sorted["Speed per core (ns/day/core)"],  # color by speed per core
    s=df_sorted["Speed per core (ns/day/core)"] * 100,  # size proportional to speed per core
    cmap='viridis',
    alpha=0.8,
    edgecolors='k'
)

# Colorbar
cbar = plt.colorbar(sc)
cbar.set_label("Speed per core (ns/day/core)")

# Labels and title
plt.xlabel("Total cores")
plt.ylabel("Performance (ns/day)")
plt.title("GROMACS Tuning Performance")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# Save plot to PNG
output_plot = "tune_summary_scatter.png"
plt.savefig(output_plot, dpi=300)
print(f"Scatter plot saved to: {output_plot}")