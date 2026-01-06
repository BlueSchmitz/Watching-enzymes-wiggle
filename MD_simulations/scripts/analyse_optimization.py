#!/usr/bin/env python3
# usage: python analyse_optimization.py path/to/tune_summary.csv

import sys
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
if "HREX" in file_path:
    output_csv = "tune_summary_HREX_sorted.csv"
else:
    output_csv = "tune_summary_sorted.csv"
df_sorted.to_csv(output_csv, index=False)
print(f"\nSorted table saved to: {output_csv}")

# Scatter plot
plt.figure(figsize=(10, 6))
if "HREX" in file_path:
    sc = plt.scatter(
        df_sorted["Cores per replica"],
        df_sorted["Performance (ns/day)"],
        c=df_sorted["MPI ranks per replica (np)"],  # color by ranks per replica
        s=df_sorted["MPI ranks per replica (np)"] * 100,  # size proportional to ranks per replica
        cmap='viridis',
        alpha=0.8,
        edgecolors='k'
    )
else:
    sc = plt.scatter(
        df_sorted["Total cores"],
        df_sorted["Performance (ns/day)"],
        c=df_sorted["MPI ranks (np)"],  # color by ranks 
        s=df_sorted["MPI ranks (np)"] * 100,  # size proportional to ranks 
        cmap='viridis',
        alpha=0.8,
        edgecolors='k'
    )    

# Colorbar
cbar = plt.colorbar(sc)
if "HREX" in file_path:
    cbar.set_label("MPI ranks per replica (np)")
else:
    cbar.set_label("MPI ranks (np)")

# Labels and title
plt.xlabel("Total cores")
plt.ylabel("Performance (ns/day)")
if "HREX" in file_path:
    plt.title("GROMACS HREX Tuning Performance for HREX with 4 Replicas")
else:
    plt.title("GROMACS Tuning Performance")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# Save plot to PNG
if "HREX" in file_path:
    output_plot = "tune_summary_HREX_scatter.png"
else:
    output_plot = "tune_summary_scatter.png"
plt.savefig(output_plot, dpi=300)
print(f"Scatter plot saved to: {output_plot}")