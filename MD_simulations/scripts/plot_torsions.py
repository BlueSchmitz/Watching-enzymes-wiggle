#!/usr/bin/env python3
# usage: python plot_torsions.py 
'''Plot the dihedrals of the tail residues to check for the minimal scaling factor needed.'''

import pandas as pd
import matplotlib.pyplot as plt
import math

lambdas = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
residues = range(250, 260)

# Read all COLVARs
data = {}  
for l in lambdas:
    path = f"./l{l}/COLVAR"

    # Extract column names from the header
    with open(path) as f:
        for line in f:
            if line.startswith("#! FIELDS"):
                columns = line.split()[2:]  # skip "#!" and "FIELDS"
                break

    # Load dataframe with correct column names
    df = pd.read_csv(path, delim_whitespace=True, comment="#", names=columns)
    data[l] = df

# Plot one figure per lambda
for l in lambdas:
    df = data[l]

    n_res = len(residues)
    ncols = 5  # adjust for aesthetics
    nrows = math.ceil(n_res / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3*ncols, 2.5*nrows),
        sharex=True,
        sharey=True
    )

    axes = axes.flatten()

    for i, r in enumerate(residues):
        ax = axes[i]

        phi = f"phi{r}"
        psi = f"psi{r}"

        if phi in df.columns:
            ax.plot(df["time"], df[phi], label="φ", linewidth=0.8)
        if psi in df.columns:
            ax.plot(df["time"], df[psi], label="ψ", linewidth=0.8)

        ax.set_title(f"Residue {r}", fontsize=9)
        ax.set_ylim(-3.2, 3.2)

    # Remove unused subplots (if any)
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    # One global legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")

    fig.suptitle(f"Torsions for λ = {l}", fontsize=14)
    fig.supxlabel("Time")
    fig.supylabel("Angle [rad]")

    plt.tight_layout()
    plt.savefig(f"torsions_lambda_{l}.png", dpi=300)
    plt.close(fig)

    print(f"Saved torsions_lambda_{l}.png")