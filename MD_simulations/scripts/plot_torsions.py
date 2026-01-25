#!/usr/bin/env python3
# usage: python plot_torsions.py 
'''Plot the dihedrals of the tail residues to check for the minimal scaling factor needed.'''

import pandas as pd
import matplotlib.pyplot as plt

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

# Plot
fig, axes = plt.subplots(len(lambdas), 2, figsize=(12, 4*len(lambdas)), sharex=True, sharey='row')

for i, l in enumerate(lambdas):
    df = data[l]

    # φ plot
    for r in residues:
        col_name = f'phi{r}'
        if col_name in df.columns:
            axes[i,0].plot(df["time"], df[col_name], label=f"{r}")
    axes[i,0].set_title(f"Lambda {l} φ")
    axes[i,0].set_ylabel("Angle (rad)")
    if i == len(lambdas)-1:
        axes[i,0].set_xlabel("Time")
    axes[i,0].legend(fontsize=8, title="Residue")

    # ψ plot
    for r in residues:
        col_name = f'psi{r}'
        if col_name in df.columns:
            axes[i,1].plot(df["time"], df[col_name], label=f"{r}")
    axes[i,1].set_title(f"Lambda {l} ψ")
    if i == len(lambdas)-1:
        axes[i,1].set_xlabel("Time")

plt.tight_layout()

# Save figure 
plt.savefig("torsions_all_lambdas.png", dpi=300)
plt.close(fig)
print(f"Saved plot: torsions_all_lambdas.png")