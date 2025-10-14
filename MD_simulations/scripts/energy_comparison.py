#!/usr/bin/env python3
# usage: python plot_energy_comparison.py energies.xvg energies_scaled.xvg > energy_diff.log
'''Compare and plot multiple energy terms from 2 GROMACS .xvg files.'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: python plot_energy_comparison.py <energies.xvg> <energies_scaled.xvg>")
    sys.exit(1)

file1, file2 = Path(sys.argv[1]), Path(sys.argv[2])
label1 = "Original"
label2 = "Scaled"

def read_xvg(fname):
    # Read multi-column .xvg file with legend info
    legends = []
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith('@ s') and 'legend' in line:
                legends.append(line.split('"')[1])
            elif line.startswith(('#', '@')):
                continue
            else:
                parts = line.split()
                if len(parts) > 1:
                    data.append([float(x) for x in parts])
    data = np.array(data)
    time = data[:,0]
    values = data[:,1:]
    return time, values, legends

# Load both files
t1, val1, leg1 = read_xvg(file1)
t2, val2, leg2 = read_xvg(file2)

if leg1 != leg2:
    print("Warning: legend mismatch between files. Matching by order.")
terms = leg1

# Plot comparison
plt.figure(figsize=(10,6))
colors = plt.cm.tab10(np.linspace(0,1,len(terms)))

print(f"\nComparing energy terms between:\n {file1.name} ({label1})\n and {file2.name} ({label2})\n")

for i, term in enumerate(terms):
    # compute mean/std
    m1, s1 = np.mean(val1[:,i]), np.std(val1[:,i])
    m2, s2 = np.mean(val2[:,i]), np.std(val2[:,i])
    diff = m2 - m1

    print(f"{term:15s} → {label1}: {m1:10.2f} ± {s1:6.2f} | {label2}: {m2:10.2f} ± {s2:6.2f} | Δ = {diff:8.2f}")

    plt.plot(t1, val1[:,i], color=colors[i], linestyle='-',  label=f"{term} ({label1})")
    plt.plot(t2, val2[:,i], color=colors[i], linestyle='--', label=f"{term} ({label2})")

plt.xlabel("Time (ps)")
plt.ylabel("Energy (kJ/mol)")
plt.title("Energy Term Comparison: Original vs Scaled Topology")
plt.legend(fontsize=8, ncol=2)
plt.grid(True, linestyle=':')
plt.tight_layout()

out_png = "energy_comparison.png"
plt.savefig(out_png, dpi=150)
plt.close()

print(f"\n Saved comparison plot: {out_png}\n")
