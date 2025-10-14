#!/usr/bin/env python3
# usage: python plot_xvg.py pressure_1000.xvg [pressure_500.xvg ...]
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Handle command-line input
if len(sys.argv) < 2:
    print("Usage: python plot_xvg.py <file1.xvg> [file2.xvg ...]")
    sys.exit(1)

plt.figure(figsize=(8,5))

# Loop over all provided .xvg files
for filename in sorted(sys.argv[1:]):
    fpath = Path(filename)
    data = []
    with open(fpath) as f:
        for line in f:
            if line.startswith(('#', '@')): # skip comments
                continue
            parts = line.split()
            if len(parts) >= 2:
                data.append((float(parts[0]), float(parts[1])))

    if not data:
        print(f"Warning: no data in {fpath}")
        continue

    time, val = zip(*data)
    val = np.array(val)

    # Compute mean and standard deviation
    mean_val = np.mean(val)
    std_val = np.std(val)
    print(f"{fpath.stem:20s}: mean = {mean_val:10.3f}, std = {std_val:10.3f}")

    # Create legend label from the filename
    # Example: pressure_1000.xvg → "Restraint = 1000 kJ/mol/nm²"
    label = fpath.stem
    if "_" in label:
        restraint = label.split("_")[-1]
        label = f"{restraint} kJ mol⁻¹ nm⁻²"
    else:
        label = label

    plt.plot(time, val, label=label, linewidth=1.5)

# Auto-label axes and add target lines
first_name = sys.argv[1].lower()
if "pressure" in first_name:
    ylabel = "Pressure (bar)"
    title = "NPT Equilibration Pressure (Restraint Scaling)"
    plt.axhline(1.0, color='gray', linestyle='--', label='Target: 1 bar')
elif "density" in first_name:
    ylabel = "Density (kg/m³)"
    title = "NPT Equilibration Density (Restraint Scaling)"
elif "temperature" in first_name:
    ylabel = "Temperature (K)"
    title = "NVT Equilibration Temperature"
    plt.axhline(298, color='gray', linestyle='--', label='Target: 298 K')
elif "potential" in first_name:
    ylabel = "Potential Energy (kJ/mol)"
    title = "Energy Minimization Potential Energy"
else:
    ylabel = "Value"
    title = "Combined NPT Property"

plt.xlabel("Time (ps)")
plt.ylabel(ylabel)
plt.title(title)
plt.legend(title="Restraint", fontsize=9)
plt.grid(True)
plt.tight_layout()

# Save plot
prop = "pressure" if "pressure" in first_name else "density" if "density" in first_name else "temperature" if "temperature" in first_name else "potential" if "potential" in first_name else "value"
out_png = f"{prop}.png"
plt.savefig(out_png, dpi=150)
plt.close()
print(f"Saved plot: {out_png}")