#!/usr/bin/env python3
# usage: python NMR_Lambert.py file.xlsx

"""
Expects an Excel workbook where:
- Sheet 1 has columns "Time" and "Concentration"
What it does:
- Fits a monoexponential decay to the data: y(t) = y0 * exp(-R0 * t)
- Extracts parameters y0 (initial ratio) and R0 (decay rate), and calculates τ = 1/R0 (characteristic time)
- Plots the data and fit, and prints the parameters with uncertainties.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import lambertw
import sys

# Load data from Excel
file_name = sys.argv[1]
sheet_name = 0  
enzyme = file_name.split(".")[0]
print(enzyme)

df = pd.read_excel(file_name, sheet_name=sheet_name)

t = (df["Time"].to_numpy(float))/60  # convert seconds to minutes
Sobs = df["Concentration_h"].to_numpy(float)

S0 = float(Sobs[0])

#def W_approx(x):
#    return np.log1p(x) * (1 - np.log1p(np.log1p(x)) / (2 + np.log1p(x)))

def model(t, Vmax, Km):
    x = (S0 / Km) * np.exp((S0 - Vmax * t) / Km)
    return np.real(Km * lambertw(x))

Vmax0 = (Sobs[0] - Sobs[-1]) / (t[-1] - t[0])
Km0 = S0 / 2
p0 = [Vmax0, Km0]

bounds = (
    [0, 1e-9],
    [np.inf, np.inf]
)

#dt = np.gradient(t)
#sigma = 1 / dt
#sigma = sigma / np.mean(sigma)

popt, pcov = curve_fit(
    model,
    t,
    Sobs,
    p0=p0,
    #sigma=sigma,
    absolute_sigma=False,
    bounds=bounds,
    maxfev=20000
)
Vmax_fit, Km_fit = popt

corr = pcov[0,1] / np.sqrt(pcov[0,0]*pcov[1,1])
print(corr)

# Timescale diagnostics
tS = (Km_fit + S0) / Vmax_fit
tQ = (27 * Km_fit * S0) / (4 * Vmax_fit)

ratio = tQ / tS

print("\n--- Timescale diagnostics ---")
print(f"tS = {tS:.4f} min")
print(f"tQ = {tQ:.4f} min")
print(f"tQ/tS = {ratio:.6e}")

# Plot
t_fit = np.linspace(t.min(), t.max(), 500)
S_fit = model(t_fit, Vmax_fit, Km_fit)

plt.figure(figsize=(6,4))
plt.plot(t, Sobs, "o", label="data")
plt.plot(t_fit, S_fit, "-", label="Lambert-W fit")
plt.xlabel("Time (min)")
plt.ylabel("Concentration (mM)")
plt.legend()
plt.tight_layout()
plt.savefig(f"lambertw_progress_curve_{enzyme}.png", dpi=300)

# Print results
print(f"Vmax = {Vmax_fit:.4f} ± {np.sqrt(pcov[0, 0]):.4f} mM/min")
print(f"Km   = {Km_fit:.4f} ± {np.sqrt(pcov[1, 1]):.4f} mM")

v = -np.gradient(Sobs, t)

plt.figure(figsize=(6,4))
plt.plot(Sobs, v, 'o')
plt.xlabel("[S]")
plt.ylabel("Rate")
plt.savefig(f"gradient.png", dpi=300)

plt.figure(figsize=(6,4))
plt.semilogy(t, Sobs, "o")
plt.savefig(f"semilogy.png", dpi=300)