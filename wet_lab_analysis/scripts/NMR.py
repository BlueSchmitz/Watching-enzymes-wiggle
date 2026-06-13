#!/usr/bin/env python3
# usage: python NMR.py file.xlsx

"""
Expects an Excel workbook where:
- Sheet 1 has columns "time" and "ratio" (intensity ratio of two peaks over time)
What it does:
- Fits a monoexponential decay to the data: y(t) = y0 * exp(-R0 * t)
- Extracts parameters y0 (initial ratio) and R0 (decay rate), and calculates τ = 1/R0 (characteristic time)
- Plots the data and fit, and prints the parameters with uncertainties.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42

# Load data from Excel
file_name = sys.argv[1]
sheet_name = 0  

df = pd.read_excel(file_name, sheet_name=sheet_name)

t = (df["time"].values)/60  # convert seconds to minutes
try:
    y = df["ratio"].values
except:
    y=df["integrals"].values
    print("min of y:", np.min(y))
    print("max of y:", np.max(y))
    # convert to percentage of max value
    #y = y / np.max(y)

# Monoexponential decay model
# y(t) = y0 * exp(-R0 * t) + C (with offset)
def monoexp_offset(t, A, R0, C):
    return A * np.exp(-R0 * t) + C

# Fit
p0 = [y[0], 1e-3, 0]  # initial guesses

params, cov = curve_fit(monoexp_offset, t, y, p0=p0)

A_fit, R0_fit, C_fit = params
perr = np.sqrt(np.diag(cov))

A_err, R0_err, C_err = perr

tau = 1.0 / R0_fit
tau_err = R0_err / (R0_fit ** 2)

# Print results
print(f"A   = {A_fit:.4f} ± {A_err:.4f}")
print(f"R0  = {R0_fit:.6e} ± {R0_err:.6e} (1/min)")
print(f"C   = {C_fit:.4f} ± {C_err:.4f}")
print(f"τ   = {tau:.2f} ± {tau_err:.2f} (min)")

# Plot
t_fit = np.linspace(t.min(), t.max(), 500)
y_fit = monoexp_offset(t_fit, A_fit, R0_fit, C_fit)
# Predicted values at data points
y_pred = monoexp_offset(t, A_fit, R0_fit, C_fit)
# R^2 calculation
ss_res = np.sum((y - y_pred) ** 2)
ss_tot = np.sum((y - np.mean(y)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

plt.figure(figsize=(7, 5))
plt.scatter(t, y, label="Data", color="black")
plt.plot(t_fit, y_fit, label=f"$y(t) = {A_fit:.3f} e^{{-{R0_fit:.3e} t}} + {C_fit:.3f}$", color="red")

plt.xlabel("Time [min]")
plt.ylabel("Intensity ratio")

# Text box with fit results
fit_text = (
    rf"$\tau = {tau:.1f} \pm {tau_err:.1f}$ min" "\n"
    rf"$R^2 = {r_squared:.4f}$"
)

plt.text(
    0.75, 0.8, fit_text,          # x, y in axes coordinates
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

output_figure = f"HD_exchange_fit_{file_name.split('.')[0]}.pdf"
plt.legend()
plt.tight_layout()
plt.savefig(output_figure)
plt.close()