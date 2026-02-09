#!/usr/bin/env python3
# usage: python Tm_calculation.py Tm_data.xlsx

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ---- Sigmoid function (4-parameter logistic) ----
def sigmoid(x, L, x0, k, b):
    """
    L  = curve maximum
    x0 = midpoint (inflection point)
    k  = slope
    b  = baseline offset
    """
    return L / (1 + np.exp(-k * (x - x0))) + b


# ---- Load Excel data ----
# Replace with your actual file name
file_path = "data/Test.xlsx"

df = pd.read_excel(file_path)

# Convert European decimal commas to dots
df["Temp"] = df["Temp"].astype(str).str.replace(",", ".").astype(float)
df["LmDERA_Rep1"] = df["LmDERA_Rep1"].astype(str).str.replace(",", ".").astype(float)

xdata = df["Temp"].values
ydata = df["LmDERA_Rep1"].values
max_idx = np.argmax(ydata)
x_fit_data = xdata[:max_idx+10]
y_fit_data = ydata[:max_idx+10]


# ---- Initial parameter guesses ----
L_guess = max(ydata)
x0_guess = x_fit_data[np.argmax(np.gradient(y_fit_data))]
k_guess = 0.5
b_guess = min(ydata)

p0 = [L_guess, x0_guess, k_guess, b_guess]


# ---- Curve fitting ----
params, covariance = curve_fit(sigmoid, x_fit_data, y_fit_data, p0=p0, maxfev=10000)

L, x0, k, b = params

print("Fitted parameters:")
print(f"L  (max)        = {L:.3f}")
print(f"x0 (midpoint)   = {x0:.3f}")
print(f"k  (slope)      = {k:.3f}")
print(f"b  (baseline)   = {b:.3f}")


# ---- Plot ----
x_fit = np.linspace(min(xdata), max(xdata), 500)
y_fit = sigmoid(x_fit, *params)

plt.figure(figsize=(8, 6))
plt.scatter(xdata, ydata, label="Data", s=25)
plt.plot(x_fit, y_fit, label="Sigmoid fit", linewidth=2)
plt.xlabel("Temperature")
plt.ylabel("LmDERA_Rep1")
plt.legend()
plt.tight_layout()
plt.show()
