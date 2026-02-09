#!/usr/bin/env python3
# usage: python Tm_calculation.py Tm_data.xlsx

'''
Use this script to calculate melting temperatures (Tm) from fluorescence data in an Excel file
by fitting a Boltzman sigmoid curve and finding the highest derivative.
'''

import sys
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import math
import matplotlib
matplotlib.use("Agg") # Use non-GUI backend (avoids errors)

# Define Boltzmann sigmoid (lower plateau, upper plateau, Tm, slope)
def sigmoid(x, L, x0, k, b):
    """
    L  = curve maximum
    x0 = midpoint (inflection point)
    k  = slope
    b  = baseline offset
    """
    return L / (1 + np.exp(-k * (x - x0))) + b

# Load Excel
file_path = sys.argv[1]
df = pd.read_excel(file_path)

xdata = df.iloc[:, 0].astype(str).str.replace(",", ".").astype(float)
enzymes = df.columns[1:] # assume columns like Enzyme1_Rep1, Enzyme1_Rep2, etc.

# Create output folder
os.makedirs("Tm_plots", exist_ok=True)

# Prepare results dictionaries
summary_rows = []
fit_rows = []

# Prepare multi-panel figure (automatic grid layout based on number of enzymes)
unique_enzymes = sorted(set(col.split("_")[0] for col in enzymes))

for enzyme in sorted(set(col.split("_")[0] for col in enzymes)):
    reps = [col for col in enzymes if col.startswith(enzyme)]
    Tm_values = []
    success_count = 0
    fail_count = 0
    
    plt.figure()
    plt.title(enzyme)
    plt.xlabel("Temperature [°C]")
    plt.ylabel("Fluorescence [RFU]") # relative fluorescence units

    cmap = plt.get_cmap("viridis", max(1, len(reps)))  # Viridis colormap for replicates
    
    for i, rep in enumerate(reps):
        color = cmap(i)
        ydata = df[rep].astype(str).str.replace(",", ".").astype(float).values

        # Plot raw data
        plt.plot(xdata, ydata, 'o', color=color, markersize=3, label=f'{rep} data')
        
        ydata_np = ydata.copy()  # ensure it's a NumPy array
        xdata_np = xdata.values  # convert Series to NumPy array

        # Keep only finite values
        mask = np.isfinite(ydata_np) & np.isfinite(xdata_np)
        xdata_np = xdata_np[mask]
        ydata_np = ydata_np[mask]

        max_idx = np.argmax(ydata_np)
        x_fit_data = xdata_np[:max_idx + 10]
        y_fit_data = ydata_np[:max_idx + 10]

        # Keep only finite values
        mask = np.isfinite(y_fit_data) & np.isfinite(x_fit_data)
        x_fit_data = x_fit_data[mask]
        y_fit_data = y_fit_data[mask]

        # Initial guesses for sigmoid parameters
        L_guess = max(ydata)
        x0_guess = x_fit_data[np.argmax(np.gradient(y_fit_data))]
        k_guess = 2
        b_guess = min(ydata)

        p0 = [L_guess, x0_guess, k_guess, b_guess]

        try:
            params, covariance = curve_fit(sigmoid, x_fit_data, y_fit_data, p0=p0, maxfev=10000)

            L, x0, k, b = params

            residuals = y_fit_data - sigmoid(x_fit_data, *params)
            rss = np.sum(residuals**2)

            # Derivative-based Tm
            x_fit = np.linspace(min(xdata), max(xdata), 500) # generate smooth temp values
            y_fit = sigmoid(x_fit, *params) # smooth fitted curve
            dydx = np.gradient(y_fit, x_fit) # numerical derivative
            Tm_derivative = x_fit[np.argmax(dydx)] # Tm from max derivative
            Tm_values.append(Tm_derivative)
            success_count += 1

            # Store fit parameters per replicate
            fit_rows.append({
                "Enzyme": enzyme,
                "Replicate": rep,
                "Baseline": b,
                "Amplitude": L,
                "Tm_fit": x0,
                "Slope": k,
                "σBaseline": np.sqrt(covariance[3,3]),
                "σAmplitude": np.sqrt(covariance[0,0]),
                "σTm": np.sqrt(covariance[1,1]),
                "σSlope": np.sqrt(covariance[2,2]),
                "Tm_derivative": Tm_derivative
            })


            # Plot fit line
            plt.plot(x_fit, y_fit, '-', color=color, label=f'{rep} fit')
            
        except Exception as e:
            print(f"Could not fit {rep}: {e}")
            fail_count += 1
            fit_rows.append({
                "Enzyme": enzyme,
                "Replicate": rep,
                "A1": np.nan,
                "A2": np.nan,
                "Tm_fit": np.nan,
                "k": np.nan,
                "σA1": np.nan,
                "σA2": np.nan,
                "σTm": np.nan,
                "σk": np.nan,
                "Tm_derivative": np.nan
            })
    
    # Plot mean Tm ± SD
    if Tm_values:
        mean_Tm = np.mean(Tm_values)
        std_Tm = np.std(Tm_values)
        plt.axvline(mean_Tm, color='red', linestyle='--', label=f'Tm = {mean_Tm:.2f} ± {std_Tm:.2f}')
    else:
        mean_Tm = np.nan
        std_Tm = np.nan

    plt.legend(fontsize="small", loc = 'upper right')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"Tm_plots/{enzyme}_Tm_plot.png", dpi=300)
    plt.close()

    summary_rows.append({
        "Enzyme": enzyme,
        "Mean_Tm": mean_Tm,
        "SD_Tm": std_Tm,
        "N_success": success_count,
        "N_failed": fail_count
    })

# Save to Excel with multiple sheets
summary_df = pd.DataFrame(summary_rows)
fits_df = pd.DataFrame(fit_rows)

output_name = f"Tm_summary_{'_'.join(summary_df['Enzyme'])}.xlsx"
with pd.ExcelWriter(output_name) as writer:
    summary_df.to_excel(writer, sheet_name="Enzyme_Summary", index=False)
    fits_df.to_excel(writer, sheet_name="Replicate_Fits", index=False)

print(f"Summary saved to {output_name}")
