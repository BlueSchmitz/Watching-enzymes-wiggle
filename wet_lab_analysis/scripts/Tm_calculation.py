#!/usr/bin/env python3
# usage: python Tm_calculation.py Tm_data.xlsx

'''
Use this script to calculate melting temperatures (Tm) from fluorescence data in an Excel file
by fitting a Boltzman sigmoid curve and finding the highest derivative.
'''

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import matplotlib
matplotlib.use("Agg") # Use non-GUI backend (avoids errors)

# Define Boltzmann sigmoid (lower plateau, upper plateau, Tm, slope)
def boltzmann(x, A1, A2, Tm, k):
    return A1 + ((A2 - A1) / (1 + np.exp((Tm - x)/k)))

# Load Excel
import sys
file_path = sys.argv[1]
df = pd.read_excel(file_path)

temperatures = df.iloc[:, 0].values
enzymes = df.columns[1:]  # assume columns like Enzyme1_Rep1, Enzyme1_Rep2, etc.

# Create output folder
os.makedirs("Tm_plots", exist_ok=True)

# Prepare results dictionaries
summary_rows = []
fit_rows = []

for enzyme in sorted(set(col.split("_")[0] for col in enzymes)):
    reps = [col for col in enzymes if col.startswith(enzyme)]
    Tm_values = []
    success_count = 0
    fail_count = 0
    
    plt.figure()
    plt.title(enzyme)
    plt.xlabel("Temperature")
    plt.ylabel("Fluorescence")
    
    N = len(reps)
    color_map = matplotlib.colormaps["viridis"].resampled(N)  # Viridis colormap for replicates
    
    for i, rep in enumerate(reps):
        color = color_map(i)
        ydata = df[rep].values
        
        # Initial guess
        p0 = [min(ydata), max(ydata), temperatures[np.argmax(np.gradient(ydata))], 2]
        max_idx = np.argmax(ydata)
        x_fit_data = temperatures[:max_idx+1]
        y_fit_data = ydata[:max_idx+1]
        
        try:
            popt, pcov = curve_fit(boltzmann, x_fit_data, y_fit_data, p0=p0, maxfev=5000)
            perr = np.sqrt(np.diag(pcov))
            
            # Derivative-based Tm
            x_fit = np.linspace(min(x_fit_data), max(x_fit_data), 500)
            y_fit = boltzmann(x_fit, *popt)
            dydx = np.gradient(y_fit, x_fit)
            Tm_derivative = x_fit[np.argmax(dydx)]
            Tm_values.append(Tm_derivative)
            success_count += 1
            
            # Store fit parameters per replicate
            fit_rows.append({
                "Enzyme": enzyme,
                "Replicate": rep,
                "A1": popt[0],
                "A2": popt[1],
                "Tm_fit": popt[2],
                "k": popt[3],
                "σA1": perr[0],
                "σA2": perr[1],
                "σTm": perr[2],
                "σk": perr[3],
                "Tm_derivative": Tm_derivative
            })
            
            # Plot points and fit line
            plt.plot(temperatures, ydata, 'o', color=color, markersize=3, label=f'{rep} data')
            plt.plot(x_fit, y_fit, '-', color=color, label=f'{rep} fit')
            
            # Shaded fit uncertainty
            y_upper = boltzmann(x_fit, 
                                popt[0] + perr[0], 
                                popt[1] + perr[1], 
                                popt[2] + perr[2], 
                                popt[3] + perr[3])
            y_lower = boltzmann(x_fit, 
                                popt[0] - perr[0], 
                                popt[1] - perr[1], 
                                popt[2] - perr[2], 
                                popt[3] - perr[3])
            plt.fill_between(x_fit, y_lower, y_upper, color=color, alpha=0.4)
            
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
        plt.axvline(mean_Tm, color='red', linestyle='--', label=f'Mean Tm = {mean_Tm:.2f} ± {std_Tm:.2f}')
    else:
        mean_Tm = np.nan
        std_Tm = np.nan
    
    plt.legend()
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
