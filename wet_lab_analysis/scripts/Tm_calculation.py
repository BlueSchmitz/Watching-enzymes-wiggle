import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

# Define Boltzmann sigmoid (lower plateau, upper plateau, Tm, slope)
def boltzmann(x, A1, A2, Tm, k):
    return A1 + ((A2 - A1) / (1 + np.exp((Tm - x)/k)))

# Load Excel
file_path = "Tm.xlsx"
df = pd.read_excel(file_path)

temperatures = df.iloc[:, 0].values
enzymes = df.columns[1:]  # assume columns like Enzyme1_Rep1, Enzyme1_Rep2, etc.

# Create output folder
os.makedirs("Tm_plots", exist_ok=True)

# Prepare results dictionary
results_summary = {"Enzyme": [], "Mean_Tm": [], "SD_Tm": []}

for enzyme in set(col.split("_")[0] for col in enzymes):  # group replicates by enzyme
    reps = [col for col in enzymes if col.startswith(enzyme)]
    Tm_values = []
    
    plt.figure()
    plt.title(enzyme)
    plt.xlabel("Temperature")
    plt.ylabel("Fluorescence")
    
    for rep in reps:
        ydata = df[rep].values
        p0 = [min(ydata), max(ydata), temperatures[np.argmax(np.gradient(ydata))], 1]
        # Find index of max fluorescence
        max_idx = np.argmax(ydata)
        # Use only data up to the maximum
        x_fit_data = temperatures[:max_idx+1]
        y_fit_data = ydata[:max_idx+1]
        
        try:
            popt, pcov = curve_fit(boltzmann, x_fit_data, y_fit_data, p0=p0) # initial guess for Boltzmann parameters
            # curve_fit: nonlinear least-squares fitting to Boltzmann function
            # popt: best-fit parameters (A1, A2, Tm, k)
            # pcov: covariance matrix (used to estimate parameter uncertainties)
            perr = np.sqrt(np.diag(pcov)) # standard errors of each parameter (square root of the diagonal of pcov)
            
            # Derivative
            x_fit = np.linspace(min(x_fit_data), max(x_fit_data), 500) # temperature grid for smooth plotting and derivative calculation
            y_fit = boltzmann(x_fit, *popt) # fitted Boltzmann curve using the best-fit parameters
            dydx = np.gradient(y_fit, x_fit) # numerical derivative of the fited curve
            Tm_derivative = x_fit[np.argmax(dydx)] # finds the index of maximum derivative
            Tm_values.append(Tm_derivative)
            
            # Plot replicate
            plt.plot(temperatures, ydata, 'o', label=f'{rep} data')
            plt.plot(x_fit, y_fit, '-', label=f'{rep} fit')
            
        except Exception as e:
            print(f"Could not fit {rep}: {e}")
    
    # Plot mean Tm
    if Tm_values:
        mean_Tm = np.mean(Tm_values)
        std_Tm = np.std(Tm_values)
        plt.axvline(mean_Tm, color='r', linestyle='--', label=f'Tm = {mean_Tm:.2f} ± {std_Tm:.2f}')

        # Add to results
        results_summary["Enzyme"].append(enzyme)
        results_summary["Mean_Tm"].append(mean_Tm)
        results_summary["SD_Tm"].append(std_Tm)
    
    plt.legend()
    plt.savefig(f"Tm_plots/{enzyme}_Tm_plot.png", dpi=300)
    plt.close()
    
    print(f"{enzyme}: Tm = {mean_Tm:.2f} ± {std_Tm:.2f} °C")

# Save summary table to Excel
results_df = pd.DataFrame(results_summary)
results_df.to_excel("Tm_summary.xlsx", index=False)
print("Summary table saved to Tm_summary.xlsx")
