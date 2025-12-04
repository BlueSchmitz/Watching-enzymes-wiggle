#!/usr/bin/env python3
# usage: python specific_activity_calculation.py activity_assay_data.xlsx

"""
Expects an Excel workbook where:
- Sheet 1 (data) contains time in the first column (seconds) and columns named like:
    Enzyme1_Rep1, Enzyme1_Rep2, Enzyme2_Rep1, ...
- Sheet 2 (lookup table) contains columns:
    Enzyme   Protein_concentration
  where `Enzyme` matches the enzyme prefix used in the data columns (e.g. "Enzyme1"),
  and `Protein_concentration` is the protein concentration to be used in the specific activity formula.

What this script does:
- Sliding window of 9 points (data measured every 15 s) -> window length = 9 points = 120 s
- For each window: linear fit (y = m x + b), compute slope m (Abs / sec) and R^2.
- For each replicate: average slopes for windows with R^2 >= 0.99 (report mean, SD, #windows)
- Per enzyme: average replicate means (report mean and SD across replicates)
- specific_activity = (dAbs * Vreaction) / (ext_coeff * dt * Vprotein * cprotein * l)
- Saves:
    - Excel: "activity_summary_<inputname>.xlsx" with sheets "Replicate_Slopes" and "Enzyme_Summary"
    - PNG plots per enzyme in "activity_plots/" showing raw replicates and accepted window fits.
"""

import sys
import os
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend

# ---------- Configurable variables ----------
SAMPLING_INTERVAL = 15.0       # seconds between data points
WINDOW_POINTS = 9              # number of points in sliding window
R2_THRESHOLD = 0.995           # threshold for accepting a window
# Specific activity constants
Vreaction = 0.0002             # liters
extinction_coefficient_NADH = 6220.0  # M^-1 cm^-1
Vprotein = 0.00001             # liters (volume of protein in reaction)
path_length_l = 0.5            # cm
# ------------------------------------------------

def linear_fit_and_r2(x, y):
    """
    Fit y = m x + b via least-squares and return slope m and R^2.
    x and y are 1D arrays.
    """
    if len(x) < 2:
        return np.nan, np.nan
    # Fit with numpy.polyfit
    m, b = np.polyfit(x, y, 1)
    y_pred = m * x + b
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    if ss_tot == 0:
        r2 = 0.0
    else:
        r2 = 1 - ss_res / ss_tot
    return float(m), float(r2)

def safe_float(x):
    """Convert to float and allow commas as decimal separators."""
    if pd.isna(x):
        return np.nan
    if isinstance(x, (float, int)):
        return float(x)
    s = str(x).strip()
    s = s.replace(",", ".")
    try:
        return float(s)
    except:
        return np.nan

def main(infile):
    p = Path(infile)
    if not p.exists():
        print("Input file not found:", infile)
        return

    # Read sheets: assume first is data, second is lookup
    xls = pd.ExcelFile(infile)
    if len(xls.sheet_names) < 2:
        raise ValueError("Excel workbook must have at least 2 sheets: data and concentration lookup")
    sheet_data = xls.sheet_names[0]
    sheet_lookup = xls.sheet_names[1]
    print("Using data sheet:", sheet_data, "and lookup sheet:", sheet_lookup)

    df_data = pd.read_excel(infile, sheet_name=sheet_data)
    df_lookup = pd.read_excel(infile, sheet_name=sheet_lookup)

    # Prepare lookup dict: key = enzyme (string), value = protein concentration (float)
    lookup_cols = {c.lower(): c for c in df_lookup.columns}
    # attempt to find sensible column names
    enzyme_col = df_lookup.columns[0]
    conc_col = df_lookup.columns[1]

    lookup = {}
    for _, row in df_lookup.iterrows():
        en = str(row[enzyme_col]).strip()
        con = safe_float(row[conc_col])
        lookup[en] = con

    # Time vector from first column
    time_col = df_data.columns[0]
    time = df_data[time_col].values.astype(float)

    # Build list of measurement columns (skip the first column)
    meas_cols = list(df_data.columns[1:])

    # Group columns by enzyme prefix before "_Rep"
    def enzyme_key(colname):
        # try to split on "_Rep" pattern
        if "_rep" in colname.lower():
            # cut at the position of "_rep" (case-insensitive)
            idx = colname.lower().index("_rep")
            return colname[:idx]
        # else split at last underscore
        if "_" in colname:
            return colname.rsplit("_", 1)[0]
        return colname

    enzymes = sorted({enzyme_key(c) for c in meas_cols})

    # Output containers
    replicate_rows = []   # per-replicate summary
    window_rows = []      # each window slope/R2 if you want to save (optional)
    enzyme_summary_rows = []

    # Create output plot folder
    plot_dir = Path("activity_plots")
    plot_dir.mkdir(exist_ok=True)

    # For each enzyme
    for enz in enzymes:
        # find replicate columns starting with this enzyme key
        reps = [c for c in meas_cols if enzyme_key(c) == enz]
        if len(reps) == 0:
            continue
        print(f"\nProcessing enzyme: {enz}  ({len(reps)} replicates)")

        replicate_means = []  # mean slope per replicate to later average across replicates
        replicate_sd = []     # per-rep SD of accepted window slopes
        replicate_counts = []

        # Prepare colors for plotting
        cmap = plt.get_cmap("viridis", max(1, len(reps)))

        # Prepare figure
        fig, ax = plt.subplots(figsize=(8,5))
        ax.set_title(enz)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Absorbance at 260 nm")

        # Collect all y for replicate-SD across temps (aligned)
        all_y_reps = []

        for i, rep in enumerate(reps):
            y = df_data[rep].values.astype(float)
            all_y_reps.append(y)
            color = cmap(i)

            # Sliding windows
            window_slopes = []
            window_r2s = []
            window_indices = []   # list of (start_idx, end_idx) for accepted windows

            n_points = len(y)
            for start in range(0, n_points - WINDOW_POINTS + 1): # range of possible windows
                end = start + WINDOW_POINTS
                x_win = time[start:end].astype(float) # collect time values
                y_win = y[start:end].astype(float) # collect absorbance values

                # Skip windows that contain NaNs (because number of measurements differ)
                if np.isnan(x_win).any() or np.isnan(y_win).any():
                    continue

                # Fit linear function and compute R2
                m, r2 = linear_fit_and_r2(x_win, y_win)
                window_rows.append({
                    "Enzyme": enz,
                    "Replicate": rep,
                    "start_time": float(time[start]),
                    "end_time": float(time[end-1]),
                    "slope_Abs_per_s": m,
                    "R2": r2,
                    "Accepted": r2 >= R2_THRESHOLD
                })
                if r2 >= R2_THRESHOLD:
                    window_slopes.append(m)
                    window_r2s.append(r2)
                    window_indices.append((start, end))

            # summarize replicate
            if len(window_slopes) > 0:
                mean_slope = float(np.mean(window_slopes))
                sd_slope = float(np.std(window_slopes, ddof=0))
                n_windows = len(window_slopes)
            else:
                mean_slope = np.nan
                sd_slope = np.nan
                n_windows = 0

            replicate_rows.append({
                "Enzyme": enz,
                "Replicate": rep,
                "N_windows_used": n_windows,
                "Mean_slope_Abs_per_s": mean_slope,
                "SD_slope_Abs_per_s": sd_slope
            })

            # Average per enzyme:
            if not np.isnan(mean_slope):
                replicate_means.append(mean_slope)
                replicate_sd.append(sd_slope)
                replicate_counts.append(n_windows)

            # --- Plotting ---
            # raw points
            ax.plot(time, y, 'o', markersize=3, color=color, alpha=0.7, label=f'{rep} data')

            # plot accepted window as slope triangles
            for (start, end) in window_indices:
                x_win = time[start:end].astype(float)
                y_win = y[start:end].astype(float)
    
                # linear fit for the window
                b = np.polyfit(x_win, y_win, 1)[1]  # get intercept
                y_fit = m * x_win + b
    
                # Plot the dashed line representing the slope
                ax.plot(x_win, y_fit, '--', color=color, linewidth=1, alpha=0.9)
    
                # Draw vertical lines to complete the triangle
                ax.vlines(x_win[0], y_win[0], y_fit[0], color=color, linewidth=1, alpha=0.9, linestyles='dashed')
                ax.vlines(x_win[-1], y_win[-1], y_fit[-1], color=color, linewidth=1, alpha=0.9, linestyles='dashed')

            # Also optionally plot the replicate-mean slope as a dashed line (across full time span)
            if not np.isnan(mean_slope):
                # compute intercept by fitting slope to middle of data for plotting only
                # b_mean chosen so that the line crosses the mean of y at mean of time
                b_mean = np.mean(y) - mean_slope * np.mean(time)
                ax.plot(time, mean_slope * time + b_mean, '-', color=color, alpha=0.7, linewidth=1)

        # End replicates loop for this enzyme

        # enzyme-level summary: average across replicate means
        if len(replicate_means) > 0:
            enzyme_mean_slope = float(np.mean(replicate_means))
            enzyme_sd_of_replicates = float(np.std(replicate_means, ddof=0))
            n_reps_success = len(replicate_means)
        else:
            enzyme_mean_slope = np.nan
            enzyme_sd_of_replicates = np.nan
            n_reps_success = 0

        # Compute specific activity using enzyme_mean_slope (slope in Abs/sec)
        # Formula: specific_activity = (dAbsorbance*Vreaction)/(ext_coeff*dt*Vprotein*cprotein*l)
        # With dAbsorbance/dt = enzyme_mean_slope (Abs / sec)
        # Therefore formula simplifies to:
        # specific_activity = (enzyme_mean_slope * 60 * Vreaction) / (extinction_coefficient_NADH * 1 * Vprotein * 1000 * cprotein * l)
        cprotein_value = lookup.get(enz, np.nan)
        if np.isnan(cprotein_value):
            print(f"Warning: concentration for enzyme '{enz}' not found in lookup table. specific_activity will be NaN.")
            specific_activity = np.nan
        else:
            specific_activity = (enzyme_mean_slope * 60 * Vreaction * 1000000) / (extinction_coefficient_NADH * Vprotein * float(cprotein_value) * 1000 * path_length_l)
            # *60 to convert sec to min, *1,000,000 to convert M to µM in numerator
            # *1000 to convert mg/mL to g/L in denominator
        # Save enzyme summary
        enzyme_summary_rows.append({
            "Enzyme": enz,
            "N_replicates_total": len(reps),
            "N_replicates_with_valid_slope": n_reps_success,
            "Mean_slope_all_replicates_Abs_per_s": enzyme_mean_slope,
            "SD_slope_across_replicates_Abs_per_s": enzyme_sd_of_replicates,
            "Protein_concentration_lookup_value": cprotein_value,
            "Specific_activity (units depend on conc units)": specific_activity
        })

        # Finalize plot
        ax.legend(fontsize="small", loc="best")
        ax.grid(alpha=0.3)
        outpng = plot_dir / f"{enz}_activity.png"
        fig.tight_layout()
        fig.savefig(outpng, dpi=300)
        plt.close(fig)
        print(f"Saved plot to {outpng}")

    # After all enzymes processed: write results to Excel
    out_xlsx = f"activity_summary_{p.stem}.xlsx"
    with pd.ExcelWriter(out_xlsx) as writer:
        pd.DataFrame(window_rows).to_excel(writer, sheet_name="Windows_all", index=False)
        pd.DataFrame(replicate_rows).to_excel(writer, sheet_name="Replicate_Slopes", index=False)
        pd.DataFrame(enzyme_summary_rows).to_excel(writer, sheet_name="Enzyme_Summary", index=False)

    print("Saved summary Excel to", out_xlsx)
    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python activity_analysis.py activity_file.xlsx")
        sys.exit(1)
    main(sys.argv[1])
