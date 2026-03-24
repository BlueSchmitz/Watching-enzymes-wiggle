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
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend

# ---------- Configurable variables ----------
SAMPLING_INTERVAL = 16.0       # seconds between data points
WINDOW_POINTS = 12             # number of points in sliding window
R2_THRESHOLD = 0.999           # threshold for accepting a window
# Specific activity constants
Vreaction = 0.0002             # liters
extinction_coefficient_NADH = 6220.0  # M^-1 cm^-1
Vprotein = 0.00001             # liters (volume of protein in reaction)
path_length_l = 1              # cm
# ------------------------------------------------

def linear_fit_and_r2(x, y):
    if len(x) < 2:
        return np.nan, np.nan
    m, b = np.polyfit(x, y, 1)
    y_pred = m * x + b
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 0.0 if ss_tot == 0 else 1 - ss_res / ss_tot
    return float(m), float(r2)


def safe_float(x):
    if pd.isna(x):
        return np.nan
    if isinstance(x, (float, int)):
        return float(x)
    s = str(x).strip().replace(",", ".")
    try:
        return float(s)
    except:
        return np.nan


def normalize_enzyme_name(x):
    try:
        return str(int(float(x)))
    except:
        return str(x).strip()


def slope_to_specific_activity(slope, cprotein):
    if np.isnan(slope) or np.isnan(cprotein):
        return np.nan
    return -1 * (
        (slope * 60 * Vreaction * 1_000_000) /
        (extinction_coefficient_NADH * Vprotein * float(cprotein) * 1000 * path_length_l)
    )

def plot_activity_boxplot(df_summary, df_replicates, tukey_df, outfile_png):
    enzymes = df_summary["Enzyme"].astype(str).tolist()
    x = np.arange(len(enzymes))

    fig, ax = plt.subplots(figsize=(8,5))

    # Data per enzyme
    data = [df_replicates[df_replicates["Enzyme"]==e]["Specific_activity"].dropna().values
            for e in enzymes]

    # Boxplot (remove green mean triangle)
    bp = ax.boxplot(data, positions=x, widths=0.6, patch_artist=True, showmeans=False)

    # Set all boxes same color
    box_color = plt.cm.viridis(0.5)
    for patch in bp['boxes']:
        patch.set_facecolor(box_color)
        patch.set_alpha(0.7)
    for element in ['whiskers', 'caps', 'medians', 'fliers']:
        plt.setp(bp[element], color='black')

    # Overlay individual points
    #for i, yvals in enumerate(data):
    #    ax.scatter(np.full_like(yvals, i) + np.random.uniform(-0.1,0.1,len(yvals)),
    #               yvals, color='black', zorder=10)

    ax.set_xticks(x)
    ax.set_xticklabels(enzymes)
    ax.set_ylabel("Specific activity (U/mg)")
    ax.set_title("Enzyme activities")

    '''
    # Add significance stars from Tukey
    legend_handles = []
    if not tukey_df.empty and "reject" in tukey_df.columns:
        y_max = max([max(d) if len(d)>0 else 0 for d in data])
        step = y_max * 0.07
        height = y_max * 1.05

        def stars(p):
            if p < 0.001: return "***"
            elif p < 0.01: return "**"
            elif p < 0.05: return "*"
            return ""

        sig_pairs = tukey_df[tukey_df["reject"] == True]
        for _, row in sig_pairs.iterrows():
            g1, g2 = str(row["group1"]), str(row["group2"])
            p = float(row["p-adj"])
            i, j = enzymes.index(g1), enzymes.index(g2)

            # bracket lines in black
            ax.plot([i, i, j, j], [height, height+step, height+step, height], lw=1, color='black')
            ax.text((i+j)/2, height+step*1.1, stars(p), ha='center', va='bottom', color='black')
            height += step * 2

        # Add legend for significance stars
        legend_handles = [
            Line2D([0], [0], color='none', marker='', label='* p<0.05'),
            Line2D([0], [0], color='none', marker='', label='** p<0.01'),
            Line2D([0], [0], color='none', marker='', label='*** p<0.001')
        ]

        ax.legend(handles=legend_handles, loc='upper right', fontsize='small', handlelength=0, handletextpad=0.2)
    '''
    fig.tight_layout()
    fig.savefig(outfile_png, dpi=300)
    plt.close(fig)

def main(infile):
    p = Path(infile)
    if not p.exists():
        print("Input file not found:", infile)
        return

    xls = pd.ExcelFile(infile)
    sheet_data = xls.sheet_names[0]
    sheet_lookup = xls.sheet_names[1]

    df_data = pd.read_excel(infile, sheet_name=sheet_data)
    df_lookup = pd.read_excel(infile, sheet_name=sheet_lookup)

    enzyme_col = df_lookup.columns[0]
    conc_col = df_lookup.columns[1]

    lookup = {}
    for _, row in df_lookup.iterrows():
        enz = normalize_enzyme_name(row[enzyme_col])
        lookup[enz] = safe_float(row[conc_col])

    time = df_data.iloc[:, 0].values.astype(float)
    meas_cols = list(df_data.columns[1:])

    def enzyme_key(colname):
        if "_rep" in colname.lower():
            idx = colname.lower().index("_rep")
            return colname[:idx]
        if "_" in colname:
            return colname.rsplit("_", 1)[0]
        return colname

    enzymes = sorted({enzyme_key(c) for c in meas_cols})

    replicate_rows = []
    window_rows = []
    enzyme_summary_rows = []
    all_activities = []   

    plot_dir = Path("activity_plots")
    plot_dir.mkdir(exist_ok=True)

    for enz in enzymes:
        reps = [c for c in meas_cols if enzyme_key(c) == enz]
        if not reps:
            continue

        replicate_activities = []

        cmap = plt.get_cmap("viridis", max(1, len(reps)))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.set_title(enz)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Absorbance")

        for i, rep in enumerate(reps):
            y = df_data[rep].values.astype(float)
            color = cmap(i)

            window_slopes = []
            window_indices = []

            for start in range(0, len(y) - WINDOW_POINTS + 1):
                end = start + WINDOW_POINTS
                x_win = time[start:end]
                y_win = y[start:end]

                if np.isnan(x_win).any() or np.isnan(y_win).any():
                    continue

                m, r2 = linear_fit_and_r2(x_win, y_win)

                window_rows.append({
                    "Enzyme": enz,
                    "Replicate": rep,
                    "start_time": float(time[start]),
                    "end_time": float(time[end - 1]),
                    "slope_Abs_per_s": m,
                    "R2": r2,
                    "Accepted": r2 >= R2_THRESHOLD
                })

                if r2 >= R2_THRESHOLD:
                    window_slopes.append(m)
                    window_indices.append((start, end))

            if window_slopes:
                mean_slope = float(np.mean(window_slopes))
                sd_slope = float(np.std(window_slopes, ddof=1))
                n_windows = len(window_slopes)
            else:
                mean_slope = np.nan
                sd_slope = np.nan
                n_windows = 0

            cprotein = lookup.get(enz, np.nan)
            replicate_activity = slope_to_specific_activity(mean_slope, cprotein)

            replicate_rows.append({
                "Enzyme": enz,
                "Replicate": rep,
                "N_windows_used": n_windows,
                "Mean_slope_Abs_per_s": mean_slope,
                "SD_slope_Abs_per_s": sd_slope,
                "Specific_activity": replicate_activity
            })

            if not np.isnan(replicate_activity):
                replicate_activities.append(replicate_activity)   # per-enzyme list
                all_activities.append((enz, replicate_activity))   # global list for stats

            # ---- PLOTTING ----
            ax.plot(time, y, '.', color=color, alpha=0.7, label=f'{rep} data')

            accepted_indices = []
            for start, end in window_indices:
                accepted_indices.extend(range(start, end))
            accepted_indices = np.unique(accepted_indices)

            if len(accepted_indices) > 0:
                ax.plot(
                    time[accepted_indices],
                    y[accepted_indices],
                    'v',
                    linestyle='None',
                    color=color,
                    markersize=6,
                    markeredgewidth=0,
                    label=f'{rep} accepted'
                )

            if not np.isnan(mean_slope) and len(accepted_indices) > 0:
                t_mid = np.mean(time[accepted_indices])
                y_mid = np.mean(y[accepted_indices])
                intercept = y_mid - mean_slope * t_mid
                y_fit_line = mean_slope * time + intercept
                ax.plot(time, y_fit_line, '-', linewidth=1, color=color, alpha=0.9, label=f'{rep} mean slope')

        if replicate_activities:
            enzyme_mean = float(np.mean(replicate_activities))
            enzyme_sd = float(np.std(replicate_activities, ddof=1))
            n_valid = len(replicate_activities)
        else:
            enzyme_mean = np.nan
            enzyme_sd = np.nan
            n_valid = 0

        enzyme_summary_rows.append({
            "Enzyme": enz,
            "N_replicates_total": len(reps),
            "N_replicates_with_valid_activity": n_valid,
            "Mean_specific_activity": enzyme_mean,
            "SD_specific_activity": enzyme_sd,
            "Protein_concentration_lookup_value": lookup.get(enz, np.nan)
        })

        ax.legend(fontsize="small")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(plot_dir / f"{enz}_activity.png", dpi=300)
        plt.close(fig)

    # ---------------- STATISTICS ----------------
    stats_rows = []

    if len(all_activities) > 0:
        df_stats = pd.DataFrame(all_activities, columns=["Enzyme", "Activity"])

        # One-way ANOVA
        grouped = [g["Activity"].values for _, g in df_stats.groupby("Enzyme")]
        f_stat, p_value = f_oneway(*grouped)

        stats_rows.append({
            "Test": "One-way ANOVA",
            "Group1": "",
            "Group2": "",
            "Statistic": f_stat,
            "p_value": p_value,
            "Significant(p<0.05)": p_value < 0.05
        })

        # Tukey HSD posthoc
        tukey = pairwise_tukeyhsd(
            endog=df_stats["Activity"],
            groups=df_stats["Enzyme"],
            alpha=0.05
        )

        tukey_df = pd.DataFrame(
            data=tukey.summary().data[1:],
            columns=tukey.summary().data[0]
        )

    else:
        tukey_df = pd.DataFrame()
    
    # ---- Save barplot as PNG ----
    df_summary = pd.DataFrame(enzyme_summary_rows)
    df_replicates = pd.DataFrame(replicate_rows)
    boxplot_png = f"enzyme_activity_boxplot_{p.stem}.png"
    plot_activity_boxplot(df_summary, df_replicates, tukey_df, boxplot_png)
    print(f"Saved boxplot to {boxplot_png}")

    out_xlsx = f"activity_summary_{p.stem}.xlsx"
    with pd.ExcelWriter(out_xlsx) as writer:
        pd.DataFrame(window_rows).to_excel(writer, sheet_name="Windows_all", index=False)
        pd.DataFrame(replicate_rows).to_excel(writer, sheet_name="Replicate_Slopes", index=False)
        pd.DataFrame(enzyme_summary_rows).to_excel(writer, sheet_name="Enzyme_Summary", index=False)

        if len(stats_rows) > 0:
            pd.DataFrame(stats_rows).to_excel(writer, sheet_name="ANOVA", index=False)
            tukey_df.to_excel(writer, sheet_name="Tukey_posthoc", index=False)

    print("Saved summary Excel to", out_xlsx)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python activity_analysis.py activity_file.xlsx")
        sys.exit(1)
    main(sys.argv[1])