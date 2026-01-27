#!/usr/bin/env python3
# usage: python selective_activity_analysis.py activity_assay_data.xlsx

"""
This script calculates specific enzyme activities from time course data,
averages per replicate and enzyme, and generates boxplots with selective
statistical comparisons.

Only the following comparisons are tested:
1-2, 1-3, 1-4, 3-4, 1-5, 5-6, 7-8, 7-1
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib

matplotlib.use("Agg")

# ---------- Configurable variables ----------
SAMPLING_INTERVAL = 16.0
WINDOW_POINTS = 12
R2_THRESHOLD = 0.999
Vreaction = 0.0002
extinction_coefficient_NADH = 6220.0
Vprotein = 0.00001
path_length_l = 1
# ------------------------------------------------

SELECTIVE_COMPARISONS = [
    ("1", "2"),
    ("1", "3"),
    ("1", "4"),
    ("3", "4"),
    ("1", "5"),
    ("5", "6"),
    ("7", "8"),
    ("7", "1"),
]

# ---------- Functions ----------

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

def plot_activity_boxplot(df_summary, df_replicates, comparisons, outfile_png):
    enzymes = df_summary["Enzyme"].astype(str).tolist()
    x = np.arange(len(enzymes))

    fig, ax = plt.subplots(figsize=(8,5))

    data = [df_replicates[df_replicates["Enzyme"]==e]["Specific_activity"].dropna().values
            for e in enzymes]

    # Boxplot
    bp = ax.boxplot(data, positions=x, widths=0.6, patch_artist=True, showmeans=False)
    box_color = plt.cm.viridis(0.5)
    for patch in bp['boxes']:
        patch.set_facecolor(box_color)
        patch.set_alpha(0.7)
    for element in ['whiskers', 'caps', 'medians', 'fliers']:
        plt.setp(bp[element], color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(enzymes)
    ax.set_ylabel("Specific activity")
    ax.set_title("Enzyme activities")

    # Perform only selected comparisons
    results = []
    for g1, g2 in comparisons:
        if g1 in enzymes and g2 in enzymes:
            d1 = df_replicates[df_replicates["Enzyme"]==g1]["Specific_activity"].dropna()
            d2 = df_replicates[df_replicates["Enzyme"]==g2]["Specific_activity"].dropna()
            if len(d1)>0 and len(d2)>0:
                t_stat, p_val = ttest_ind(d1, d2)
                results.append({"group1": g1, "group2": g2, "p_raw": p_val})

    if results:
        pvals = [r["p_raw"] for r in results]
        reject, pvals_corrected, _, _ = multipletests(pvals, method="holm", alpha=0.05)
        for i, r in enumerate(results):
            r["p_adj"] = pvals_corrected[i]
            r["reject"] = reject[i]

    # Add significance stars
    def stars(p):
        if p < 0.001: return "***"
        elif p < 0.01: return "**"
        elif p < 0.05: return "*"
        return ""

    y_max = max([max(d) if len(d)>0 else 0 for d in data])
    step = y_max * 0.05
    height = y_max * 1.05

    if results:
        for r in reversed(results):
            if r["reject"]:
                g1, g2 = r["group1"], r["group2"]
                i, j = enzymes.index(g1), enzymes.index(g2)
                ax.plot([i,i,j,j],[height,height+step,height+step,height], color='black', lw=1)
                ax.text((i+j)/2, height+step*1.1, stars(r["p_adj"]), ha='center', va='bottom', color='black')
                height += step*2
    
    legend_handles = [
        Line2D([0], [0], color='none', marker='', label='* p<0.05'),
        Line2D([0], [0], color='none', marker='', label='** p<0.01'),
        Line2D([0], [0], color='none', marker='', label='*** p<0.001')
    ]

    ax.legend(handles=legend_handles, loc='upper right', fontsize='small', handlelength=0, handletextpad=0.2)

    fig.tight_layout()
    fig.savefig(outfile_png, dpi=300)
    plt.close(fig)

# ---------- Main ----------

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

    lookup = {normalize_enzyme_name(row[enzyme_col]): safe_float(row[conc_col])
              for _, row in df_lookup.iterrows()}

    time = df_data.iloc[:,0].values.astype(float)
    meas_cols = list(df_data.columns[1:])

    def enzyme_key(colname):
        if "_rep" in colname.lower():
            idx = colname.lower().index("_rep")
            return colname[:idx]
        if "_" in colname:
            return colname.rsplit("_",1)[0]
        return colname

    enzymes = sorted({enzyme_key(c) for c in meas_cols})

    replicate_rows = []
    window_rows = []
    enzyme_summary_rows = []

    plot_dir = Path("activity_plots")
    plot_dir.mkdir(exist_ok=True)

    # Compute replicate-specific slopes and specific activities
    all_activities = []

    for enz in enzymes:
        reps = [c for c in meas_cols if enzyme_key(c)==enz]
        replicate_activities = []

        for rep in reps:
            y = df_data[rep].values.astype(float)

            window_slopes = []
            for start in range(0, len(y)-WINDOW_POINTS+1):
                end = start+WINDOW_POINTS
                x_win = time[start:end]
                y_win = y[start:end]
                if np.isnan(x_win).any() or np.isnan(y_win).any():
                    continue
                m, r2 = linear_fit_and_r2(x_win, y_win)
                window_rows.append({
                    "Enzyme": enz,
                    "Replicate": rep,
                    "start_time": float(time[start]),
                    "end_time": float(time[end-1]),
                    "slope_Abs_per_s": m,
                    "R2": r2,
                    "Accepted": r2>=R2_THRESHOLD
                })
                if r2>=R2_THRESHOLD:
                    window_slopes.append(m)

            if window_slopes:
                mean_slope = float(np.mean(window_slopes))
            else:
                mean_slope = np.nan

            cprotein = lookup.get(enz, np.nan)
            replicate_activity = slope_to_specific_activity(mean_slope, cprotein)

            replicate_rows.append({
                "Enzyme": enz,
                "Replicate": rep,
                "Mean_slope_Abs_per_s": mean_slope,
                "Specific_activity": replicate_activity
            })

            if not np.isnan(replicate_activity):
                replicate_activities.append(replicate_activity)
                all_activities.append((enz, replicate_activity))

        # Enzyme summary
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

    # Save Excel
    df_summary = pd.DataFrame(enzyme_summary_rows)
    df_replicates = pd.DataFrame(replicate_rows)
    out_xlsx = f"activity_summary_{p.stem}.xlsx"
    with pd.ExcelWriter(out_xlsx) as writer:
        pd.DataFrame(window_rows).to_excel(writer, sheet_name="Windows_all", index=False)
        df_replicates.to_excel(writer, sheet_name="Replicate_Slopes", index=False)
        df_summary.to_excel(writer, sheet_name="Enzyme_Summary", index=False)
    print("Saved summary Excel to", out_xlsx)

    # Plot boxplot with selective comparisons
    plot_activity_boxplot(df_summary, df_replicates, SELECTIVE_COMPARISONS, "enzyme_activity_boxplot_setup.png")
    print("Saved boxplot PNG to enzyme_activity_boxplot_setup.png")


if __name__=="__main__":
    if len(sys.argv)<2:
        print("Usage: python selective_activity_analysis.py activity_file.xlsx")
        sys.exit(1)
    main(sys.argv[1])
