import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
cmap = plt.get_cmap("viridis")

dfs = []
files = sorted(glob.glob("direct_geometry_*.csv"))
if len(files) == 0:
    raise FileNotFoundError("No direct_geometry_*.csv files found.")

for rep, file in enumerate(files):

    df = pd.read_csv(file)

    os.makedirs(str(rep), exist_ok=True)

    df["replicate"] = rep

    df["pair"] = (
        df["resname"] + "_"
        + df["resid"].astype(str) + "_"
        + df["atom"] + "_"
        + df["hydrogen"]
    )

    dfs.append(df)  

print(f"Loaded {len(dfs)} replicate(s)")
df_full = pd.concat(dfs, ignore_index=True)

for pair in df_full["pair"].unique():
    os.makedirs(pair, exist_ok=True)

print(df_full["pair"].unique())

# helper function to extract lifetimes of True values in a boolean series
def get_lifetimes(boolean_series):

    arr = boolean_series.to_numpy(dtype=bool)

    lifetimes = []
    current = 0

    for x in arr:
        if x:
            current += 1
        elif current > 0:
            lifetimes.append(current)
            current = 0

    if current > 0:
        lifetimes.append(current)

    return lifetimes

# loop through all pairs and plot distance, angle, free energy landscape, and occupancy
for rep, df in enumerate(dfs):
    for pair, g in df.groupby("pair"):

        plt.figure(figsize=(5,4))
        plt.scatter(g["time_ps"], g["distance_HA"], s=5, alpha=0.7, rasterized=True)
        plt.xlabel("Time (ps)")
        plt.ylabel("H-A distance (Å)")
        plt.title(pair)
        plt.tight_layout()
        plt.savefig(f"./{rep}/{pair}_distance.pdf")
        plt.close()

        plt.figure(figsize=(5,4))
        plt.scatter(g["time_ps"], g["angle_CHA"], s=5, alpha=0.7, rasterized=True)
        plt.xlabel("Time (ps)")
        plt.ylabel("C2-H-A Angle (deg)")
        plt.title(pair)
        plt.tight_layout()
        plt.savefig(f"./{rep}/{pair}_angle.pdf")
        plt.close()

        lifetimes = get_lifetimes(g["within_geometry"])
        if len(lifetimes) == 0:
            continue
        lifetimes_ps = np.array(lifetimes) * 2 # convert from frames to ps (2 ps per frame)
        plt.figure(figsize=(5,4))
        plt.hist(
            lifetimes_ps,
            bins=np.arange(
                1,
                max(lifetimes_ps)+2
            ) - 0.5
        )
        plt.xlabel("Lifetime (ps)")
        plt.ylabel("Count")
        plt.title(pair)
        plt.tight_layout()
        plt.savefig(f"./{rep}/{pair}_lifetimes.pdf")
        plt.close()


# use all data to plot free energy landscape and occupancy bar chart
groups = list(df_full.groupby("pair"))
for pair, g in groups:
    kB = 0.008314  # kJ/mol/K
    T = 298  # K
    H, xedges, yedges = np.histogram2d(
        g["distance_HA"],
        g["angle_CHA"],
        bins=50
    )
    P = H / np.max(H) # normalize to max bin count --> most populated bin F=0, others F>0
    F = -kB * T * np.log(P + 1e-12)  # avoid log(0)
    F[P == 0] = np.nan               # hide unsampled regions
    plt.figure(figsize=(5, 4))
    plt.imshow(F.T, origin='lower',
           extent=[xedges[0], xedges[-1],
                   yedges[0], yedges[-1]],
           aspect='auto')
    plt.colorbar(label="Free Energy (kJ/mol)")
    plt.xlabel("Distance (Å)")
    plt.ylabel("C2-H-A Angle (deg)")
    plt.title(pair)
    plt.tight_layout()
    plt.savefig(f"{pair}_free_energy.pdf")
    plt.close()
    
# plot occupancy bar chart
occ_list = []
n_list = []

for rep, df in enumerate(dfs):
    g = df.groupby("pair")["within_geometry"]
    occ_list.append(g.mean())
    n_list.append(g.count())

occ_df = pd.concat(occ_list, axis=1)
n_df   = pd.concat(n_list, axis=1)

# pooled estimator
k_total = (occ_df * n_df).sum(axis=1)
n_total = n_df.sum(axis=1)

occ_mean = k_total / n_total

# binomial error
occ_std_binom = np.sqrt(occ_mean * (1 - occ_mean) / n_total)

# replicate variability
occ_std_rep = occ_df.std(axis=1, ddof=1)

plt.figure(figsize=(5,4))
plt.bar(
    occ_mean.index,
    occ_mean * 100,
    yerr=occ_std_binom * 100,
    capsize=3
)
plt.ylabel("Occupancy (%)")
plt.xlabel("A-H Pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("occupancy.pdf")
plt.close()

# Lifetimes in stacked histograms
groups = list(df_full.groupby("pair"))

all_lifetimes = []
labels = []

for pair, g in groups:
    lt_frames = get_lifetimes(g["within_geometry"])
    if len(lt_frames):
        lt_ps = np.array(lt_frames) * 2 # convert to ps
        all_lifetimes.append(lt_ps)
        labels.append(pair)

max_lt = max(max(x) for x in all_lifetimes)
bins = np.arange(1, max_lt + 2) - 0.5

colors = plt.cm.viridis(
    np.linspace(1, 0, len(all_lifetimes))
)

plt.figure(figsize=(8,5))

plt.hist(
    all_lifetimes,
    bins=bins,
    density=True,      
    stacked=True,
    color=colors,
    edgecolor="black",
    linewidth=1,
    label=labels
)

plt.xlabel("Lifetime (ps)")
plt.ylabel("Probability density")
plt.legend()
plt.tight_layout()
plt.savefig("lifetimes_histogram_stacked.pdf")
plt.close()