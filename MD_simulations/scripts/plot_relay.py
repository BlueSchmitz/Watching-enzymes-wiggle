import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
cmap = plt.get_cmap("viridis")


dfs = []
files = sorted(glob.glob("relay_geometry_*.csv"))
if len(files) == 0:
    raise FileNotFoundError("No relay_geometry_*.csv files found.")

for rep, file in enumerate(files):

    df = pd.read_csv(file)

    df["replicate"] = rep

    df["pair"] = (
        df["resname"] + "_"
        + df["resid"].astype(str) + "_"
        + df["atom"] + "_"
        + df["hydrogen"]
    )

    dfs.append(df)  
print(f"Loaded {len(dfs)} replicate(s)")

for rep, df in enumerate(dfs):
    df["frame_global"] = df["frame"] + rep * 100000  
    df["water_uid"] = (
        df["replicate"].astype(str)
        + "_"
        + df["water_resid"].astype(str)
    )

df_full = pd.concat(dfs, ignore_index=True)

print(df_full["pair"].unique())

# helper function to extract lifetimes from a frame of events
def get_lifetimes(frame_series):

    frames = np.sort(frame_series.unique())

    lifetimes = []
    current = 1

    for i in range(1, len(frames)):

        if frames[i] == frames[i-1] + 1:
            current += 1

        else:
            lifetimes.append(current)
            current = 1

    if len(frames) > 0:
        lifetimes.append(current)

    return lifetimes


# loop through all pairs and plot distance, angle, free energy landscape, and occupancy
for pair, g in df_full.groupby("pair"):

    plt.figure(figsize=(5,4))
    plt.scatter(g["distance_HOw"], g["distance_OwA"], s=5, alpha=0.7, rasterized=True)
    plt.xlabel("Ow-H distance (Å)")
    plt.ylabel("Ow-A distance (Å)")
    plt.title(pair)
    plt.tight_layout()
    os.makedirs(f"./relay/{pair}", exist_ok=True)
    plt.savefig(f"./relay/{pair}/{pair}_distances_relay.pdf")
    plt.close()

    plt.figure(figsize=(5,4))
    plt.hist2d(g["distance_HOw"], g["distance_OwA"], bins=50, density=True, cmap=cmap)
    plt.xlabel("Ow–H distance (Å)")
    plt.ylabel("Ow–A distance (Å)")
    plt.title(pair)
    plt.colorbar(label="Probability density")
    plt.tight_layout()
    plt.savefig(f"./relay/{pair}/{pair}_distances_relay_2d.pdf")
    plt.close()

    plt.figure(figsize=(5,4))
    plt.scatter(g["angle_CHOw"], g["angle_HOwA"], s=5, alpha=0.7, rasterized=True)
    plt.xlabel("C2-H-Ow Angle (deg)")
    plt.ylabel("H-Ow-A Angle (deg)")
    plt.title(pair)
    plt.tight_layout()
    plt.savefig(f"./relay/{pair}/{pair}_angle_relay.pdf")
    plt.close()

    plt.figure(figsize=(5,4))
    plt.hist2d(g["angle_CHOw"], g["angle_HOwA"], bins=50, density=True, cmap=cmap)
    plt.xlabel("C2-H-Ow Angle (deg)")
    plt.ylabel("H-Ow-A Angle (deg)")
    plt.title(pair)
    plt.colorbar(label="Probability density")
    plt.tight_layout()
    plt.savefig(f"./relay/{pair}/{pair}_angle_relay_2d.pdf")
    plt.close()

    # water probabilities
    counts = g["water_resid"].value_counts()
    p_i = counts / counts.sum()
    prob_df = pd.DataFrame({
        "water_resid": p_i.index,
        "probability": p_i.values
    })
    prob_df.to_csv(
        f"./relay/{pair}/{pair}_water_probabilities.csv",
        index=False
    )
    plt.figure(figsize=(5,4))
    plt.bar(
        p_i.index.astype(str),
        p_i.values
    )
    plt.xlabel("Water resid")
    plt.ylabel("Probability")
    plt.title(pair)
    plt.tight_layout()
    plt.savefig(f"./relay/{pair}/{pair}_water_probabilities.pdf")
    plt.close()

    # lifetime per water per pair
    water_stats = []
    for water, wdf in g.groupby("water_uid"):
        lifetimes = get_lifetimes(wdf["frame"])
        if len(lifetimes) == 0:
            continue
        water_stats.append({
            "water_resid": water,
            "probability": len(wdf) / len(g),
            "mean_lifetime_ps": np.mean(lifetimes) * 2,
            "median_lifetime_ps": np.median(lifetimes) * 2,
            "max_lifetime_ps": np.max(lifetimes) * 2,
            "n_events": len(lifetimes)
        })
    water_df = pd.DataFrame(water_stats)
    water_df.to_csv(
        f"./relay/{pair}/{pair}_water_statistics.csv",
        index=False
    )
    # plot lifetime and probability
    fig, ax1 = plt.subplots(figsize=(5,4))
    ax1.bar(
        water_df["water_resid"].astype(str),
        water_df["probability"],
        label='Probability'
    )
    ax1.set_ylabel("Probability")
    ax2 = ax1.twinx()
    ax2.plot(
        water_df["water_resid"].astype(str),
        water_df["mean_lifetime_ps"],
        marker="o",
        color='black',
        label='Mean lifetime (ps)'
    )
    ax2.set_ylabel("Mean lifetime (ps)")
    plt.legend()
    plt.title(pair)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"./relay/{pair}/{pair}_water_prob_lifetime.pdf")
    plt.close()

# Lifetimes for A-H pairs in stacked histograms
groups = list(df_full.groupby("pair"))

all_lifetimes = []
labels = []

for pair, g in groups:
    lt_frames = get_lifetimes(g["frame_global"])
    if len(lt_frames):
        lt_ps = np.array(lt_frames) * 2 # convert to ps
        all_lifetimes.append(lt_ps)
        labels.append(pair)

if len(all_lifetimes):
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
plt.savefig("./relay/lifetimes_histogram_stacked_relay.pdf")
plt.close()

# Lifetime of A-H pairs
lifetime_data = []

for rep, df in enumerate(dfs):

    for pair, g in df.groupby("pair"):

        lifetimes = get_lifetimes(g["frame"])
        if len(lifetimes) == 0:
            continue
        lifetimes_ps = np.array(lifetimes) * 2 # convert from frames to ps (2 ps per frame)
        lifetime_data.append({
            "replicate": rep,
            "pair": pair,
            "mean_lifetime_ps": np.mean(lifetimes_ps)
        })

lifetime_df = pd.DataFrame(lifetime_data)

lifetime_summary = (
    lifetime_df.groupby("pair")["mean_lifetime_ps"]
    .agg(["mean","std"])
    .reset_index()
)

plt.figure(figsize=(5,4))
plt.bar(
    lifetime_summary["pair"],
    lifetime_summary["mean"],
    yerr=lifetime_summary["std"],
    capsize=3
)
plt.ylabel("Mean lifetime (ps)")
plt.xlabel("A-H pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("./relay/mean_lifetime_relay.pdf")
plt.close()

lifetime_summary.to_csv(
    "./relay/mean_lifetime_relay_summary.csv",
    index=False
)

# plot occupancy bar chart for all pairs
occ_data = []
n_data = []

for rep, df in enumerate(dfs):

    # number of events per pair (k)
    k = df.groupby("pair")["frame"].nunique()

    # total frames per replicate (n)
    n_frames = 50000

    for pair, k_val in k.items():
        occ_data.append({
            "replicate": rep,
            "pair": pair,
            "k": k_val,
            "n": n_frames
        })

occ_df = pd.DataFrame(occ_data)

occ_df["p"] = occ_df["k"] / occ_df["n"]

# binomial variance per replicate
occ_df["var"] = occ_df["p"] * (1 - occ_df["p"]) / occ_df["n"]

# Pooled estimator
pooled = occ_df.groupby("pair").agg(
    k_total=("k", "sum"),
    n_total=("n", "sum")
)

pooled["p"] = pooled["k_total"] / pooled["n_total"]
pooled["std"] = np.sqrt(
    pooled["p"] * (1 - pooled["p"]) / pooled["n_total"]
)

# replicate variability
rep_stats = occ_df.groupby("pair")["p"].agg(["mean", "std"]).rename(
    columns={"mean": "rep_mean", "std": "rep_std"}
)

summary = pooled.join(rep_stats)

plt.figure(figsize=(5, 4))
plt.bar(
    summary.index,
    summary["p"] * 100,
    yerr=summary["std"] * 100,
    capsize=3
)
plt.ylabel("Occupancy (%)")
plt.xlabel("A-H Pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("./relay/occupancy_relay.pdf")
plt.close()

summary.reset_index().to_csv(
    "./relay/occupancy_relay.csv",
    index=False
)

entropy_data = []

for rep, df in enumerate(dfs):

    for pair, g in df.groupby("pair"):

        counts = g["water_resid"].value_counts()

        p_i = counts / counts.sum()

        # Shannon entropy
        S = -(p_i * np.log(p_i)).sum()

        # normalized entropy
        N = len(p_i)
        S_norm = S / np.log(N) if N > 1 else 0.0

        entropy_data.append({
            "replicate": rep,
            "pair": pair,
            "entropy": S,
            "norm_entropy": S_norm
        })

entropy_df = pd.DataFrame(entropy_data)

# Aggregate across replicates
entropy_summary = (
    entropy_df.groupby("pair")["entropy"]
    .agg(["mean", "std", "count"])
    .reset_index()
)

entropy_norm_summary = (
    entropy_df.groupby("pair")["norm_entropy"]
    .agg(["mean", "std", "count"])
    .reset_index()
)

# Missing values
entropy_summary["std"] = entropy_summary["std"].fillna(0.0)
entropy_norm_summary["std"] = entropy_norm_summary["std"].fillna(0.0)

entropy_summary.to_csv("./relay/water_entropy_per_pair.csv", index=False)
entropy_norm_summary.to_csv("./relay/water_norm_entropy_per_pair.csv", index=False)

# Plot entropy
entropy_summary_sorted = entropy_summary.sort_values("mean")

plt.figure(figsize=(5, 4))
plt.bar(
    entropy_summary_sorted["pair"],
    entropy_summary_sorted["mean"],
    yerr=entropy_summary_sorted["std"],
    capsize=3
)
plt.ylabel("Water entropy S")
plt.xlabel("A-H pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("./relay/water_entropy_barplot.pdf")
plt.close()

# normalized entropy
entropy_norm_sorted = entropy_norm_summary.sort_values("mean")

plt.figure(figsize=(5, 4))
plt.bar(
    entropy_norm_sorted["pair"],
    entropy_norm_sorted["mean"],
    yerr=entropy_norm_sorted["std"],
    capsize=3
)
plt.ylabel("Normalized water entropy S")
plt.xlabel("A-H pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("./relay/water_norm_entropy_barplot.pdf")
plt.close()