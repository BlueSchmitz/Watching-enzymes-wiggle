import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
cmap = plt.get_cmap("viridis")


# Unique identifier for each acceptor-proton pair
df=pd.read_csv("relay_geometry.csv")

df["pair"] = (
    df["resname"] + "_"
    + df["resid"].astype(str) + "_"
    + df["atom"] + "_"
    + df["hydrogen"]
)

print(df["pair"].unique())

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
for pair, g in df.groupby("pair"):

    plt.figure(figsize=(5,4))
    plt.scatter(g["distance_HOw"], g["distance_OwA"], s=5, alpha=0.7, rasterized=True)
    plt.xlabel("Ow-H distance (Å)")
    plt.ylabel("Ow-A distance (Å)")
    plt.title(pair)
    plt.tight_layout()
    plt.savefig(f"{pair}_distances_relay.pdf")
    plt.close()

    plt.figure(figsize=(5,4))
    plt.hist2d(g["distance_HOw"], g["distance_OwA"], bins=50, density=True, cmap=cmap)
    plt.xlabel("Ow–H distance (Å)")
    plt.ylabel("Ow–A distance (Å)")
    plt.title(pair)
    plt.colorbar(label="Probability density")
    plt.tight_layout()
    plt.savefig(f"{pair}_distances_relay_2d.pdf")
    plt.close()

    plt.figure(figsize=(5,4))
    plt.scatter(g["angle_CHOw"], g["angle_HOwA"], s=5, alpha=0.7, rasterized=True)
    plt.xlabel("C2-H-Ow Angle (deg)")
    plt.ylabel("H-Ow-A Angle (deg)")
    plt.title(pair)
    plt.tight_layout()
    plt.savefig(f"{pair}_angle_relay.pdf")
    plt.close()

    plt.figure(figsize=(5,4))
    plt.hist2d(g["angle_CHOw"], g["angle_HOwA"], bins=50, density=True, cmap=cmap)
    plt.xlabel("C2-H-Ow Angle (deg)")
    plt.ylabel("H-Ow-A Angle (deg)")
    plt.title(pair)
    plt.colorbar(label="Probability density")
    plt.tight_layout()
    plt.savefig(f"{pair}_angle_relay_2d.pdf")
    plt.close()

    lifetimes = get_lifetimes(g["frame"])
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
    plt.savefig(f"{pair}_lifetimes_relay.pdf")
    plt.close()

    
# plot occupancy bar chart for all pairs
occupancies = []
occ = df.groupby("pair")["frame"].nunique() / 50000 #nframes CHANGE FOR FULL TRAJ

plt.figure(figsize=(5, 4))
plt.bar(
    occ.index,
    occ*100,
    capsize=3
)
plt.ylabel("Occupancy (%)")
plt.xlabel("A-H Pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("occupancy_relay.pdf")
plt.close()

# Lifetimes for A-H pairs in stacked histograms
groups = list(df.groupby("pair"))

all_lifetimes = []
labels = []

for pair, g in groups:
    lt_frames = get_lifetimes(g["frame"])
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
plt.savefig("lifetimes_histogram_stacked_relay.pdf")
plt.close()

### Water plots ###
# Calculate entropy of waters per A-H pair
entropy_data = []
for pair, g in df.groupby("pair"):
    relay_frames = g["frame"].nunique()
    counts = g["water_resid"].value_counts()
    p_i = counts / counts.sum()
    S = -(p_i * np.log(p_i)).sum()
    # normalized entropy (increases with number of waters observed)
    N = len(p_i)
    if N > 1:
        S_norm = S / np.log(N)
    else:
        S_norm = 0

    entropy_data.append({
        "pair": pair,
        "n_events": counts.sum(),
        "n_waters": len(counts),
        "entropy": S,
        "norm_entropy": S_norm
    })

entropy_df = pd.DataFrame(entropy_data)
entropy_df.to_csv("water_entropy_per_pair.csv", index=False)

# entropy barplot
plt.figure(figsize=(5,4))
entropy_df = entropy_df.sort_values("entropy")
plt.bar(
    entropy_df["pair"],
    entropy_df["entropy"]
)
plt.ylabel("Water entropy S")
plt.xlabel("A-H pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("water_entropy_barplot.pdf")
plt.close()

plt.figure(figsize=(5,4))
entropy_df = entropy_df.sort_values("norm_entropy")
plt.bar(
    entropy_df["pair"],
    entropy_df["norm_entropy"]
)
plt.ylabel("Water entropy S")
plt.xlabel("A-H pair")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("water_norm_entropy_barplot.pdf")
plt.close()

# Probability plots for waters per pair 

for pair, g in df.groupby("pair"):

    counts = g["water_resid"].value_counts()

    p_i = counts / counts.sum()

    prob_df = pd.DataFrame({
        "water_resid": p_i.index,
        "probability": p_i.values
    })

    prob_df.to_csv(
        f"{pair}_water_probabilities.csv",
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
    plt.savefig(f"{pair}_water_probabilities.pdf")
    plt.close()

# Calculate average lifetime of waters per A-H pair 
dt = 2  # ps

for pair, g in df.groupby("pair"):

    water_stats = []

    for water, wdf in g.groupby("water_resid"):

        lifetimes = get_lifetimes(wdf["frame"])

        if len(lifetimes) == 0:
            continue

        water_stats.append({
            "water_resid": water,
            "probability": len(wdf) / len(g),
            "mean_lifetime_ps": np.mean(lifetimes) * dt,
            "median_lifetime_ps": np.median(lifetimes) * dt,
            "max_lifetime_ps": np.max(lifetimes) * dt,
            "n_events": len(lifetimes)
        })

    water_df = pd.DataFrame(water_stats)

    water_df.to_csv(
        f"{pair}_water_statistics.csv",
        index=False
    )

    # plot lifetime and probability
    fig, ax1 = plt.subplots(figsize=(5,4))

    ax1.bar(
        water_df["water_resid"].astype(str),
        water_df["probability"],
    )
    ax1.set_ylabel("Probability")

    ax2 = ax1.twinx()

    ax2.plot(
        water_df["water_resid"].astype(str),
        water_df["mean_lifetime_ps"],
        marker="o",
        color='black'
    )
    ax2.set_ylabel("Mean lifetime (ps)")

    plt.title(pair)
    plt.tight_layout()
    plt.savefig(f"{pair}_water_prob_lifetime.pdf")

# Calculate average occupancy of waters per A-H pair 
