import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------- CONFIG ----------------
INPUT_FILE = "mutagenesis_activity_Tm.xlsx"

within_group_spacing = 0.6
gap_between_groups = 1

np.random.seed(0)  # only used for point jitter

# ---------------- LOAD DATA ----------------
df = pd.read_excel(INPUT_FILE)

for col in ["Specific_activity", "Tm"]:
    df[col] = (
        df[col]
        .astype(str)
        .str.replace(",", ".")
        .astype(float)
    )

wildtypes = df["Wildtype"].drop_duplicates().tolist()

# Sample one color per wildtype from viridis
wt_colors = cm.viridis(np.linspace(0, 1, len(wildtypes)))


# ---------------- PLOT FUNCTION ----------------
def make_boxplot(value_col, ylabel, output_file):

    fig, ax = plt.subplots(figsize=(8, 6))

    box_data = []
    positions = []
    labels = []
    box_colors = []

    x = 0

    for i, wt in enumerate(wildtypes):

        sub = df[df["Wildtype"] == wt]
        enzymes = sub["Enzyme"].drop_duplicates().tolist()

        wt_color = wt_colors[i]

        for enz in enzymes:

            vals = sub[sub["Enzyme"] == enz][value_col].values

            box_data.append(vals)
            positions.append(x)
            labels.append(enz)

            # Same color for all variants of a WT
            box_colors.append(wt_color)

            x += within_group_spacing

        x += gap_between_groups

    # ---------------- BOXPLOT ----------------
    bp = ax.boxplot(
        box_data,
        positions=positions,
        patch_artist=True,
        showfliers=False,
        widths=0.5
    )

    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)

    for element in ["whiskers", "caps", "medians"]:
        plt.setp(bp[element], color="black")

    # ---------------- OVERLAY DATAPOINTS ----------------
    for pos, vals in zip(positions, box_data):
        jitter = np.random.normal(pos, 0.05, size=len(vals))

        ax.scatter(
            jitter,
            vals,
            facecolors="none",
            edgecolors="black",
            s=30,
            linewidth=1,
            alpha=0.8,
            zorder=3
        )

    # ---------------- FORMATTING ----------------
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right")

    ax.set_ylabel(ylabel)
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_file}")


# ---------------- GENERATE BOTH PLOTS ----------------
make_boxplot(
    value_col="Specific_activity",
    ylabel="Specific activity (U/mg)",
    output_file="enzyme_activity_grouped_boxplot.pdf"
)

make_boxplot(
    value_col="Tm",
    ylabel="Tm (°C)",
    output_file="enzyme_Tm_grouped_boxplot.pdf"
)