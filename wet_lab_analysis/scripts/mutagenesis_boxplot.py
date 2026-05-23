import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors


# ---------------- CONFIG ----------------
INPUT_FILE = "mutagenesis_activity.xlsx"
OUTPUT_FILE = "enzyme_activity_grouped_boxplot.pdf"

within_group_spacing = 0.6
gap_between_groups = 1

np.random.seed(0)  # reproducible jitter + color shifts


# ---------------- LOAD DATA ----------------
df = pd.read_excel(INPUT_FILE)

df["Specific_activity"] = (
    df["Specific_activity"]
    .astype(str)
    .str.replace(",", ".")
    .astype(float)
)

wildtypes = sorted(df["Wildtype"].unique())


# ---------------- COLOR HELPERS ----------------
def perturb_color(color, strength=0.5):
    rgb = np.array(color)

    if len(rgb) == 4:
        rgb = rgb[:3]

    noise = np.random.uniform(-strength, strength, 3)

    new = rgb + noise
    return np.clip(new, 0, 1)


# ---------------- PLOT SETUP ----------------
fig, ax = plt.subplots(figsize=(8, 6))

box_data = []
positions = []
labels = []
box_colors = []

x = 0

cmap = cm.viridis


# ---------------- BUILD DATA ----------------
for i, wt in enumerate(wildtypes):
    sub = df[df["Wildtype"] == wt]
    enzymes = sorted(sub["Enzyme"].unique())

    # anchor viridis color per wildtype
    base_color = cmap(i / max(len(wildtypes) - 1, 1))

    for j, enz in enumerate(enzymes):
        vals = sub[sub["Enzyme"] == enz]["Specific_activity"].values

        box_data.append(vals)
        positions.append(x)
        labels.append(enz)

        # mutant-specific slight color variation around WT color
        box_colors.append(perturb_color(base_color, strength=0.20))

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


# ---------------- OVERLAY DATAPOINTS (OUTLINE ONLY) ----------------
for pos, vals in zip(positions, box_data):
    jitter = np.random.normal(pos, 0.05, size=len(vals))

    ax.scatter(
        jitter,
        vals,
        facecolors='none',     
        edgecolors='black',    
        s=30,
        linewidth=1,
        alpha=0.8,
        zorder=3
    )


# ---------------- FORMATTING ----------------
ax.set_xticks(positions)
ax.set_xticklabels(labels, rotation=45, ha="right")

ax.set_ylabel("Specific activity (U/mg)")

ax.grid(False) 


# ---------------- SAVE ----------------
fig.tight_layout()
fig.savefig(OUTPUT_FILE, format="pdf", bbox_inches="tight")

plt.close(fig)

print(f"Saved plot to: {OUTPUT_FILE}")