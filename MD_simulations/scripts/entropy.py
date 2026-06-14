import glob
import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# Load all replicates
# ------------------------------------------------------------------

dfs = []
files = sorted(glob.glob("relay_geometry_*.csv"))

if len(files) == 0:
    raise FileNotFoundError("No relay_geometry_*.csv files found.")

for rep, file in enumerate(files):

    df = pd.read_csv(file)

    df["replicate"] = rep

    # Residue identifier
    df["residue"] = (
        df["resname"] + "_"
        + df["resid"].astype(str) + "_"
        +df["atom"].astype(str)
    )

    dfs.append(df)

print(f"Loaded {len(dfs)} replicate(s)")

# ------------------------------------------------------------------
# Compute residue-level entropy per replicate
# ------------------------------------------------------------------

entropy_data = []

for rep, df in enumerate(dfs):

    for residue, g in df.groupby("residue"):

        counts = g["water_resid"].value_counts()

        p_i = counts / counts.sum()

        # Shannon entropy
        S = -(p_i * np.log(p_i)).sum()

        # Normalized entropy
        N = len(p_i)
        S_norm = S / np.log(N) if N > 1 else 0.0

        entropy_data.append({
            "replicate": rep,
            "residue": residue,
            "entropy": S,
            "norm_entropy": S_norm,
            "n_waters": N
        })

entropy_df = pd.DataFrame(entropy_data)

# ------------------------------------------------------------------
# Aggregate across replicates
# ------------------------------------------------------------------

entropy_summary = (
    entropy_df.groupby("residue")
    .agg(
        entropy_mean=("entropy", "mean"),
        entropy_std=("entropy", "std"),
        norm_entropy_mean=("norm_entropy", "mean"),
        norm_entropy_std=("norm_entropy", "std"),
        replicates=("replicate", "count")
    )
    .reset_index()
)

entropy_summary["entropy_std"] = entropy_summary["entropy_std"].fillna(0.0)
entropy_summary["norm_entropy_std"] = entropy_summary["norm_entropy_std"].fillna(0.0)

entropy_summary = entropy_summary.sort_values(
    "norm_entropy_mean",
    ascending=False
)

# ------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------

print("\nResidue-level normalized entropy:")
print("-" * 80)

for _, row in entropy_summary.iterrows():
    print(
        f"{row['residue']:15s}  "
        f"S_norm = {row['norm_entropy_mean']:.3f} ± "
        f"{row['norm_entropy_std']:.3f}  "
        f"(S = {row['entropy_mean']:.3f} ± "
        f"{row['entropy_std']:.3f})"
    )