#!/usr/bin/env python3

import pandas as pd
import MDAnalysis as mda
import os
import numpy as np

# ---------------------------------------------------------
# Input files
# ---------------------------------------------------------

tpr = "../rep1.00/topol.tpr"
xtc = "../rep1.00/traj_comp.xtc"

medoid_table = "medoids_curated.csv"

# ---------------------------------------------------------
# Load trajectory
# ---------------------------------------------------------

u = mda.Universe(tpr, xtc)

# ---------------------------------------------------------
# Read medoid table
# ---------------------------------------------------------

df = pd.read_csv(medoid_table)

# ---------------------------------------------------------
# Precompute trajectory times
# ---------------------------------------------------------

print("Reading trajectory times...")

times = np.array([ts.time for ts in u.trajectory])

print(f"Loaded {len(times)} frames")

# ---------------------------------------------------------
# Output directory
# ---------------------------------------------------------

os.makedirs("medoid_structures", exist_ok=True)

# ---------------------------------------------------------
# Extract medoids
# ---------------------------------------------------------

for _, row in df.iterrows():

    method = row["method"]
    selection = row["selection"]
    cutoff = row["cutoff"]
    cluster = row["cluster"]

    target_time = row["time_ps"]

    # -----------------------------------------------------
    # Find closest frame to target time
    # -----------------------------------------------------

    frame_idx = np.argmin(
        np.abs(times - target_time)
    )

    # Move trajectory to selected frame
    u.trajectory[frame_idx]

    # -----------------------------------------------------
    # Output filename
    # -----------------------------------------------------

    outfile = (
        f"medoid_structures/"
        f"{method}_"
        f"{selection}_"
        f"{cutoff}_"
        f"cluster{cluster}.pdb"
    )

    # -----------------------------------------------------
    # Write structure
    # -----------------------------------------------------

    u.select_atoms("protein").write(outfile)

    actual_time = u.trajectory.time

    print(
        f"Wrote {outfile} "
        f"(target={target_time:.2f} ps, "
        f"actual={actual_time:.2f} ps, "
        f"frame={frame_idx})"
    )

print("Done.")