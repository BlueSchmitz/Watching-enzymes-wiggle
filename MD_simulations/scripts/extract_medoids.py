#!/usr/bin/env python3

import pandas as pd
import MDAnalysis as mda
import os
import numpy as np

# ---------------------------------------------------------
# Input files
# ---------------------------------------------------------

tpr = "topol.tpr"
xtc = "md_fit.xtc"

medoid_table = "medoids.csv"

# ---------------------------------------------------------
# Load trajectory
# ---------------------------------------------------------

u = mda.Universe(tpr, xtc)

# ---------------------------------------------------------
# Read medoid table
# ---------------------------------------------------------

df = pd.read_csv(medoid_table)

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

    # Find closest frame to target time
    frame_idx = np.argmin(
        np.abs(u.trajectory.times - target_time)
    )

    u.trajectory[frame_idx]

    outfile = (
        f"medoid_structures/"
        f"{method}_"
        f"{selection}_"
        f"{cutoff}_"
        f"cluster{cluster}.pdb"
    )

    u.select_atoms("protein").write(outfile)

    print(
        f"Wrote {outfile} "
        f"(time={target_time} ps, frame={frame_idx})"
    )