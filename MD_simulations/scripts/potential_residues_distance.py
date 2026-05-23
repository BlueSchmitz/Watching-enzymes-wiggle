import MDAnalysis as mda
import numpy as np
import pandas as pd

from MDAnalysis.analysis.distances import distance_array

# ---------------- LOAD ----------------
u = mda.Universe("../rep1.00/topol.tpr", "md_fit.xtc")

# adjust selection depending on your system
protein = u.select_atoms("protein")
k166 = u.select_atoms("resid 166 and protein")

cutoff = 6.0  # Å

# ---------------- STORAGE ----------------
contacts = []

# ---------------- TRAJECTORY LOOP ----------------
for ts in u.trajectory:
    # distance matrix: protein atoms vs K166 atoms
    dists = distance_array(protein.positions, k166.positions)

    # atoms within cutoff
    close_atoms = np.where(dists < cutoff)

    close_resids = np.unique(protein.atoms[close_atoms[0]].resids)

    for resid in close_resids:
        contacts.append((ts.frame, resid))

# ---------------- SUMMARY ----------------
df = pd.DataFrame(contacts, columns=["frame", "resid"])

summary = df.groupby("resid").size().reset_index(name="counts")
summary["occupancy"] = summary["counts"] / len(u.trajectory)

summary = summary.sort_values("occupancy", ascending=False)

print(summary)

summary.to_csv("K166_contact_residues.csv", index=False)