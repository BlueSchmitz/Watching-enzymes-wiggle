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

# ---------------- TRAJECTORY LOOP ----------------s
for ts in u.trajectory:

    dists = distance_array(protein.positions, k166.positions)

    close_atom_indices = np.where(dists < cutoff)[0]

    close_atoms = protein.atoms[close_atom_indices]

    # unique residues only
    seen = set()

    for atom in close_atoms:

        resid = atom.resid
        resname = atom.resname

        if resid not in seen and resid != 166:
            seen.add(resid)

            contacts.append({
                "frame": ts.frame,
                "time_ps": ts.time,
                "resid": resid,
                "resname": resname
            })

# ---------------- DATAFRAME ----------------
df = pd.DataFrame(contacts)

# occupancy statistics
summary = (
    df.groupby(["resid", "resname"])
      .size()
      .reset_index(name="n_frames")
)

summary["occupancy"] = summary["n_frames"] / len(u.trajectory)

summary = summary.sort_values(
    "occupancy",
    ascending=False
)

print(summary)

summary.to_csv(
    "K166_nearby_residues.csv",
    index=False
)