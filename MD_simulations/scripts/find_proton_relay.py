#!/usr/bin/env python3
"""
Optimized proton abstraction + water-mediated relay detection pipeline
"""

import sys
import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import distance_array, calc_angles
from scipy.spatial import cKDTree


###############################################################################
# INPUT
###############################################################################

TOP = sys.argv[1]
TRAJ = sys.argv[2]

INTERMEDIATE_RESNAME = "KPS"
C2_NAME = "C1"

C2_SCREEN_CUTOFF = 5.0
FRAME_SAVE_CUTOFF = 4.5

DIRECT_DIST_CUTOFF = 3.0
WATER_DIST_CUTOFF = 3.5


###############################################################################
# LOAD SYSTEM
###############################################################################

u = mda.Universe(TOP, TRAJ)

lig = u.select_atoms(f"resname {INTERMEDIATE_RESNAME}")
C2 = lig.select_atoms(f"name {C2_NAME}")

if len(C2) != 1:
    raise ValueError("C2 not found uniquely")

C2_pos0 = C2.positions[0].copy()

waters = u.select_atoms("resname WAT and name O")

acceptor_sel = (
    "(resname ASP and name OD1 OD2) or "
    "(resname GLU and name OE1 OE2) or "
    "(resname HIS HIE HID HIP HSD HSE HSP and name ND1 NE2) or "
    "(resname TYR and name OH) or "
    "(resname CYS and name SG) or "
    "(resname SER and name OG) or "
    "(resname THR and name OG1) or "
    "(resname WAT and name O)"
)

acceptors = u.select_atoms(acceptor_sel)


###############################################################################
# HYDROGEN SELECTION (cached once)
###############################################################################

hydrogens = lig.select_atoms("name H*")

bonded_H = []
for H in hydrogens:
    if np.linalg.norm(H.position - C2_pos0) < 1.25:
        bonded_H.append(H)

if len(bonded_H) == 0:
    raise ValueError("No bonded hydrogens found")

print(f"Bonded H atoms: {len(bonded_H)}")
print(f"Acceptors: {len(acceptors)}")
print(f"Water oxygens: {len(waters)}")


###############################################################################
# STORAGE
###############################################################################

results = []
per_frame_data = []
water_events = []


###############################################################################
# MAIN LOOP (OPTIMIZED)
###############################################################################

n_frames = len(u.trajectory)

for iframe, ts in enumerate(u.trajectory):

    # --------------------------------------------------
    # CACHE POSITIONS ONCE PER FRAME
    # --------------------------------------------------
    water_pos = waters.positions.copy()
    acc_pos = acceptors.positions.copy()

    # KDTree for fast water proximity queries
    water_tree = cKDTree(water_pos)

    for acc_i, acc in enumerate(acceptors):
        A_pos = acc_pos[acc_i].reshape(1, 3)

        for H in bonded_H:
            H_pos = H.position.reshape(1, 3)

            # --------------------------------------------------
            # CORE GEOMETRY
            # --------------------------------------------------
            dist_HA = distance_array(H_pos, A_pos)[0, 0]

            if dist_HA > FRAME_SAVE_CUTOFF:
                continue

            angle = np.degrees(
                calc_angles(C2_pos0.reshape(1, 3), H_pos, A_pos)[0]
            )

            direct = (dist_HA < DIRECT_DIST_CUTOFF) and (angle > 140)

            # --------------------------------------------------
            # WATER BRIDGE (FAST KD-TREE QUERY)
            # --------------------------------------------------
            bridge_found = False

            # only nearby waters to H
            hw_idx = water_tree.query_ball_point(H_pos[0], WATER_DIST_CUTOFF)

            if hw_idx:
                hw_pos = water_pos[hw_idx]

                d_wa = distance_array(hw_pos, A_pos).flatten()

                if np.any(d_wa < WATER_DIST_CUTOFF):
                    bridge_found = True

                    for i, wi in enumerate(hw_idx):
                        if d_wa[i] < WATER_DIST_CUTOFF:
                            water_events.append({
                                "frame": iframe,
                                "time_ps": ts.time,
                                "hydrogen": H.name,
                                "water_index": waters[wi].index,
                                "acceptor": acc.resname + str(acc.resid),
                                "H_Ow": np.linalg.norm(H_pos[0] - water_pos[wi]),
                                "Ow_A": d_wa[i],
                                "C_H_A_angle": angle
                            })

            # --------------------------------------------------
            # FRAME DATA
            # --------------------------------------------------
            per_frame_data.append({
                "frame": iframe,
                "time_ps": ts.time,
                "hydrogen": H.name,
                "acceptor_resname": acc.resname,
                "acceptor_resid": acc.resid,
                "acceptor_atom": acc.name,
                "H_A_distance": dist_HA,
                "CHA_angle": angle,
                "direct_contact": direct,
                "water_bridge": bridge_found
            })

            # --------------------------------------------------
            # SCORE TRACKING
            # --------------------------------------------------
            if direct:
                results.append((acc.resname, acc.resid, H.name, 1.0, 0.0, 1.0))
            elif bridge_found:
                results.append((acc.resname, acc.resid, H.name, 0.0, 1.0, 0.5))


###############################################################################
# OUTPUT
###############################################################################

df = pd.DataFrame(results, columns=[
    "acceptor", "resid", "H",
    "direct_occupancy", "water_occupancy", "relay_score"
])

per_frame_df = pd.DataFrame(per_frame_data)
water_df = pd.DataFrame(water_events)

df = df.groupby(["acceptor", "resid", "H"], as_index=False).mean()
df = df.sort_values("relay_score", ascending=False)

df.to_csv("relay_candidates.csv", index=False)
per_frame_df.to_csv("proton_abstraction_per_frame.csv", index=False)
water_df.to_csv("water_bridge_events.csv", index=False)

print("\nTop relay candidates:\n")
print(df.head(20))