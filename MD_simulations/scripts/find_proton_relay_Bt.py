#!/usr/bin/env python3
"""
Proton abstraction + water-mediated relay detection pipeline
"""

import sys
import numpy as np
import pandas as pd
import mdanalysis as mda

from collections import defaultdict
from MDAnalysis.lib.distances import distance_array, calc_angles


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
ANGLE_CUTOFF = 130.0


###############################################################################
# LOAD SYSTEM
###############################################################################

u = mda.Universe(TOP, TRAJ)

C2 = u.select_atoms(f"resname {INTERMEDIATE_RESNAME} and name {C2_NAME}")

lig = u.select_atoms(f"resname {INTERMEDIATE_RESNAME}")
waters = u.select_atoms("resname WAT and name O")


if len(C2) != 1:
    raise ValueError("C2 not found uniquely")


###############################################################################
# HYDROGENS ON C2
###############################################################################

u.trajectory[0]

hydrogens = lig.select_atoms("name H*")

bonded_H = []
for H in hydrogens:
    if np.linalg.norm(H.position - C2.positions[0]) < 1.25:
        bonded_H.append(H)


###############################################################################
# ACCEPTORS
###############################################################################

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
# STORAGE
###############################################################################

results = []
per_frame_data = []

water_events = []


###############################################################################
# MAIN LOOP
###############################################################################

for acc in acceptors:

    for H in bonded_H:

        direct_frames = []
        water_frames = []

        for iframe, ts in enumerate(u.trajectory):

            H_pos = H.position.reshape(1, 3)
            C_pos = C2.positions
            A_pos = acc.position.reshape(1, 3)

            dist_HA = distance_array(H_pos, A_pos)[0, 0]

            if dist_HA > FRAME_SAVE_CUTOFF:
                continue

            angle = np.degrees(
                calc_angles(C_pos, H_pos, A_pos)[0]
            )

            # --------------------------------------------------
            # DIRECT INTERACTION
            # --------------------------------------------------
            direct = (dist_HA < DIRECT_DIST_CUTOFF) and (angle > 140)

            # --------------------------------------------------
            # WATER BRIDGE DETECTION
            # --------------------------------------------------
            bridge_found = False

            for w in waters:

                w_pos = w.position.reshape(1, 3)

                d_hw = distance_array(H_pos, w_pos)[0, 0]
                if d_hw > WATER_DIST_CUTOFF:
                    continue

                d_wa = distance_array(w_pos, A_pos)[0, 0]
                if d_wa > WATER_DIST_CUTOFF:
                    continue

                bridge_found = True

                water_events.append({
                    "frame": iframe,
                    "time_ps": ts.time,
                    "hydrogen": H.name,
                    "water_index": w.index,
                    "acceptor": acc.resname + str(acc.resid),
                    "H_Ow": d_hw,
                    "Ow_A": d_wa,
                    "C_H_A_angle": angle
                })

            # --------------------------------------------------
            # FRAME OUTPUT
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

            if direct:
                direct_frames.append(iframe)

            if bridge_found:
                water_frames.append(iframe)

        # ------------------------------------------------------
        # SCORE DIRECT RELAY
        # ------------------------------------------------------
        occupancy = len(direct_frames) / len(u.trajectory)

        score = occupancy + 0.5 * (len(water_frames) / len(u.trajectory))

        results.append({
            "acceptor": acc.resname,
            "resid": acc.resid,
            "H": H.name,
            "direct_occupancy": occupancy,
            "water_occupancy": len(water_frames) / len(u.trajectory),
            "relay_score": score
        })


###############################################################################
# OUTPUT
###############################################################################

df = pd.DataFrame(results).sort_values("relay_score", ascending=False)
per_frame_df = pd.DataFrame(per_frame_data)
water_df = pd.DataFrame(water_events)


df.to_csv("relay_candidates.csv", index=False)
per_frame_df.to_csv("proton_abstraction_per_frame.csv", index=False)
water_df.to_csv("water_bridge_events.csv", index=False)

print("\nTop relay candidates:\n")
print(df.head(20))