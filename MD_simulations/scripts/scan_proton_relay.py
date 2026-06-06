#!/usr/bin/env python3

"""
Fast proton abstraction / water-relay screening

Outputs
-------
direct_candidates.csv
relay_candidates.csv

Scoring
--------
direct_score =
    direct_occupancy * mean_direct_quality

relay_score =
    relay_occupancy * mean_relay_quality
"""

import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm

from scipy.spatial import cKDTree
from MDAnalysis.lib.distances import calc_angles


###############################################################################
# INPUT
###############################################################################

TOP = sys.argv[1]
TRAJ = sys.argv[2]

INTERMEDIATE_RESNAME = "KPS"
C2_NAME = "C1"

###############################################################################
# GEOMETRY CUTS
###############################################################################

# Nearby acceptor search
ACCEPTOR_SCREEN = 6.0

# Direct abstraction
DIRECT_DIST_CUTOFF = 3.0
DIRECT_ANGLE_CUTOFF = 140.0

# Water relay
HOW_DIST_CUTOFF = 3.0
OWA_DIST_CUTOFF = 3.5

CHOw_ANGLE_CUTOFF = 120.0
HOwA_ANGLE_CUTOFF = 100.0


###############################################################################
# HELPER
###############################################################################

def clamp01(x):
    return max(0.0, min(1.0, x))


###############################################################################
# LOAD
###############################################################################

u = mda.Universe(TOP, TRAJ)

lig = u.select_atoms(
    f"resname {INTERMEDIATE_RESNAME}"
)

C2 = lig.select_atoms(
    f"name {C2_NAME}"
)

if len(C2) != 1:
    raise ValueError(
        f"Expected exactly one {C2_NAME}"
    )

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
    "(resname THR and name OG1)"
)

acceptors = u.select_atoms(
    acceptor_sel
)

waters = u.select_atoms(
    "resname WAT and name O"
)

###############################################################################
# DONOR H IDENTIFICATION
###############################################################################

hydrogens = lig.select_atoms(
    "name H*"
)

bonded_H = []

C2_pos0 = C2.positions[0]

for H in hydrogens:

    d = np.linalg.norm(
        H.position - C2_pos0
    )

    if d < 1.25:
        bonded_H.append(H)

if len(bonded_H) == 0:
    raise ValueError(
        "No donor hydrogens found near C2"
    )

print()
print(f"Donor H atoms: {len(bonded_H)}")
print(f"Acceptors: {len(acceptors)}")
print(f"Waters: {len(waters)}")
print()

###############################################################################
# STORAGE
###############################################################################

stats = defaultdict(lambda: {
    "direct_frames": 0,
    "relay_frames": 0,
    "direct_quality_sum": 0.0,
    "relay_quality_sum": 0.0,
})

###############################################################################
# MAIN LOOP
###############################################################################

n_frames = len(u.trajectory)

for iframe, ts in enumerate(
    tqdm(u.trajectory, total=n_frames)
):

    # Keep track of which acceptors had direct/relay in this frame for occupancy counting
    frame_direct = set()
    frame_relay = set() 

    ###########################################################################
    # CACHE POSITIONS
    ###########################################################################

    C2_pos = C2.positions[0].reshape(1, 3)

    acc_pos = acceptors.positions.copy()
    water_pos = waters.positions.copy()

    acc_tree = cKDTree(acc_pos)
    water_tree = cKDTree(water_pos)

    ###########################################################################
    # Distance screen
    ###########################################################################

    nearby_acc = acc_tree.query_ball_point(
        C2_pos[0],
        ACCEPTOR_SCREEN
    )

    if not nearby_acc:
        continue

    nearby_waters = water_tree.query_ball_point(
        C2_pos[0],
        ACCEPTOR_SCREEN
    )

    # Loop through nearby acceptors
    for acc_idx in nearby_acc:

        acc = acceptors[acc_idx]

        A_pos = acc_pos[acc_idx].reshape(1, 3)
        
        # Loop through donor hydrogens
        for H in bonded_H:

            H_pos = H.position.reshape(1, 3)

            key = (
                    acc.resname,
                    acc.resid,
                    acc.name,
                    H.name
                )
            

            ###########################################################################
            # Direct transfer 
            ###########################################################################
            
            # Calculate distance and angle for direct abstraction
            dist_HA = np.linalg.norm(
                H_pos[0] - A_pos[0]
            )

            direct_quality = 0.0
            angle_CHA = 0.0

            if dist_HA < DIRECT_DIST_CUTOFF:

                angle_CHA = np.degrees(
                    calc_angles(
                        C2_pos,
                        H_pos,
                        A_pos
                    )[0]
                )

                if angle_CHA > DIRECT_ANGLE_CUTOFF:

                    distance_score = clamp01(
                        (DIRECT_DIST_CUTOFF - dist_HA)
                        / DIRECT_DIST_CUTOFF
                    )

                    angle_score = clamp01(
                        (angle_CHA - DIRECT_ANGLE_CUTOFF)
                        / (180.0 - DIRECT_ANGLE_CUTOFF)
                    )

                    direct_quality = (
                        distance_score +
                        angle_score
                    ) / 2.0

                    frame_direct.add(key)
                    stats[key]["direct_quality_sum"] += direct_quality


            ###################################################################
            # RELAY
            ###################################################################

            bridge_found = False
            best_relay_quality = 0.0
            relay_quality = 0.0

            best_relay_HOw = np.inf
            best_relay_OwA = np.inf
            best_relay_CHOw = 0.0
            best_relay_HOwA = 0.0

            if nearby_waters:
                
                # Loop over all nearby waters and check for relay geometry
                for wi in nearby_waters:

                    Ow_pos = water_pos[wi].reshape(1, 3)

                    dist_HOw = np.linalg.norm(
                        H_pos[0] - Ow_pos[0]
                    )

                    if dist_HOw > HOW_DIST_CUTOFF:
                        continue

                    dist_OwA = np.linalg.norm(
                        Ow_pos[0] - A_pos[0]
                    )

                    if dist_OwA > OWA_DIST_CUTOFF:
                        continue

                    angle_CHOw = np.degrees(
                        calc_angles(
                            C2_pos,
                            H_pos,
                            Ow_pos
                        )[0]
                    )

                    if angle_CHOw < CHOw_ANGLE_CUTOFF:
                        continue

                    angle_HOwA = np.degrees(
                        calc_angles(
                            H_pos,
                            Ow_pos,
                            A_pos
                        )[0]
                    )

                    if angle_HOwA < HOwA_ANGLE_CUTOFF:
                        continue

                    # Score the relay geometry
                    dist_score_HOw = clamp01(
                        (HOW_DIST_CUTOFF - dist_HOw)
                        / HOW_DIST_CUTOFF
                    )

                    dist_score_OwA = clamp01(
                        (OWA_DIST_CUTOFF - dist_OwA)
                        / OWA_DIST_CUTOFF
                    )

                    angle_score_CHOw = clamp01(
                        (angle_CHOw - CHOw_ANGLE_CUTOFF)
                        / (180.0 - CHOw_ANGLE_CUTOFF)
                    )

                    angle_score_HOwA = clamp01(
                        (angle_HOwA - HOwA_ANGLE_CUTOFF)
                        / (180.0 - HOwA_ANGLE_CUTOFF)
                    )

                    current_quality = np.mean([
                        dist_score_HOw,
                        dist_score_OwA,
                        angle_score_CHOw,
                        angle_score_HOwA
                    ])

                    # Update best relay geometry if this one is better
                    if current_quality > best_relay_quality:

                        bridge_found = True

                        best_relay_quality = current_quality

                        best_relay_HOw = dist_HOw
                        best_relay_OwA = dist_OwA
                        best_relay_CHOw = angle_CHOw
                        best_relay_HOwA = angle_HOwA

                if bridge_found:
                    frame_relay.add(key)
                    stats[key]["relay_quality_sum"] += best_relay_quality           


###############################################################################
# CALCULATE SCORE AND BUILD SUMMARY
###############################################################################

rows_direct = []
rows_relay = []

for key, s in stats.items():

    direct_occ = s["direct_frames"] / max(1, n_frames)
    relay_occ = s["relay_frames"] / max(1, n_frames)

    if s["direct_frames"] > 0:
        rows_direct.append({
            "resname": key[0],
            "resid": key[1],
            "atom": key[2],
            "hydrogen": key[3],
            "occupancy": direct_occ,
            "mean_quality": s["direct_quality_sum"] / s["direct_frames"],
            "score": direct_occ * (s["direct_quality_sum"] / s["direct_frames"])
        })

    if s["relay_frames"] > 0:
        rows_relay.append({
            "resname": key[0],
            "resid": key[1],
            "atom": key[2],
            "hydrogen": key[3],
            "occupancy": relay_occ,
            "mean_quality": s["relay_quality_sum"] / s["relay_frames"],
            "score": relay_occ * (s["relay_quality_sum"] / s["relay_frames"])
        })

df_direct = pd.DataFrame(rows_direct).sort_values("score", ascending=False)
df_relay = pd.DataFrame(rows_relay).sort_values("score", ascending=False)

df_direct.to_csv("direct_candidates.csv", index=False)
df_relay.to_csv("relay_candidates.csv", index=False)

print(f"Found {len(df_direct)} direct candidates and {len(df_relay)} relay candidates.")