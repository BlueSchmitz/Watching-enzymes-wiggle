#!/usr/bin/env python3

"""
Characterize direct and relay proton-transfer candidates by extracting
distances and angles per frame for later statistical and persistence analysis.
"""

import sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm

from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.lib.distances import distance_array

###############################################################################

# INPUT

###############################################################################

TOP = sys.argv[1]
TRAJ = sys.argv[2]
direct_candidates_csv = sys.argv[3]
relay_candidates_csv = sys.argv[4]

INTERMEDIATE_RESNAME = "KPS"
C2_NAME = "C1"

ACCEPTOR_SCREEN = 6.0

HOW_DIST_CUTOFF = 3.0
OWA_DIST_CUTOFF = 3.5
CHOw_ANGLE_CUTOFF = 120.0
HOwA_ANGLE_CUTOFF = 100.0

###############################################################################

# HELPER

###############################################################################

def angle(v1, v2):
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom < 1e-12:
        return np.nan
    return np.degrees(
        np.arccos(
            np.clip(np.dot(v1, v2) / denom, -1.0, 1.0)
        )
    )


###############################################################################

# LOAD SYSTEM

###############################################################################

u = mda.Universe(TOP, TRAJ)

lig = u.select_atoms(f"resname {INTERMEDIATE_RESNAME}")
C2 = lig.select_atoms(f"name {C2_NAME}")

if len(C2) != 1:
    raise ValueError("Expected exactly one C2 atom")

bonded_H = lig.select_atoms("name H1 H2")

if len(bonded_H) != 2:
    raise ValueError("Expected exactly two bonded hydrogens")

waters = u.select_atoms("resname WAT and name O")

###############################################################################

# LOAD CANDIDATES

###############################################################################

direct_df = pd.read_csv(direct_candidates_csv)
relay_df = pd.read_csv(relay_candidates_csv)

direct_candidates = [
    (r.resname, r.resid, r.atom, r.hydrogen)
    for _, r in direct_df.iterrows()
]

relay_candidates = [
    (r.resname, r.resid, r.atom, r.hydrogen)
    for _, r in relay_df.iterrows()
]

# hydrogen lookup (fast)

hydrogen_lookup = {h.name: h for h in bonded_H}

# acceptor lookup (direct + relay)

acceptor_lookup = {}

for df in [direct_df, relay_df]:
    for _, r in df.iterrows():

        key = (r.resname, r.resid, r.atom)

        if key in acceptor_lookup:
            continue

        ag = u.select_atoms(
            f"resname {r.resname} and resid {r.resid} and name {r.atom}"
        )

        if len(ag) != 1:
            raise ValueError(f"Could not select atom: {key}")

        acceptor_lookup[key] = ag[0]

# Map: resid → indices of H atoms (in waters only)
water_H_indices = {}

for res in u.select_atoms("resname WAT").residues:
    H_atoms = res.atoms.select_atoms("name H*")
    water_H_indices[res.resid] = H_atoms.indices

###############################################################################

# STORAGE

###############################################################################

direct_records = []
relay_records = []

###############################################################################

# MAIN LOOP

###############################################################################

n_frames = len(u.trajectory)

for ts in tqdm(u.trajectory, total=n_frames):

    all_positions = ts.positions
    C2_pos = C2.positions[0]
    H_pos = bonded_H.positions

    # -----------------------------
    # DIRECT
    # -----------------------------
    for r in direct_candidates:

        resname, resid, atom_name, h_name = r

        acc = acceptor_lookup[(resname, resid, atom_name)]
        H = hydrogen_lookup[h_name]

        d = distance_array(
            H.position.reshape(1, 3),
            acc.position.reshape(1, 3),
            box=ts.dimensions
        )[0, 0]

        ang = angle(
            C2_pos - H.position,
            acc.position - H.position
        )

        direct_records.append({
            "frame": ts.frame,
            "time_ps": ts.time,
            "resname": resname,
            "resid": resid,
            "atom": atom_name,
            "hydrogen": h_name,
            "distance_HA": d,
            "angle_CHA": ang,
        })

    # -----------------------------
    # RELAY (interaction-based storage)
    # -----------------------------

    ns_wat = FastNS(ACCEPTOR_SCREEN, waters.positions, waters.dimensions)
    results = ns_wat.search(C2_pos.reshape(1, 3))
    nearby_wat = results.get_pairs()[:, 1]

    if len(nearby_wat) == 0:
        continue

    wat_pos_nb = waters.positions[nearby_wat]

    for r in relay_candidates:

        resname, resid, atom_name, h_name = r

        acc = acceptor_lookup[(resname, resid, atom_name)]
        H = hydrogen_lookup[h_name]

        # distances: ALL waters in shell
        d_HOw = distance_array(
            H.position.reshape(1, 3),
            wat_pos_nb,
            box=ts.dimensions
        )[0]

        d_OwA = distance_array(
            wat_pos_nb,
            acc.position.reshape(1, 3),
            box=ts.dimensions
        )[:, 0]

        for w in range(len(nearby_wat)):

            water_O = waters[nearby_wat[w]]
            Ow_pos = water_O.position

            ang_CHOw = angle(C2_pos - H.position, Ow_pos - H.position)
            ang_HOwA = angle(H.position - Ow_pos, acc.position - Ow_pos)

            # water hydrogens (optional geometry refinement)
            H_idx = water_H_indices[water_O.resid]
            Hw_pos = all_positions[H_idx]

            ang_hw1 = np.nan
            ang_hw2 = np.nan
            ang_hw_max = np.nan

            if len(Hw_pos) >= 2:
                ang_hw1 = angle(Hw_pos[0] - Ow_pos, acc.position - Ow_pos)
                ang_hw2 = angle(Hw_pos[1] - Ow_pos, acc.position - Ow_pos)
                ang_hw_max = max(ang_hw1, ang_hw2)

            # compute “reactivity score” (NOT a filter)
            score = np.nanmean([
                np.clip((HOW_DIST_CUTOFF - d_HOw[w]) / HOW_DIST_CUTOFF, 0, 1),
                np.clip((OWA_DIST_CUTOFF - d_OwA[w]) / OWA_DIST_CUTOFF, 0, 1),
                np.clip((ang_CHOw - CHOw_ANGLE_CUTOFF) / (180 - CHOw_ANGLE_CUTOFF), 0, 1),
                np.clip((ang_HOwA - HOwA_ANGLE_CUTOFF) / (180 - HOwA_ANGLE_CUTOFF), 0, 1),
            ])

            relay_records.append({
                "frame": ts.frame,
                "time_ps": ts.time,

                "resname": resname,
                "resid": resid,
                "atom": atom_name,
                "hydrogen": h_name,

                "water_resid": water_O.resid,

                "distance_HOw": d_HOw[w],
                "distance_OwA": d_OwA[w],

                "angle_CHOw": ang_CHOw,
                "angle_HOwA": ang_HOwA,

                "angle_Hw1OwA": ang_hw1,
                "angle_Hw2OwA": ang_hw2,
                "angle_HwOwA_max": ang_hw_max,

                # NEW:
                "relay_score": score,
                "within_relay_geometry": (
                    (d_HOw[w] < HOW_DIST_CUTOFF) &
                    (d_OwA[w] < OWA_DIST_CUTOFF) &
                    (ang_CHOw >= CHOw_ANGLE_CUTOFF) &
                    (ang_HOwA >= HOwA_ANGLE_CUTOFF)
                )
            })

###############################################################################

# OUTPUT

###############################################################################

pd.DataFrame(direct_records).to_csv("direct_geometry.csv", index=False)
pd.DataFrame(relay_records).to_csv("relay_geometry.csv", index=False)

print("Done.")
