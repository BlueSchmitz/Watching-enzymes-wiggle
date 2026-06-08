#!/usr/bin/env python3
"""
Fast geometric characterization of proton-transfer candidates

Outputs
-------
direct_geometry.csv
relay_geometry.csv
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
DIRECT_CSV = sys.argv[3]
RELAY_CSV = sys.argv[4]

INTERMEDIATE_RESNAME = "KPS"
C2_NAME = "C1"

ACCEPTOR_SCREEN = 6.0
DIRECT_DIST_CUTOFF = 3.0
DIRECT_ANGLE_CUTOFF = 140.0
HOW_DIST_CUTOFF = 3.0
OWA_DIST_CUTOFF = 3.5
CHOw_ANGLE_CUTOFF = 120.0
HOwA_ANGLE_CUTOFF = 100.0

###############################################################################
# HELPERS
###############################################################################

def angle(v1, v2):
    denom = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
    cos = np.einsum("ij,ij->i", v1, v2) / np.clip(denom, 1e-12, None)
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))

###############################################################################
# LOAD SYSTEM
###############################################################################

u = mda.Universe(TOP, TRAJ)

lig = u.select_atoms(f"resname {INTERMEDIATE_RESNAME}")
C2 = lig.select_atoms(f"name {C2_NAME}")
bonded_H = lig.select_atoms("name H1 H2")
waters = u.select_atoms("resname WAT and name O")

if len(C2) != 1:
    raise ValueError("Expected exactly one C2 atom")
if len(bonded_H) != 2:
    raise ValueError("Expected exactly two bonded hydrogens")

###############################################################################
# LOAD CANDIDATES
###############################################################################

direct_df = pd.read_csv(DIRECT_CSV)
relay_df = pd.read_csv(RELAY_CSV)

direct_keys = list(zip(direct_df.resname, direct_df.resid, direct_df.atom, direct_df.hydrogen))
relay_keys = list(zip(relay_df.resname, relay_df.resid, relay_df.atom, relay_df.hydrogen))

###############################################################################
# FAST LOOKUPS (IMPORTANT OPTIMIZATION)
###############################################################################

hydrogen = {h.name: h for h in bonded_H}

acceptor_atoms = {}
for df in [direct_df, relay_df]:
    for r in df.itertuples(index=False):
        key = (r.resname, r.resid, r.atom)
        if key not in acceptor_atoms:
            ag = u.select_atoms(f"resname {r.resname} and resid {r.resid} and name {r.atom}")
            acceptor_atoms[key] = ag[0]

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

    C2_pos = C2.positions[0]
    H_pos = bonded_H.positions
    box = ts.dimensions

    # ============================================================
    # DIRECT (vectorised per candidate batch)
    # ============================================================

    for resname, resid, atom_name, h_name in direct_keys:

        H = hydrogen[h_name]
        A = acceptor_atoms[(resname, resid, atom_name)]

        d = distance_array(
            H.position.reshape(1, 3),
            A.position.reshape(1, 3),
            box=box
        )[0, 0]

        v1 = C2_pos - H.position
        v2 = A.position - H.position

        ang = np.degrees(
            np.arccos(
                np.clip(
                    np.dot(v1, v2) /
                    (np.linalg.norm(v1) * np.linalg.norm(v2)),
                    -1, 1
                )
            )
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
            "within_geometry": bool(
                (d < DIRECT_DIST_CUTOFF) &
                (ang < DIRECT_ANGLE_CUTOFF)
                )
            })

    # ============================================================
    # RELAY (fast prefilter)
    # ============================================================

    ns_wat = FastNS(ACCEPTOR_SCREEN, waters.positions, box)
    nearby = ns_wat.search(C2_pos.reshape(1, 3)).get_pairs()[:, 1]

    if len(nearby) == 0:
        continue

    wat_pos = waters.positions[nearby]

    for resname, resid, atom_name, h_name in relay_keys:

        H = hydrogen[h_name]
        A = acceptor_atoms[(resname, resid, atom_name)]

        d_HOw = distance_array(
            H.position.reshape(1, 3),
            wat_pos,
            box=box
        )[0]

        d_OwA = distance_array(
            wat_pos,
            A.position.reshape(1, 3),
            box=box
        )[:, 0]

        # angles
        vCH = C2_pos - H.position
        vHOw = wat_pos - H.position
        vHA = A.position - wat_pos

        den1 = np.linalg.norm(vCH) * np.linalg.norm(vHOw, axis=1)
        ang_CHOw = np.degrees(
            np.arccos(
                np.clip(np.einsum("ij,j->i", vHOw, vCH) / np.clip(den1, 1e-12, None), -1, 1)
            )
        )

        den2 = np.linalg.norm(vHOw, axis=1) * np.linalg.norm(vHA, axis=1)
        ang_HOwA = np.degrees(
            np.arccos(
                np.clip(np.einsum("ij,ij->i", vHOw, vHA) / np.clip(den2, 1e-12, None), -1, 1)
            )
        )

        # ============================================================
        # GEOMETRY FILTER (THIS IS THE KEY PART)
        # ============================================================
        geom_mask = (
            (d_HOw < HOW_DIST_CUTOFF) &
            (d_OwA < OWA_DIST_CUTOFF) &
            (ang_CHOw >= CHOw_ANGLE_CUTOFF) &
            (ang_HOwA >= HOwA_ANGLE_CUTOFF)
        )

        if not np.any(geom_mask):
            continue

        idx = np.where(geom_mask)[0]

        # ============================================================
        # QUALITY SCORE
        # ============================================================

        dist_score = (
            (HOW_DIST_CUTOFF - d_HOw[idx]) / HOW_DIST_CUTOFF +
            (OWA_DIST_CUTOFF - d_OwA[idx]) / OWA_DIST_CUTOFF
        )

        angle_score = (
            (ang_CHOw[idx] - CHOw_ANGLE_CUTOFF) / (180.0 - CHOw_ANGLE_CUTOFF) +
            (ang_HOwA[idx] - HOwA_ANGLE_CUTOFF) / (180.0 - HOwA_ANGLE_CUTOFF)
        )

        quality = np.clip(dist_score, 0, 1) + np.clip(angle_score, 0, 1)

        best = idx[np.argmax(quality)]

        relay_records.append({
            "frame": ts.frame,
            "time_ps": ts.time,
            "resname": resname,
            "resid": resid,
            "atom": atom_name,
            "hydrogen": h_name,
            "water_resid": waters[nearby[best]].resid,
            "distance_HOw": d_HOw[best],
            "distance_OwA": d_OwA[best],
            "angle_CHOw": ang_CHOw[best],
            "angle_HOwA": ang_HOwA[best],
            "quality": float(np.max(quality)),
        })

###############################################################################
# OUTPUT
###############################################################################

pd.DataFrame(direct_records).to_csv("direct_geometry.csv", index=False)
pd.DataFrame(relay_records).to_csv("relay_geometry.csv", index=False)

print("Done.")