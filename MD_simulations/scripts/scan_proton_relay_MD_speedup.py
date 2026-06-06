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
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.transformations import fit_rot_trans, center_in_box, wrap


###############################################################################
# INPUT
###############################################################################

TOP = sys.argv[1]
TRAJ = sys.argv[2]

INTERMEDIATE_RESNAME = "KPS"
C2_NAME = "C1"

DEBUG = True

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

# Reference group for fitting (TIM barrel backbone)
ref = u.select_atoms("resid 1-221 and name CA")

u.trajectory.add_transformations(
    wrap(u.atoms),
    center_in_box(u.select_atoms("protein or resname KPS")),
    fit_rot_trans(u.select_atoms("protein"), ref)
)

lig = u.select_atoms(
    f"resname {INTERMEDIATE_RESNAME}"
)

C2 = lig.select_atoms(
    f"name {C2_NAME}"
)

if len(C2) != 1:
    raise ValueError(
        f"Expected exactly one C2"
    )

bonded_H = lig.select_atoms(
    "name H1 H2"
)

if len(bonded_H) != 2:
    raise ValueError(
        f"Expected exactly two H atoms on C2"
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

for iframe, ts in enumerate(tqdm(u.trajectory, total=n_frames)):

    if DEBUG and iframe % 1000 == 0:
        print(f"\n[Frame {iframe}/{n_frames}]")

    C2_pos = C2.positions[0]

    H_pos = bonded_H.positions               # (2, 3)
    acc_pos = acceptors.positions            # (nA, 3)
    wat_pos = waters.positions               # (nW, 3)

    frame_direct = set()
    frame_relay = set()

    stats_local_direct = {}
    stats_local_relay = {}

    ###########################################################################
    # NEIGHBOR SEARCH (FastNS instead of KDTree)
    ###########################################################################

    if len(acc_pos) == 0:
        continue

    ns_acc = FastNS(ACCEPTOR_SCREEN, acc_pos, acceptors.dimensions)
    results = ns_acc.search(C2_pos.reshape(1, 3))
    nearby_acc = results.get_indices()

    if DEBUG and iframe % 1000 == 0:
        print(f"  Nearby acceptors: {len(nearby_acc)}")

    if len(nearby_acc) == 0:
        continue

    acc_pos_nb = acc_pos[nearby_acc]

    ###########################################################################
    # VECTORIZED H ↔ A DISTANCES
    ###########################################################################

    # (nH, nA)
    d_HA = np.linalg.norm(
        H_pos[:, None, :] - acc_pos_nb[None, :, :],
        axis=-1
    )

    direct_mask = d_HA < DIRECT_DIST_CUTOFF

    if DEBUG and iframe % 1000 == 0:
        print(f"  Direct pairs (H-A < cutoff): {np.sum(direct_mask)}")

    if not np.any(direct_mask):
        continue

    ###########################################################################
    # DIRECT ANGLES (vectorized)
    ###########################################################################

    # broadcast only valid pairs
    H_idx, A_idx = np.where(direct_mask)

    vCH = C2_pos - H_pos[H_idx]
    vHA = acc_pos_nb[A_idx] - H_pos[H_idx]

    cos_theta = np.einsum("ij,ij->i", vCH, vHA) / (
        np.linalg.norm(vCH, axis=1) * np.linalg.norm(vHA, axis=1)
    )

    angles = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

    good = angles > DIRECT_ANGLE_CUTOFF

    if DEBUG and iframe % 1000 == 0:
        print(f"  Direct angle-filtered pairs: {np.sum(good)}")

    H_idx = H_idx[good]
    A_idx = A_idx[good]
    angles = angles[good]

    dists = d_HA[H_idx, A_idx]

    distance_score = (DIRECT_DIST_CUTOFF - dists) / DIRECT_DIST_CUTOFF
    angle_score = (angles - DIRECT_ANGLE_CUTOFF) / (180.0 - DIRECT_ANGLE_CUTOFF)

    quality = 0.5 * (np.clip(distance_score, 0, 1) + np.clip(angle_score, 0, 1))

    ###########################################################################
    # RECORD DIRECT RESULTS
    ###########################################################################

    for h, a, q in zip(H_idx, A_idx, quality):

        acc = acceptors[nearby_acc[a]]

        key = (acc.resname, acc.resid, acc.name, bonded_H[h].name)

        frame_direct.add(key)

        stats[key]["direct_quality_sum"] += q
        stats[key]["direct_frames"] += 1

    ###############################################################################
    # RELAY (CLEAN VECTORISED VERSION)
    ###############################################################################

    if len(wat_pos) == 0:
        continue

    ns_wat = FastNS(ACCEPTOR_SCREEN, wat_pos, waters.dimensions)
    results = ns_wat.search(C2_pos.reshape(1, 3))
    nearby_wat = results.get_indices()

    if DEBUG and iframe % 1000 == 0:
        print(f"  Nearby waters: {len(nearby_wat)}")

    if len(nearby_wat) == 0:
        continue

    wat_pos_nb = wat_pos[nearby_wat]

    # distances
    d_HOw = np.linalg.norm(
        H_pos[:, None, :] - wat_pos_nb[None, :, :],
        axis=-1
    )

    d_OwA = np.linalg.norm(
        wat_pos_nb[:, None, :] - acc_pos_nb[None, :, :],
        axis=-1
    )

    ###############################################################################
    # STEP 1: VALID (H, Ow) PAIRS
    ###############################################################################

    H_idx, W_idx = np.where(d_HOw < HOW_DIST_CUTOFF)

    if DEBUG and iframe % 1000 == 0:
        print(f"  H-W candidate pairs: {len(H_idx)}")

    if len(H_idx) == 0:
        continue

    ###############################################################################
    # STEP 2: EXPAND TO (H, Ow, A)
    ###############################################################################

    H_exp = []
    W_exp = []
    A_exp = []

    for h, w in zip(H_idx, W_idx):

        valid_A = np.where(d_OwA[w] < OWA_DIST_CUTOFF)[0]

        if len(valid_A) == 0:
            continue

        H_exp.extend([h] * len(valid_A))
        W_exp.extend([w] * len(valid_A))
        A_exp.extend(valid_A)

    if DEBUG and iframe % 1000 == 0:
            print(f"  Expanded H-W-A triples: {len(H_exp)}")

    if len(H_exp) == 0:
        continue

    H_exp = np.array(H_exp)
    W_exp = np.array(W_exp)
    A_exp = np.array(A_exp)

    ###############################################################################
    # STEP 3: GEOMETRY
    ###############################################################################

    Hv = H_pos[H_exp]
    Owv = wat_pos_nb[W_exp]
    Av = acc_pos_nb[A_exp]

    vCH = C2_pos - Hv
    vHOw = Owv - Hv
    vOwA = Av - Owv

    def angle(v1, v2):
        return np.degrees(
            np.arccos(
                np.clip(
                    np.einsum("ij,ij->i", v1, v2) /
                    (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)),
                    -1, 1
                )
            )
        )

    CHOw = angle(vCH, vHOw)
    HOwA = angle(vHOw, vOwA)

    mask = (
        (CHOw >= CHOw_ANGLE_CUTOFF) &
        (HOwA >= HOwA_ANGLE_CUTOFF)
    )

    if DEBUG and iframe % 1000 == 0:
        print(f"  Valid relay triples after geometry: {np.sum(mask)}")

    if not np.any(mask):
        continue

    H_exp = H_exp[mask]
    W_exp = W_exp[mask]
    A_exp = A_exp[mask]
    CHOw = CHOw[mask]
    HOwA = HOwA[mask]

    ###############################################################################
    # STEP 4: SCORING
    ###############################################################################

    d1 = d_HOw[H_exp, W_exp]
    d2 = d_OwA[W_exp, A_exp]

    quality = np.mean([
        (HOW_DIST_CUTOFF - d1) / HOW_DIST_CUTOFF,
        (OWA_DIST_CUTOFF - d2) / OWA_DIST_CUTOFF,
        (CHOw - CHOw_ANGLE_CUTOFF) / (180 - CHOw_ANGLE_CUTOFF),
        (HOwA - HOwA_ANGLE_CUTOFF) / (180 - HOwA_ANGLE_CUTOFF),
    ], axis=0)

    ###############################################################################
    # STEP 5: BEST PER (H, A)
    ###############################################################################

    best = {}

    for h, a, q in zip(H_exp, A_exp, quality):

        acc = acceptors[nearby_acc[a]]
        h_atom = H[h]

        key = (acc.resname, acc.resid, acc.name, h_atom.name)

        if key not in best or q > best[key]:
            best[key] = q

    for key, q in best.items():
        stats[key]["relay_quality_sum"] += q
        stats[key]["relay_frames"] += 1


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