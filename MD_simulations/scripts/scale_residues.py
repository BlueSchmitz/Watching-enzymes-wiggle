#!/usr/bin/env python3
# usage: python scale_residues.py > processed_scaled.top
# possible usage: python scale_residues.py Protein_chain_A > processed_scaled.top
'''Adds "_" to specific atom names in processed.top'''

import sys

# select residues to be scaled
residues_to_scale = list(range(19,26)) + list(range(75,83)) + \
                    list(range(168,179)) + list(range(202,209)) + \
                    list(range(249,260))

chain_filter = sys.argv[1] if len(sys.argv) == 2 else None
section = None
current_chain = None
scale_active = True  # default: scale everything

with open("processed.top") as f:
    for line in f:
        stripped = line.strip()
        
        # Always print comments and blank lines
        if stripped == "" or stripped.startswith(";"):
            print(line, end="")
            continue

        # Section headers
        if stripped.startswith("["):
            section = stripped.lower()
            print(line, end="")
            continue

        # Detect moleculetype name
        if section == "[ moleculetype ]":
            parts = stripped.split()
            if parts:
                current_chain = parts[0]
                scale_active = (chain_filter is None or current_chain == chain_filter)
            print(line, end="")
            continue

        parts = line.split()

        # Modify only inside atoms section and when scaling is active
        if section == "[ atoms ]" and scale_active:
            try:
                resnum = int(parts[2])
                if resnum in residues_to_scale:
                    parts[1] = parts[1] + "_"
            except (IndexError, ValueError):
                pass  # keep line unchanged

        print("  ".join(parts))