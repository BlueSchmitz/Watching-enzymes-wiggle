#!/usr/bin/env python3
# usage: python scale_residues.py > processed_scaled.top
'''Adds "_" to specific atom names in processed.top'''

# list with all residues named in the paper to be scaled
residues_to_scale = list(range(19,26)) + list(range(75,83)) + \
                    list(range(168,179)) + list(range(202,209)) + \
                    list(range(249,260))

section = None

with open("processed.top") as f:
    for line in f:
        stripped = line.strip()
        # Track section
        if stripped.startswith("["):
            section = stripped.lower()
            print(line, end="")
            continue
        if stripped == "" or stripped.startswith(";"):
            print(line, end="")
            continue
        parts = line.split()
        # Only process atoms section
        if section == "[ atoms ]":
            try:
                resnum = int(parts[2])
            except:
                print(line, end="")
                continue
            if resnum in residues_to_scale:
                parts[1] = parts[1] + "_"
        print("  ".join(parts))
