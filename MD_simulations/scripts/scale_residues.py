#!/usr/bin/env python3
# usage: python scale_residues.py > processed_scaled.top
# possible usage: python scale_residues.py Protein_chain_A > processed_scaled.top
'''Adds "_" to specific atom names in processed.top If a chain identifier is given as an additional
argument, only that chain will be processed. Otherwise the residue count of chain A will be used to
calculate the global residue numbers of all chains.'''

import sys

# Residues to scale (per chain, local numbering;
# in case of multiple chains the residue number of chain A will be added automatically)
residues_to_scale = {
    *range(19,26),
    *range(75,83),
    *range(168,179),
    *range(202,209),
    *range(249,260),
}

chain_filter = sys.argv[1] if len(sys.argv) == 2 else None
section = None
current_chain = None
scale_active = True  # default: scale everything (only no if additional argument is given)
last_resnum = None # to automatically detect residue number of chain A
res_offset = 0  # keeps track of global --> local mapping

with open("processed.top") as f:
    for line in f:
        stripped = line.strip()

        if stripped == "" or stripped.startswith(";"):
            print(line, end="")
            continue

        # Section headers
        if stripped.startswith("["):
            section = stripped.lower()
            print(line, end="")
            continue

        # Molecule name
        if section == "[ moleculetype ]":
            current_chain = stripped.split()[0]
            scale_active = (chain_filter is None or current_chain == chain_filter)
            last_resnum = None
            res_offset = 0
            print(line, end="")
            continue

        parts = line.split()

        # Atoms section
        if section == "[ atoms ]" and scale_active:
            try:
                global_res = int(parts[2])

                # Detect new chain by residue number reset
                if last_resnum is not None and global_res < last_resnum:
                    res_offset += last_resnum

                local_res = global_res - res_offset
                last_resnum = global_res

                if local_res in residues_to_scale:
                    parts[1] += "_"

            except (IndexError, ValueError):
                pass

        print("  ".join(parts))