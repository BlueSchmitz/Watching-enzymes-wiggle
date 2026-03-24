#!/usr/bin/env python3
# usage: python scale_residues.py > processed_scaled.top
# possible usage: python scale_residues.py Protein_chain_A > processed_scaled.top
'''Adds "_" to specific atom names in processed.top If a chain identifier is given as an additional
argument, only that chain will be processed. Otherwise the residue count of chain A will be used to
calculate the global residue numbers of all chains.'''

import sys

# Residues to scale (per chain, local numbering;
# in case of multiple chains the residue number of chain A will be added automatically)
residues_to_scale = set().union(
    range(1,6),
    range(18,25),
    range(69,77),
    range(167,178),
    range(196,203),
    range(220,230)
)

chain_filter = sys.argv[1] if len(sys.argv) == 2 else None
section = None
current_chain = None
scale_active = True  # default: scale everything (only no if additional argument is given)
first_resnum = None  # first residue number of current chain (for local numbering)

with open("processed.top") as f:
    for line in f:
        stripped = line.strip()

        # Preserve blank lines and comments
        if stripped == "" or stripped.startswith(";"):
            print(line, end="")
            continue

        # Track section headers
        if stripped.startswith("["):
            section = stripped.lower()
            print(line, end="")
            continue

        # Detect chain name from moleculetype
        if section == "[ moleculetype ]":
            current_chain = stripped.split()[0]
            scale_active = (chain_filter is None or current_chain == chain_filter)
            first_resnum = None  # reset for each chain
            print(line, end="")
            continue

        parts = line.split()

        # Modify only atoms in selected chains
        if section == "[ atoms ]" and scale_active:
            try:
                global_res = int(parts[2])

                # First residue in this chain defines local index 1
                if first_resnum is None:
                    first_resnum = global_res

                local_res = global_res - first_resnum + 1

                if local_res in residues_to_scale:
                    parts[1] += "_"

            except (IndexError, ValueError):
                pass  # leave line unchanged

        print("  ".join(parts))