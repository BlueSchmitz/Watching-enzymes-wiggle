#!/usr/bin/env python3
# usage: python model_missing_residues.py alnfile.ali part_sequence full_sequence #models #loop_refinements residue1 residue2
# e.g. python model_missing_residues.py 5EKY_fill.ali 5EKY 5EKY_fill 2 5 A167 A259 (chain and residue number)

'''
Use Modeller to model missing residues in a PDB structure based on its own sequence.
part_sequence and full_sequence should be the terms in the alignment file (e.g. 5EKY and 5EKY_fill).
The alignment file can be generated with get_ali_file.py.
The script also computes the distance between two CA atoms in the generated models to help select the best model.
'''

import sys
from modeller import *
from modeller.automodel import *    # LoopModel is in automodel
import math
import os
import csv

# Handle inputs
if len(sys.argv) != 8:
    print(f"Usage: {sys.argv[0]} <alnfile.ali> <part_sequence> <full_sequence> <#models> <#loop_refinements> <residue1> <residue2>")
    sys.exit(1)

alnfile, part_sequence, full_sequence, models, loop_refinements, res1, res2 = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7]

# parse residues
def parse_residue(res_str):
    """Parse residue like 'A167' → ('A', 167)"""
    chain = res_str[0]
    resnum = int(res_str[1:])
    return chain, resnum

chain1, resnum1 = parse_residue(res1)
chain2, resnum2 = parse_residue(res2)

# Run Modeller 
log.verbose()
env = Environ()
env.io.atom_files_directory = ['.'] # directories for input atom files

# LoopModel to rebuild missing segments and refine the loops
a = LoopModel(env,
              alnfile = alnfile,
              knowns = part_sequence,
              sequence = full_sequence)

# generate this many models
a.starting_model = 1
a.ending_model   = models

# do this many loop refinements
a.loop.starting_model = 1
a.loop.ending_model   = loop_refinements
a.loop.md_level       = refine.fast  # use refine.slow for higher quality

a.make()

# Compute distances (to choose best model)
pdb_folder = '.' # Folder with MODELLER output PDBs
output_csv = f'{part_sequence}_{res1}_{res2}_CA_distance.csv'

# List all PDB files matching MODELLER naming
pdb_files = sorted([f for f in os.listdir(pdb_folder) if f.startswith(full_sequence) and f.endswith('.pdb')])

# Function to get coordinates of CA atom for a given residue number and chain
def get_ca_coordinates(pdb_file, res_num, chain='A'):
    with open(pdb_file, 'r') as f:
        for line in f:
            if line[:4] == 'ATOM':
                # ATOM type (e.g. C-alpha)
                atom_type = line[12:16].strip()
                # AMINO ACID type (e.g. alanine)
                aa_type = line[17:20].strip()
                # residue number
                res_number = int(line[22:26])
                # Protein chain
                chain_id = line[21]
                # return coordinates if it's the CA atom of the specified residue and chain
                # coordinates
                if atom_type == 'CA' and chain_id == chain and res_number == res_num:
                    xcoord = float(line[30:38])
                    ycoord = float(line[38:46])
                    zcoord = float(line[46:54])
                    return (xcoord, ycoord, zcoord)
    return None  # if not found

# Open CSV for writing
with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Model', f'{res1}_to_{res2}_Distance_Angstrom'])
    
    for pdb_file in pdb_files:
        ca1 = get_ca_coordinates(pdb_file, resnum1, chain1)
        ca2 = get_ca_coordinates(pdb_file, resnum2, chain2)
        if ca1 is None or ca2 is None:
            dist = 'NA'
        else:
            dist = math.sqrt((ca1[0]-ca2[0])**2 + (ca1[1]-ca2[1])**2 + (ca1[2]-ca2[2])**2)
        writer.writerow([pdb_file, dist])

print(f'Distances saved to {output_csv}')
