from modeller import *
from modeller.automodel import *    # LoopModel is in automodel
import math

log.verbose()
env = Environ()

# directories for input atom files ('.' is current working dir)
env.io.atom_files_directory = ['.']

# LoopModel to rebuild missing segments and refine the loops
a = LoopModel(env,
              alnfile = 'alignment_5EKY.ali',
              knowns = '5EKY',
              sequence = '5EKY_fill')

# generate this many models
a.starting_model = 1
a.ending_model   = 2   # generate 2 models (should be good enough due to homology)

# loop refinement: refine loops for each model in the range [1..5]
a.loop.starting_model = 1
a.loop.ending_model   = 5
a.loop.md_level       = refine.fast  # use refine.slow for higher quality

a.make()

# script to print K167−Y259 distance (to choose best model)
import os
import math
import csv

# Folder with MODELLER output PDBs
pdb_folder = '.'  # current folder
output_csv = 'k167_y259_distances.csv'

# List all PDB files matching MODELLER naming
pdb_files = sorted([f for f in os.listdir(pdb_folder) if f.startswith('5EKY_fill') and f.endswith('.pdb')])

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
                chain = line[21]
                # coordinates
                xcoord = float(line[30:38])
                ycoord = float(line[38:46])
                zcoord = float(line[46:54])

                # return coordinates if it's the CA atom of the specified residue and chain
                if atom_type == 'CA' and res_num == res_number and chain == chain:
                    return (xcoord, ycoord, zcoord)
    return None  # if not found

# Open CSV for writing
with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Model', 'K167_CA_to_Y259_CA_Distance_Angstrom'])
    
    for pdb_file in pdb_files:
        ca1 = get_ca_coordinates(pdb_file, 167, 'A')
        ca2 = get_ca_coordinates(pdb_file, 259, 'A')
        if ca1 is None or ca2 is None:
            dist = 'NA'
        else:
            dist = math.sqrt((ca1[0]-ca2[0])**2 + (ca1[1]-ca2[1])**2 + (ca1[2]-ca2[2])**2)
        writer.writerow([pdb_file, dist])

print(f'Distances saved to {output_csv}')
