#!/usr/bin/env python3
# usage: python get_ali_file.py pdb_ID
'''Generate a sequence alignment file (.seq) from a PDB file using Modeller.'''

import sys
from modeller import *

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <pdb_ID>")
    sys.exit(1)
# Get the sequence of the PDB file, and write to an alignment file
pdb = sys.argv[1]

e = Environ()
m = Model(e, file=pdb)
aln = Alignment(e)
aln.append_model(m, align_codes=pdb)
aln.write(file=pdb+'.seq')