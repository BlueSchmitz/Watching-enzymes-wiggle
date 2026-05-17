import parmed as pmd

prmtop = "1JCL_KPS.prmtop"
inpcrd = "1JCL_KPS.inpcrd"

structure = pmd.load_file(prmtop, inpcrd)

# sanity check
print(structure)

# write outputs
structure.save("1JCL_KPS.top", format="gromacs")
structure.save("1JCL_KPS.gro", format="gro")