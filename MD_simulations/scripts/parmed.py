import MD_simulations.scripts.run_parmed as pmd

prmtop = "1JCL_pH7_KPS.prmtop"
inpcrd = "1JCL_pH7_KPS.inpcrd"

structure = pmd.load_file(prmtop, inpcrd)

# sanity check
print(structure)

# write outputs
structure.save("1JCL_pH7_KPS.top", format="gromacs")
structure.save("1JCL_pH7_KPS.gro", format="gro")