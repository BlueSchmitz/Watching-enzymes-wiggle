#!/usr/bin/env python3
# usage: compare_protonation.py pdb1 pdb2 differences.txt pdb3 renamed_pdb.pdb
'''
Compare two PDB files (pdb1 and pdb2) for differences in hydrogen counts per residue,
and rename residues in a third PDB file (pdb3) based on these differences.
The output is saved in renamed_pdb.pdb and can be used by GROMACS.
'''
import sys

def count_hydrogens(pdb_file):
    """Count hydrogens per residue from a PDB file."""
    counts = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom = line[12:16].strip()
            resname = line[17:20].strip()
            resid = line[22:26].strip()

            if atom.startswith("H"):  # only hydrogens
                key = (resname, resid)
                counts[key] = counts.get(key, 0) + 1
    return counts


def compare_hydrogens(pdb1, pdb2, outfile):
    """Write residues with differing H counts and return dict of differences."""
    counts1 = count_hydrogens(pdb1)
    counts2 = count_hydrogens(pdb2)

    diffs = {}
    with open(outfile, "w") as out:
        out.write("# resname resid H_count_in_1 H_count_in_2\n")
        for key in sorted(set(counts1.keys()) | set(counts2.keys()), key=lambda x: int(x[1])):
            h1 = counts1.get(key, 0)
            h2 = counts2.get(key, 0)
            if h1 != h2:
                resname, resid = key
                out.write(f"{resname.lower()} {resid} {h1} {h2}\n")
                diffs[key] = (h1, h2)
    return diffs, counts1


def map_protonation(resname, hcount):
    """
    Map residue to protonation state based on hydrogen count.
    Names are based on GROMACS conventions.
    """
    resname = resname.upper()
    if resname == "HIS":
        if hcount >= 8:
            return "HIP"  # fully protonated
        else:
            return "HIE"  # default partially protonated
    elif resname == "ASP":
        return "ASH" if hcount > 4 else "ASP"
    elif resname == "GLU":
        return "GLH" if hcount > 5 else "GLU"
    elif resname == "LYS":
        return "LYN" if hcount < 9 else "LYS"
    elif resname == "ARG":
        return "ARG" if hcount >= 11 else "ARN" 
    else:
        return resname


def rewrite_pdb(pdb_in, pdb_out, diffs, counts_ref):
    """Rewrite pdb3 with updated residue names into renamed_pdb."""
    with open(pdb_in, 'r') as f_in, open(pdb_out, 'w') as f_out:
        for line in f_in:
            if not line.startswith("ATOM"):
                f_out.write(line)
                continue

            resname = line[17:20].strip()
            resid = line[22:26].strip()
            key = (resname, resid)

            if key in diffs: # change name in pdb3 if different hydrogen count
                new_resname = map_protonation(resname, counts_ref.get(key, 0))
                line = line[:17] + f"{new_resname:>3}" + line[20:]
            f_out.write(line)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} pdb1 pdb2 differences.txt pdb3 renamed_pdb.pdb")
        sys.exit(1)

    pdb1, pdb2, outfile, pdb3, renamed_pdb = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

    # Step 1: Compare & save differences
    diffs, counts1 = compare_hydrogens(pdb1, pdb2, outfile)

    # Step 2: Apply renaming to pdb3
    rewrite_pdb(pdb3, renamed_pdb, diffs, counts1)
