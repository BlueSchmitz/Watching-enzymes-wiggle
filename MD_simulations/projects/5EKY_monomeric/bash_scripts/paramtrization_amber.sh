#!/usr/bin/env bash
set -euo pipefail

##### user input #####
PDB_IN="../../inputs/5EKY_fill.BL00440001.pdb"
PH=7.0
BOX_PADDING_A=15.0     # Angstrom
FF="ff14SB"
WATER="tip3p"
######################

# derived variables
PDB_FILE=$(basename "$PDB_IN")          # 5EKY_fill.BL00440001.pdb
BASE_NOEXT="${PDB_FILE%.pdb}"           # 5EKY_fill.BL00440001
BASENAME="${BASE_NOEXT%%_*}"            # 5EKY
PDB_PH="${BASENAME}_pH${PH}.pdb"        # 5EKY_pH7.0.pdb
PRMTOP="${BASENAME}.prmtop"             # 5EKY.prmtop
INPCRD="${BASENAME}.inpcrd"             # 5EKY.inpcrd

# check dependencies
for cmd in pdb2pqr tleap python; do
    command -v $cmd >/dev/null 2>&1 || {
        echo "ERROR: $cmd not found in PATH"
        exit 1
    }
done

# import parmed
python - <<'EOF'
import parmed
EOF

# 1 Protonation
echo "Protonating at pH ${PH}"
pdb2pqr \
    --ff=AMBER \
    --titration-state-method=propka \
    --with-ph=${PH} \
    --drop-water \
    --keep-chain \
    --pdb-output ${PDB_PH} \
    ${PDB_IN} \
    ${BASENAME}_pH${PH}.pqr

# 2 Strip hydrogens added by pdb2pqr
echo "Stripping extra hydrogens from PDB"

pdb4amber \
  -i ${PDB_PH} \
  -o ${BASENAME}_pH${PH}_noH.pdb \
  --nohyd \
  --dry

# 3 Parametrization with AMBER
echo "Running tleap (${FF} + ${WATER})"

cat > tleap.in <<EOF
source leaprc.protein.${FF}
source leaprc.water.${WATER}

mol = loadpdb ${BASENAME}_pH${PH}_noH.pdb
check mol

solvatebox mol TIP3PBOX ${BOX_PADDING_A} iso
addions mol Na+ 0

saveamberparm mol ${PRMTOP} ${INPCRD}
savepdb mol ${BASENAME}_amber.pdb
quit
EOF

tleap -f tleap.in > tleap.log

# 4 Convert to GROMACS with parmed
echo "Converting to GROMACS with ParmEd"

cat > amber_to_gmx.py <<EOF
import parmed as pmd

amber = pmd.load_file("${PRMTOP}", "${INPCRD}")
amber.save("conf.gro")
amber.save("topol.top")
EOF

python amber_to_gmx.py

echo "Parametrization complete. Generated files:"
echo " - ${PDB_PH} (protonated PDB)"
echo " - ${PRMTOP} (AMBER topology)"
echo " - ${INPCRD} (AMBER coordinates)"
echo " - conf.gro (GROMACS coordinates)"
echo " - topol.top (GROMACS topology)"