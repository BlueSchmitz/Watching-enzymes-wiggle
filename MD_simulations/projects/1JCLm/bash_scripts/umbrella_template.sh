#!/bin/bash

### Project-specific settings ###
project_dir=1JCLm
output_dir=8_umbrella2
pH=7

export GMXLIB=$TMPDIR/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts
pdb=$TMPDIR/MD_simulations/projects/$project_dir/inputs/*.pdb

# Short equilibration
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/umbrella_npt.mdp -c ../confXXX.gro -r ../confXXX.gro -p ../topol.top -n ../index.ndx -o nptXXX.tpr
if [ $? -ne 0 ]; then
    echo "FAIL: grompp equilibration failed" > FAIL
    exit 1
fi

apptainer exec $GROMACS_CONTAINER mpirun -np 8 gmx_mpi mdrun -deffnm nptXXX -v
if [ $? -ne 0 ]; then
    echo "FAIL: mdrun equilibration failed" > FAIL
    exit 1
fi

# Umbrella run
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/umbrella_md.mdp -c nptXXX.gro -r nptXXX.gro -t nptXXX.cpt -p ../topol.top -n ../index.ndx -o umbrellaXXX.tpr
if [ $? -ne 0 ]; then
    echo "FAIL: grompp umbrella failed" > FAIL
    exit 1
fi

apptainer exec $GROMACS_CONTAINER mpirun -np 8 gmx_mpi mdrun -deffnm umbrellaXXX -v
if [ $? -ne 0 ]; then
    echo "FAIL: mdrun umbrella failed" > FAIL
    exit 1
fi
