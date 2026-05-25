#!/bin/bash

#SBATCH -J EcDERA_dump_closed
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task=2
#SBATCH --gpus=0
#SBATCH --output=./EcDERA_dump_closed_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

# Exit immediately on errors
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=EcDERA
output_dir=6_HREX/analysis
pH=7

export GMXLIB=$TMPDIR/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts
pdb=$TMPDIR/MD_simulations/projects/$project_dir/inputs/*.pdb

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_TASKS=$SLURM_NTASKS

### Copy project to scratch ###
echo "=== Copying project to scratch ==="
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/ "$TMPDIR"
cd $TMPDIR/MD_simulations/projects/$project_dir

### Load modules ###
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module list

### Prepare output directory ###
mkdir -p ./outputs/$output_dir
cd ./outputs/$output_dir

# Function to copy back results when error occurs and before the script exits
function copy_back_results {
    set +e +u # Disable exit on error for this function
    echo "=== Copying results back to home at $(date). ==="
    if [[ -d "$TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir" ]]; then
        rsync -av --partial --inplace \
          "$TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/"
        echo "=== Copy complete at $(date) ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
}
trap copy_back_results EXIT

# Downsample and save trajectory
echo "============= Downsizing and Exporting trajectory ============="
#echo -e "q" | apptainer exec $GROMACS_CONTAINER gmx_mpi make_ndx -f topol.tpr -o index.ndx # make index file with default groups
# downsample trajectory 
#echo 1 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s topol.tpr -f traj_comp.xtc -o md_center_mol.xtc -center -pbc mol -ur compact
# fit trajectory to reference (TIM barrel backbone) to remove overall rotation and translation, keep whole system
#echo 4 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s topol.tpr -f md_center_mol.xtc -o md_fit.xtc -fit rot+trans -n index.ndx
#rm md_center_mol.xtc
#echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s topol.tpr -n index.ndx -o md_1000.xtc -dt 1000
#echo "Trajectory saved as md_1000.xtc"

# extract closed conformations
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_1.pdb -dump 1
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_2.pdb -dump 2
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_3.pdb -dump 3
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_4.pdb -dump 4
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_5.pdb -dump 5
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_6.pdb -dump 6
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_7.pdb -dump 7
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_8.pdb -dump 8
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_9.pdb -dump 9
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed.xtc -s ../rep1.00/topol.tpr -n index.ndx -o closed_10.pdb -dump 10

# extract cluster conf
echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s ../rep1.00/topol.tpr -n index.ndx -o cluster6.pdb -dump 22100