#!/bin/bash

#SBATCH -J 1JCLm_HREX
#SBATCH -t 48:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 120
#SBATCH --cpus-per-task=1
#SBATCH --gpus=0
#SBATCH --output=./1JCLm_HREX_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

# Exit immediately on errors
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=1JCLm
output_dir=6_HREX
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

### Run HREX ###
apptainer exec $GROMACS_CONTAINER mpirun -np 120 \
    gmx_mpi mdrun \
    -multidir rep* \
    -replex 2000 \
    -plumed ../plumed.dat \
    -cpt 15 \
    -ntomp 1 \
    -hrex \
    -maxh 46.5 \
    -cpi state.cpt

echo "mdrun exited at $(date)"