#!/bin/bash  

#SBATCH -J umbrella_1JCLm_2windows
#SBATCH -t 04:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=umbrella_2windows_%j.out

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=1JCLm
output_dir=8_umbrella2
pH=7

export GMXLIB=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/inputs/*.pdb

# Create temporary directory on scratch for this job
tmpdir=$(mktemp -d /gpfs/scratch1/shared/rleveson/Blue/tmp.XXXXXX)
mkdir -p "$tmpdir/MD_simulations/projects/"

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module list

### Copy project to scratch ###
echo "=== Copying project to scratch ==="
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir "$tmpdir/MD_simulations/projects/"
cd $tmpdir/MD_simulations/projects/$project_dir

# Function to copy back results when error occurs and before the script exits
function copy_back_results {
    set +e +u # Disable exit on error for this function
    echo "=== Copying results back to home at $(date). ==="
    if [[ -d "$tmpdir/MD_simulations/projects/$project_dir/outputs/$output_dir" ]]; then
        rsync -av --partial --inplace \
          "$tmpdir/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/"
        echo "=== Copy complete at $(date) ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
    # Clean up temporary directory
    rm -rf "$tmpdir"
}
trap copy_back_results EXIT

mkdir -p ./outputs/$output_dir
cd ./outputs/$output_dir/

# 4 Run umbrella sampling windows

cd ./COM_1.084
apptainer exec $GROMACS_CONTAINER mpirun -np 16 gmx_mpi mdrun -deffnm umbrella173 -v -px umbrella173_pullx -pf umbrella173_pullf -cpi umbrella173.cpt

cd ../COM_1.151
apptainer exec $GROMACS_CONTAINER mpirun -np 16 gmx_mpi mdrun -deffnm umbrella146 -v -px umbrella146_pullx -pf umbrella146_pullf -cpi umbrella146.cpt