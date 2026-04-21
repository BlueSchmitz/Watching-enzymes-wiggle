#!/bin/bash  

#SBATCH -J BtDERA_plumed_sanity_checks
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16 
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=BtDERA_plumed_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_TASKS=$SLURM_NTASKS

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=BtDERA
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
    if [[ -d "$tmpdir/MD_simulations/projects/$project_dir/outputs" ]]; then
        rsync -av --partial --inplace \
          "$tmpdir/MD_simulations/projects/$project_dir/outputs" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir"
        echo "=== Copy complete at $(date) ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
    # Clean up temporary directory
    rm -rf "$tmpdir"
    echo "=== Temporary directory cleaned up at $(date). ==="
}
trap copy_back_results EXIT

# mkdir outputs
mkdir -p ./outputs/5_sanity_checks ./outputs/6_HREX

### 6 Set up HREX-MD with PLUMED ###
echo "============= Scale topologies with PLUMED =============" 
cd ./outputs/6_HREX
# 4. Scale the Hamiltonian of the selected atoms by the factors 1.00, 0.95, 0.91, 0.87, 0.83, 0.79, 0.76, 0.72, 0.69, 0.66, 0.63, and 0.60
for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
do 
    mkdir -p ./rep${i} # create directory for each replica
    apptainer exec $PLUMED_CONTAINER plumed partial_tempering ${i} < processed_scaled.top  > ./rep${i}/scaled_${i}.top
    echo "Generated scaled_${i}.top in ./rep${i} with scaling factor ${i}."
    cd ./rep${i}
    apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/hrex.mdp -c ../npt_5.gro -p scaled_${i}.top -o topol.tpr -maxwarn 1
    echo "Generated topol.tpr from scaled_${i}.top in ./rep${i}."
    cd ..
done 
: > plumed.dat # empty plumed file
echo "Generated all scaled topology files and tpr files for HREX replicas."