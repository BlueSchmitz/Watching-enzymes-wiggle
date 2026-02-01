#!/bin/bash  

#SBATCH -J 1JCLm_HREX
#SBATCH -t 120:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 120
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLm_HREX_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --signal=B:USR1@3600

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
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/ "$TMPDIR"
cd $TMPDIR/MD_simulations/projects/$project_dir

### Trap: stop GROMACS cleanly when walltime approaches ###
function handle_usr1 {
    echo "=== USR1 received: stopping mdrun cleanly at $(date) ==="
    touch STOP_MDRUN
}
trap handle_usr1 USR1

### Load modules ###
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module list

### Prepare output directory ###
mkdir -p ./outputs/$output_dir
cd ./outputs/$output_dir

echo "============= HREX-MD with GROMACS and PLUMED ============="
echo "Starting mdrun at $(date)"

### Run HREX ###
apptainer exec $GROMACS_CONTAINER mpirun -np 120 \
    gmx_mpi mdrun \
    -multidir rep* \
    -replex 2000 \
    -plumed ../plumed.dat \
    -cpt 15 \
    -ntomp 1 \
    -hrex \
    -cpi state.cpt \
    -stop STOP_MDRUN

echo "mdrun finished at $(date)"

### Copy results back ###
echo "============= Copying project outputs back to home ============="
rsync -av --partial --inplace \
    --exclude 'rleveson.*' \
    --exclude '*.out' \
    "$TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
    "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/"

echo "============= Copy complete at $(date) ============="
