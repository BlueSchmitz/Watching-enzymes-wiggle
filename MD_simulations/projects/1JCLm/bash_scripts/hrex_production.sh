#!/bin/bash

#SBATCH -J 1JCLm_HREX
#SBATCH -t 48:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 120
#SBATCH --cpus-per-task=1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLm_HREX_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --signal=B:USR1@1800

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

### Trap: stop GROMACS cleanly when walltime approaches ###
handle_usr1() {
    echo "=== USR1 received at $(date): requesting clean stop ==="
    touch STOP_MDRUN
}
trap handle_usr1 USR1

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

### Detect restart vs fresh start ###
shopt -s nullglob
rep_dirs=(rep*/)
shopt -u nullglob

n_rep=${#rep_dirs[@]}
if [[ "$n_rep" -eq 0 ]]; then
    echo "ERROR: No replica directories (rep*) found"
    exit 1
fi

n_cpt=0
for r in "${rep_dirs[@]}"; do
    [[ -f "$r/state.cpt" ]] && ((n_cpt++))
done

echo "Found $n_cpt / $n_rep replica checkpoints"

MD_ARGS=()

if [[ "$n_cpt" -eq 0 ]]; then
    echo "=== Fresh start ==="
elif [[ "$n_cpt" -eq "$n_rep" ]]; then
    echo "=== Restarting from checkpoints ==="
    MD_ARGS+=(-cpi state.cpt)
else
    echo "ERROR: Partial checkpoints detected ($n_cpt / $n_rep)" >&2
    exit 1
fi

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
    "${MD_ARGS[@]}"

echo "mdrun exited at $(date)"

### Copy results back ###
echo "============= Copying project outputs back to home ============="

echo "=== Copying non-replica files ==="
rsync -av --partial --inplace \
    --exclude 'rep*' \
    "$TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
    "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/"

echo "=== Copying replicas ==="
for rep in "$TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir/"/rep*; do
    [ -d "$rep" ] || continue
    repname=$(basename "$rep")
    echo ">>> Syncing $repname"
    rsync -av --partial --inplace \
        "$rep/" "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/$repname/"
done

echo "============= Copy complete at $(date) ============="

# Decide whether to requeue
FINAL_STEP=110000000
all_done=true

for r in rep*/md.log; do
    if ! grep -q "Writing checkpoint, step ${FINAL_STEP} " "$r"; then
        all_done=false
        break
    fi
done

if $all_done; then
    echo "=== Simulation fully completed (step $FINAL_STEP reached) ==="
    echo "=== No requeue ==="
else
    echo "=== Simulation stopped before final step ==="
    echo "=== Requeueing job ==="
    scontrol requeue "$SLURM_JOB_ID"
fi