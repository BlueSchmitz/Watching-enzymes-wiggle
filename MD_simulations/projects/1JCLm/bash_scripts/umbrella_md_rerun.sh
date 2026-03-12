#!/bin/bash  

#SBATCH -J umbrella_1JCLm 
#SBATCH -t 10:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --array=1-11
#SBATCH --ntasks=8
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=umbrella_MD_%A_%a.out

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
# Array task ID
TASK_ID=$SLURM_ARRAY_TASK_ID

# COM windows that still need to run
target_windows=(
COM_0.564
COM_0.661
COM_1.038
COM_1.084
COM_1.120
COM_1.151
COM_1.689
COM_1.745
COM_1.795
COM_2.570
COM_2.683
)

scripts=()

for win in "${target_windows[@]}"; do
    for s in "$win"/frame-*"_run-umbrella.sh"; do
        [[ -f "$s" ]] && scripts+=("$s")
    done
done

# Ensure deterministic order
IFS=$'\n' scripts=($(printf "%s\n" "${scripts[@]}" | sort))
unset IFS

echo "Scripts that will run:"
printf '%s\n' "${scripts[@]}"

# check if tasks and scripts match up
num_scripts=${#scripts[@]}

if (( TASK_ID > num_scripts )); then
    echo "Error: TASK_ID $TASK_ID exceeds number of scripts $num_scripts"
    exit 1
fi

# Select script for this task (TASK_ID is 1-based, array is 0-based)
script=${scripts[$((TASK_ID-1))]}

if [[ -z "$script" ]]; then
    echo "No script found for task $TASK_ID"
    exit 1
fi

# Extract directory and frame number
script_dir=$(dirname "$script")
frame_number=$(echo "$script" | grep -oP 'frame-\K[0-9]+')

# 1 thread per rank
export OMP_NUM_THREADS=1

# Run script 
echo "Current directory: $(pwd)"
echo "Task $TASK_ID: Running $script..."

export OMP_NUM_THREADS=1

# Run inside subshell to ensure we return to the original directory after execution
(
  cd "$script_dir" || exit 1
  echo "Working directory: $(pwd)"
  bash "$(basename "$script")"
)

echo "Task $TASK_ID: Finished $script."