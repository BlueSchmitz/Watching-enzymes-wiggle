#!/bin/bash  

#SBATCH -J umbrella_1JCLm_analysis 
#SBATCH -t 00:30:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task 2
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=umbrella_analysis_%j.out

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
module load tqdm/4.66.1-GCCcore-12.3.0
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

# Output files
tpr_output="umbrella_tpr_paths.dat"
pullf_output="umbrella_pullf_paths.dat"

# Iterate over matching directories
for dir in ./COM_*; do
    # Check if the directory exists
    if [[ -d "$dir" ]]; then
        echo "Checking folder: $dir"

        # Find the umbrella tpr and pullf files (ignoring backups with '#' in their name)
        tpr_file=$(ls "$dir"/umbrella*.tpr 2>/dev/null | grep -v '#' | head -n 1)
        pullf_file=$(ls "$dir"/umbrella*_pullf.xvg 2>/dev/null | grep -v '#' | head -n 1)
        xtc_file=$(ls "$dir"/umbrella*.xtc 2>/dev/null | grep -v '#' | head -n 1)

        # Debug output
        if [[ -n "$tpr_file" ]]; then
            echo "Found TPR: $tpr_file"
        else
            echo "No umbrella*.tpr found in $dir"
        fi

        if [[ -n "$pullf_file" ]]; then
            echo "Found PULLF: $pullf_file"
        else
            echo "No umbrella*_pullf.xvg found in $dir"
        fi
        
        if [[ -n "$xtc_file" ]]; then
            echo "Found XTC: $xtc_file"
        else
            echo "No umbrella*.xtc found in $dir"
        fi

        # Run gmx energy for each energy term separately
        if [[ -n "$tpr_file" ]]; then
            edr_file="${tpr_file%.tpr}.edr"

            # COM-Pull-En
            energy_out_com="${tpr_file%.tpr}_com_pull_en.xvg"
            echo "Extracting COM-Pull-En energy..."
            echo COM-Pull-En | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -s "$tpr_file" -f "$edr_file" -o "$energy_out_com"

            # Potential
            energy_out_pot="${tpr_file%.tpr}_potential.xvg"
            echo "Extracting Potential energy..."
            echo Potential | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -s "$tpr_file" -f "$edr_file" -o "$energy_out_pot"

            # Total-Energy
            energy_out_tot="${tpr_file%.tpr}_total_energy.xvg"
            echo "Extracting Total-Energy..."
            echo Total-Energy | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -s "$tpr_file" -f "$edr_file" -o "$energy_out_tot"
        
            # Calculate distance 
            echo "Calculating distance between Lys167 NZ and Tyr259 OH..."
            apptainer exec $GROMACS_CONTAINER gmx_mpi distance -s "$tpr_file" -f "$xtc_file" -oall $dir/lys167_tyr259_dist.xvg -select "resid 167 and name NZ plus resid 259 and name OH"
        fi

        # Append to output files if both exist
        if [[ -n "$tpr_file" && -n "$pullf_file" ]]; then
            echo "$tpr_file" >> "$tpr_output"
            echo "$pullf_file" >> "$pullf_output"
        fi
    fi
done

echo "Paths saved to $tpr_output and $pullf_output"

apptainer exec $GROMACS_CONTAINER gmx_mpi wham \
        -it "$tpr_output" \
        -if "$pullf_output" \
        -o profile.xvg \
        -hist histo.xvg \
        -bsres profile_bootstrap.xvg \
        -nBootstrap 200 \
        -bs-method b-hist \
        -unit kJ \
        -xvg xmgrace \
        -temp 298

# Plot 
python $scripts/umbrella_analysis.py