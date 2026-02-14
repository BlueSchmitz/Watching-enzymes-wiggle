#!/bin/bash  

#SBATCH -J umbrella_1JCLm 
#SBATCH -t 02:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task 2
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./umbrella_%j.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_TASKS=$SLURM_NTASKS

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=1JCLm
output_dir=8_umbrella
pH=7

export GMXLIB=$TMPDIR/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts
pdb=$TMPDIR/MD_simulations/projects/$project_dir/inputs/*.pdb

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module list

### Copy project to scratch ###
echo "=== Copying project to scratch ==="
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/ "$TMPDIR"
cd $TMPDIR/MD_simulations/projects/$project_dir

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

# 1 Steered MD to pull tail into sctive site 
# remove restraints from npt_5.gro 
#mkdir -p ./outputs/$output_dir
#cp ./outputs/4_equilibration/npt_5.gro ./outputs/$output_dir/npt.gro
#cp ./outputs/4_equilibration/topol_5.top ./outputs/$output_dir/topol_5.top
cd ./outputs/$output_dir/
#sed '/#ifdef POSRES/,/#endif/ s|#include "posre_.*\.itp"|#include "posre_core_CA.itp"|' topol_5.top > topol.top # Change the protein restraint block (POSRES)
#echo "Changed position restraints from topol.top to restrain core_CA."
# define groups for pulling: tail_COM and active_site_COM
#apptainer exec $GROMACS_CONTAINER gmx_mpi make_ndx -f npt.gro -o index.ndx << EOF
#r 257-259
#name 18 tail_COM 
#r 166-168
#name 19 active_site_COM
#r 1-247 & a CA
#name 20 core_CA
#q
#EOF
# restrain the core CA atoms during pulling
#echo 20 | apptainer exec $GROMACS_CONTAINER gmx_mpi genrestr -f npt.gro -n index.ndx -o posre_core_CA.itp -fc 10 10 10
# run steered MD
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/pull.mdp -c ./pull1/conf128.gro -n index.ndx -p topol.top -o pull.tpr -r ./pull1/conf128.gro
apptainer exec $GROMACS_CONTAINER mpirun -np 8 gmx_mpi mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg -v -ntomp 2 -maxh 1.9 \

# 2 Assemble COM distances indexed by frame number
echo 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep
# compute distances
#nframes=$(ls conf*.gro | wc -l)
#for (( i=0; i<${nframes}; i++ ))
#do
#    apptainer exec $GROMACS_CONTAINER gmx_mpi distance -s pull.tpr \
#        -f conf${i}.gro \
#        -n index.ndx \
#        -select 'com of group "tail_COM" plus com of group "active_site_COM"' \
#        -oall dist${i}.xvg 
#done
# compile summary
#touch summary_distances.dat
#for (( i=0; i<${nframes}; i++ ))
#do
#    d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
#    echo "${i} ${d}" >> summary_distances.dat
#    rm dist${i}.xvg
#done

# 3 Prepare umbrella sampling windows
#python $scripts/setupUmbrella.py summary_distances.dat 0.1 ../../bash_scripts/umbrella_template.sh &> caught-output.txt #edit!!!

#################################################
# 4 Run umbrella sampling windows
# Collect all .sh scripts in COM_* directories
#scripts=()
#failed_only=false
#for script in COM_*/frame-*_run-umbrella.sh; do
#    if [ -f "$script" ]; then
#        script_dir=$(dirname "$script")
#        fail_file="$script_dir/FAIL"
#        
#        # If --failed-only is set, check for FAIL file
#        if $failed_only && [ ! -f "$fail_file" ]; then
#            continue
#        fi
#        
#        scripts+=("$script")
#    fi
#done
# If no scripts found, exit
#if [ ${#scripts[@]} -eq 0 ]; then
#    echo "No shell scripts found to execute."
#    exit 1
#fi
#echo "Found ${#scripts[@]} umbrella windows"

#for script in "${scripts[@]}"; do
#    script_dir=$(dirname "$script")  # Extract the directory path
        
    # Extract the frame number (XXX) from the filename
#    frame_number=$(echo "$script" | grep -oP 'frame-\K[0-9]+')
    
    # Define the log file name as frame-XXX.log
#    log_file="$script_dir/frame-${frame_number}.log"

#    echo "Running $script... (Logging to $log_file)"
        
    # Execute the script and log both stdout and stderr
#    bash "$script" > "$log_file" 2>&1

    # Check for execution success or failure
#    if [ $? -ne 0 ]; then
#        echo "Execution of $script failed. Check $log_file for details."
#    else
#        echo "Execution of $script completed successfully."
#            
#        # Remove the FAIL file if execution was successful
#        fail_file="$script_dir/FAIL"
#        if [ -f "$fail_file" ]; then
#            rm "$fail_file"
#            echo "Removed FAIL file from $script_dir."
#        fi
#    fi
#done
#echo "All scripts executed."