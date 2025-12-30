#!/bin/bash  

#SBATCH -J umbrella 
#SBATCH -t 10:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./outputs/umbrella_%j.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_TASKS=$SLURM_NTASKS

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

# Trap errors and print line number + command
trap 'echo "Error in ${BASH_SOURCE[0]} at line ${LINENO}: ${BASH_COMMAND}"' ERR

# path to gromacs apptainer container
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif

# Load modules:  
module load 2023
#module load OpenMPI/4.1.5-GCC-12.3.0
module load matplotlib/3.7.2-gfbf-2023a

# Copy input files to scratch
#mpicopy $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/
#echo "Contents of $TMPDIR:"
#tree -a -L 10 $TMPDIR
#cd $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs

# Set paths for mdp_templates and pyscripts
mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/inputs # Input PDB file (with correct protonation states)
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts

# 1 Steered MD to pull tail into sctive site 
# remove restraints from npt_5.gro 
mkdir -p ./outputs/8_umbrella
cp ./outputs/4_equilibration/npt_5.gro ./outputs/8_umbrella/npt.gro
cp ./outputs/4_equilibration/topol_5.top ./outputs/8_umbrella/topol_5.top
sed '/#ifdef POSRES/,/#endif/ s|#include "posre_.*\.itp"|#include "posre_core_CA.itp"|' topol_5.top > topol.top # Change the protein restraint block (POSRES)
echo "Changed position restraints from topol.top to restrain core_CA."
cd ./outputs/8_umbrella/
# define groups for pulling: tail_COM and active_site_COM
gmx_mpi make_ndx -f npt.gro -o index.ndx << EOF
r 257-259
name 19 tail_COM 
r 166-168
name 20 active_site_COM
r 1-247 & a CA
name 21 core_CA
q
EOF
# restrain the core CA atoms during pulling
echo 21 | gmx genrestr -f npt.gro -n index.ndx -o posre_core_CA.itp -fc 10 10 10
# run steered MD
gmx_mpi grompp -f $mdp/pull.mdp -c npt.gro -n index.ndx -p topol.top -o pull.tpr
gmx_mpi mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg -v

# 2 Assemble COM distances indexed by frame number
echo 0 | gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep
# compute distances
nframes=$(ls conf*.gro | wc -l)
for (( i=0; i<${nframes}; i++ ))
do
    gmx_mpi distance -s pull.tpr \
        -f conf${i}.gro \
        -n index.ndx \
        -select 'distance between com of group "tail_COM" and com of group "active_site_COM"' \
        -oall dist${i}.xvg 
done
# compile summary
touch summary_distances.dat
for (( i=0; i<${nframes}; i++ ))
do
    d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
    echo "${i} ${d}" >> summary_distances.dat
    rm dist${i}.xvg
done

# 3 Prepare umbrella sampling windows
python setupUmbrella.py summary_distances.dat 0.1 umbrella_template.sh &> caught-output.txt #edit!!!

# 4 Run umbrella sampling windows
# Collect all .sh scripts in COM_* directories
scripts=()
for script in COM_*/frame-*_run-umbrella.sh; do
    if [ -f "$script" ]; then
        script_dir=$(dirname "$script")
        fail_file="$script_dir/FAIL"
        
        # If --failed-only is set, check for FAIL file
        if $failed_only && [ ! -f "$fail_file" ]; then
            continue
        fi
        
        scripts+=("$script")
    fi
done
# If no scripts found, exit
if [ ${#scripts[@]} -eq 0 ]; then
    echo "No shell scripts found to execute."
    exit 1
fi
echo "Found ${#scripts[@]} umbrella windows"

for script in "${scripts[@]}"; do
    script_dir=$(dirname "$script")  # Extract the directory path
        
    # Extract the frame number (XXX) from the filename
    frame_number=$(echo "$script" | grep -oP 'frame-\K[0-9]+')
    
    # Define the log file name as frame-XXX.log
    log_file="$script_dir/frame-${frame_number}.log"

    echo "Running $script... (Logging to $log_file)"
        
    # Execute the script and log both stdout and stderr
    bash "$script" > "$log_file" 2>&1

    # Check for execution success or failure
    if [ $? -ne 0 ]; then
        echo "Execution of $script failed. Check $log_file for details."
    else
        echo "Execution of $script completed successfully."
            
        # Remove the FAIL file if execution was successful
        fail_file="$script_dir/FAIL"
        if [ -f "$fail_file" ]; then
            rm "$fail_file"
            echo "Removed FAIL file from $script_dir."
        fi
    fi
done
echo "All scripts executed."
