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
