#!/bin/bash

#SBATCH -J copy_over_project
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --cpus-per-task=8
#SBATCH --gpus=0
#SBATCH --output=copy_over_project_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

# Load modules:  
module load 2023
module list

project_dir=BtDERA
output_dir=6_HREX

### Copy project to scratch ###
echo "=== Copying project to scratch ==="
mkdir -p /gpfs/work1/0/prjs2080/$project_dir/outputs/$output_dir/
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir /gpfs/work1/0/prjs2080/$project_dir/outputs/$output_dir/
echo "=== Copy complete at $(date). ==="