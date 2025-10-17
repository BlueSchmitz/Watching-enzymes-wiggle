#!/bin/bash  

#SBATCH --job-name="hrex production run monomer 5EKY"   
#SBATCH --time=48:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu-a100
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=Research-AS-BN
#SBATCH --output=/scratch/blueschmitz/5EKY_monomeric_auto/hrex.out
#SBATCH --mail-type=ALL ##you can also set BEGIN/END

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail

# Load modules:  
module load 2025
module load openmpi/4.1.7
module load cuda
module load gromacs 

# Print modules 
module list

echo "Production HREX-MD with GROMACS and PLUMED"
srun gmx_mpi mdrun -multidir rep* -replex 2000 -f ../hrex.mdp -plumed ../plumed.dat -cpt 15 # 4 ps/0.002 ps = exchange every 2000 steps