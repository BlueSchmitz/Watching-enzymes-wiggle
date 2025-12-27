#!/bin/bash  

#SBATCH -J sanity_checks_Ec5EKYm  
#SBATCH -t 00:20:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./outputs/sanity_checks_Ec5EKYm%j.out

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

function copy_back_results {
    echo "=== Copying results back to home ==="
    rsync -av \
          --exclude 'rleveson.*' \
          --exclude 'sanity_checks_Ec5EKYm*.out' \
        $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/ \
        $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/
    echo "=== Copy complete ==="
}

# Trap errors and print line number + command
trap 'echo "Error in ${BASH_SOURCE[0]} at line ${LINENO}: ${BASH_COMMAND}"' ERR
trap copy_back_results EXIT

# Load modules:  
module load 2024
module load mpicopy/4.2-gompi-2024a
module load Miniconda3/24.7.1-0

# Copy input files to scratch
mpicopy $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/
#echo "Contents of $TMPDIR:"
#tree -a -L 10 $TMPDIR
cd $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs
#cd $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs

# activate environment with GROMACS installed
unset CONDA_SHLVL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gromacs2020
# Force use of Conda MPI
export PATH=$CONDA_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v "$CONDA_PREFIX/lib" | paste -sd:)

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
export GMXLIB=$TMPDIR/MD_simulations/force_fields # make sure this is correct
project=$TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs
mdp=$TMPDIR/MD_simulations/mdp_templates
pdb=$TMPDIR/MD_simulations/projects/5EKY_monomeric/inputs # Input PDB file (with correct protonation states)
pyscripts=$TMPDIR/MD_simulations/scripts

### Sanity checks ###
echo "Sanity checks"
cd 5_sanity_checks
# 1. Sanity check of created topologies: compare energies between original and 1.00 scaled system
# produce tpr + short trajectory (use processed.top)
gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p processed.top -o make_traj.tpr -maxwarn 1
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm make_traj
cp make_traj.xtc traj.xtc
# produce tpr for original topology
gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p processed.top -o topol_pro.tpr -maxwarn 1
# compute energies using original
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -rerun traj.xtc -s topol_pro.tpr -e ener_pro.edr -g rerun_pro.log
# produce tpr from scaled topology (use same mdp / coords)
gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p scaled_1.00.top -o scaled_1.00.tpr -maxwarn 1
# Recompute energies using the scaled topology and the previously made trajectory
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -rerun traj.xtc -s scaled_1.00.tpr -e ener_scaled_1.00.edr -g rerun_scaled_1.00.log
# compare energies
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | gmx_mpi energy -f ener_pro.edr -o energies_pro.xvg -xvg none
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | gmx_mpi energy -f ener_scaled_1.00.edr -o energies_scaled_1.00.xvg -xvg none
python $pyscripts/energy_comparison.py energies_pro.xvg energies_scaled_1.00.xvg > energy_diff_1.00.log
# 2. Sanity check of created topologies: compare energies between 1.00 scaled system and 0.5 scaled system
gmx_mpi grompp -f $mdp/sanity_check.mdp -c ./npt_5.gro -p ./scaled_0.5_all.top -o scaled_0.5_all.tpr -maxwarn 2
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -rerun traj.xtc -s scaled_0.5_all.tpr -e ener_scaled_0.5_all.edr -g rerun_scaled_0.5_all.log
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | gmx_mpi energy -f ener_scaled_0.5_all.edr -o energies_scaled_0.5_all.xvg -xvg none
python $pyscripts/energy_comparison.py energies_scaled_1.00.xvg energies_scaled_0.5_all.xvg > energy_diff_0.5_all.log
# 3. Sanity check of replica-exchange implementation
# Run a short HREX with two equivalent topology files (topol.tpr and scaled_1.00.tpr)
mkdir -p ./rep0 ./rep1
: > plumed.dat # empty plumed file
cp plumed.dat ./rep0/plumed.dat
cp plumed.dat ./rep1/plumed.dat
# copy mdp files and set different random seeds for both replicas
cp $mdp/sanity_check.mdp ./rep0/sanity_check_rep0.mdp
cp $mdp/sanity_check.mdp ./rep1/sanity_check_rep1.mdp
sed -i 's/^gen_seed.*/gen_seed = 12345/' ./rep0/sanity_check_rep0.mdp
sed -i 's/^gen_seed.*/gen_seed = 67890/' ./rep1/sanity_check_rep1.mdp
# Prepare tpr for rep0 (original topology)
cd rep0
gmx_mpi grompp -f sanity_check_rep0.mdp -c ../npt_5.gro -p ../processed.top -o topol.tpr -maxwarn 2
# Prepare tpr for rep1 (scaled topology)
cd ../rep1
gmx_mpi grompp -f ./sanity_check_rep1.mdp -c ../npt_5.gro -p ../scaled_1.00.top -o topol.tpr -maxwarn 2
cd ..
# Run a short HREX simulation with 2 replicas
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -multidir rep0 rep1 -replex 200 -hrex -plumed plumed.dat -ntomp $OMP_NUM_THREADS -v -nb cpu -dlb no
# Inspect logs
echo "Replica-exchange acceptance"
grep "Repl" rep0/md.log

# save enerything back to home
#echo "=== Inspecting $TMPDIR/outputs_tmp before copying back ==="
#ls -lR $TMPDIR/outputs_tmp
#echo "=== End of inspection ==="

#echo "=== Copying project outputs back to home ==="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude 'sanity_checks_Ec5EKYm*.out' \
      $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/
#echo "=== Copy complete ==="