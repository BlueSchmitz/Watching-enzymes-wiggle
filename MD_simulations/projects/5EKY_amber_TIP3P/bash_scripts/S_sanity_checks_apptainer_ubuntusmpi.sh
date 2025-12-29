#!/bin/bash  

#SBATCH -J sanity_checks_Ec5EKYm_apptainer  
#SBATCH -t 00:05:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./outputs/sanity_checks_Ec5EKYm_apptainer_%j.out

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
cd $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
export GMXLIB=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/force_fields # make sure this is correct
project=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs
mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/inputs # Input PDB file (with correct protonation states)
pyscripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts

### Sanity checks ###
echo "Sanity checks"
cd 5_sanity_checks
# 1. Sanity check of created topologies: compare energies between original and 1.00 scaled system
# produce tpr + short trajectory (use processed.top)
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p processed.top -o make_traj.tpr -maxwarn 1
apptainer exec $GROMACS_CONTAINER gmx_mpi mdrun -deffnm make_traj
cp make_traj.xtc traj.xtc
# produce tpr for original topology
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p processed.top -o topol_pro.tpr -maxwarn 1
# compute energies using original
apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -rerun traj.xtc -s topol_pro.tpr -e ener_pro.edr -g rerun_pro.log -ntomp $OMP_NUM_THREADS
# produce tpr from scaled topology (use same mdp / coords)
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p scaled_1.00.top -o scaled_1.00.tpr -maxwarn 1
# Recompute energies using the scaled topology and the previously made trajectory
apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -rerun traj.xtc -s scaled_1.00.tpr -e ener_scaled_1.00.edr -g rerun_scaled_1.00.log -ntomp $OMP_NUM_THREADS
# compare energies
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_pro.edr -o energies_pro.xvg -xvg none
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_scaled_1.00.edr -o energies_scaled_1.00.xvg -xvg none
python $pyscripts/energy_comparison.py energies_pro.xvg energies_scaled_1.00.xvg > energy_diff_1.00.log
# 2. Sanity check of created topologies: compare energies between 1.00 scaled system and 0.5 scaled system
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c ./npt_5.gro -p ./scaled_0.5_all.top -o scaled_0.5_all.tpr -maxwarn 2
apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -rerun traj.xtc -s scaled_0.5_all.tpr -e ener_scaled_0.5_all.edr -g rerun_scaled_0.5_all.log -ntomp $OMP_NUM_THREADS
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_scaled_0.5_all.edr -o energies_scaled_0.5_all.xvg -xvg none
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
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f sanity_check_rep0.mdp -c ../npt_5.gro -p ../processed.top -o topol.tpr -maxwarn 2
# Prepare tpr for rep1 (scaled topology)
cd ../rep1
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f ./sanity_check_rep1.mdp -c ../npt_5.gro -p ../scaled_1.00.top -o topol.tpr -maxwarn 2
cd ..
# Run a short HREX simulation with 2 replicas
apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -multidir rep0 rep1 -replex 200 -hrex -plumed plumed.dat -dlb no -ntomp $OMP_NUM_THREADS -v -nb cpu
# Inspect logs
echo "Replica-exchange acceptance"
grep "Repl" rep0/md.log

# save enerything back to home
#echo "=== Inspecting $TMPDIR/outputs_tmp before copying back ==="
#ls -lR $TMPDIR/outputs_tmp
#echo "=== End of inspection ==="

#echo "=== Copying project outputs back to home ==="
#rsync -av --exclude 'rleveson.*' \
#      $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/ \
#      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/
#echo "=== Copy complete ==="