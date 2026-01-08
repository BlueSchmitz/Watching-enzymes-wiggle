#!/bin/bash  

#SBATCH -J sanity_checks_Ec5EKYm_apptainer  
#SBATCH -t 01:00:00
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

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts

### Prepare HREX ###
echo "Sanity checks"
cp ./outputs/4_equilibration/npt_5.gro ./outputs/6_HREX/npt_5.gro
cp ./outputs/4_equilibration/topol_5.top ./outputs/6_HREX/topol_5.top
cd ./outputs/6_HREX
# Files we need: npt_5.gro as the starting structure, topol.top without pointers to posre.itp (no more restraints)
# 1. Remove constraints from topol.top by commenting out the line that includes posre_5.itp
cp topol_5.top topol_prod.top
sed -i '/#ifdef POSRES/,/#endif/d' topol_prod.top # Remove the protein restraint block (POSRES)
echo "Removed position restraints from topol_prod.top."
# 2. Re-name topology file (already self-cotained)
#gmx_mpi grompp -f tempering.mdp -c npt_5.gro -p topol_prod.top -pp processed.top -o dummy.tpr -r npt_5.gro -maxwarn 2
mv topol_prod.top processed.top
echo "processed.top is self-contained."
# 3. Edit the processed.top file to indicate which atoms we want to scale (marked with an _ after the residue name)
python $scripts/scale_residues.py > processed_scaled.top
echo "Generated processed_scaled.top with selected residues for scaling."
# 4. Scale the Hamiltonian of the selected atoms by the factors 1.00, 0.95, 0.91, 0.87, 0.83, 0.79, 0.76, 0.72, 0.69, 0.66, 0.63, and 0.60
for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
do 
  mkdir -p ./rep${i} # create directory for each replica
  plumed partial_tempering ${i} < processed_scaled.top  > ./rep${i}/scaled_${i}.top
  echo "Generated scaled_${i}.top in ./rep${i} with scaling factor ${i}."
  cd ./rep${i}
  gmx_mpi grompp -f ../$mdp/hrex.mdp -c ../npt_5.gro -p scaled_${i}.top -o topol.tpr -maxwarn 1
  echo "Generated topol.tpr from scaled_${i}.top in ./rep${i}."
  cd ..
done 
echo "Generated all scaled topology files and tpr files for HREX replicas."
# 5. Empty plumed.dat file (no additional CVs, only HREX)
: > plumed.dat # empty plumed file

### Sanity checks ###
cp ./npt_5.gro ../5_sanity_checks/npt_5.gro
cp ./processed.top ../5_sanity_checks/processed.top
cp ./rep1.00/scaled_1.00.top ../5_sanity_checks/scaled_1.00.top
cd ../5_sanity_checks

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
python $scripts/energy_comparison.py energies_pro.xvg energies_scaled_1.00.xvg > energy_diff_1.00.log

# 2. Sanity check of created topologies: compare energies between 1.00 scaled system and 0.5 scaled system
# scale all residues to 0.5
python $scripts/scale_residues.py > processed_scaled_all.top
plumed partial_tempering 0.5 < processed_scaled_all.top  > scaled_0.5_all.top
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c ./npt_5.gro -p ./scaled_0.5_all.top -o scaled_0.5_all.tpr -maxwarn 2
apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -rerun traj.xtc -s scaled_0.5_all.tpr -e ener_scaled_0.5_all.edr -g rerun_scaled_0.5_all.log -ntomp $OMP_NUM_THREADS
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_scaled_0.5_all.edr -o energies_scaled_0.5_all.xvg -xvg none
python $scripts/energy_comparison.py energies_scaled_1.00.xvg energies_scaled_0.5_all.xvg > energy_diff_0.5_all.log

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