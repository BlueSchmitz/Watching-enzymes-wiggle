#!/bin/bash  

#SBATCH --job-name="monomeric 5EKY set up"   
#SBATCH --time=48:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu-a100
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=Research-AS-BN
#SBATCH --output=/scratch/blueschmitz/5EKY_monomeric_auto/setup.out
#SBATCH --mail-type=ALL ##you can also set BEGIN/END

: '
Folder structure:
.
└── MD_simulations/
    ├── scripts /
    │   ├── model_missing_residues.py
    │   ├── get_ali_file.py
    │   ├── compare_protonation.py
    │   ├── plot_xvg.py
    │   ├── scale_all_residues.py
    │   ├── scale_residues.py
    │   └── energy_comparison.py
    ├── mdp_templates/
    │   ├── minim.mdp
    │   ├── nvt.mdp
    │   ├── npt.mdp
    │   ├── sanity_check.mdp
    │   └── hrex.mdp
    ├── force_fields/
    │   ├── amber99sb-star-ildnp.ff/
    └── projects/
        └── 5EKY_monomeric/
            ├── bash_scripts
            ├── inputs/
            │   └── pdb
            └── outputs/
                ├── 1_protonation
                ├── 2_parametrization
                ├── 3_minimization
                ├── 4_equilibration
                ├── 5_sanity_checks/
                │   ├── rep0
                │   └── rep1
                └── 6_HREX/
                    ├── rep0.60
                    ├── rep0.63
                    └── ...

Run this script from the projects/5EKY_monomeric/ directory, it contains relative paths
'
# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail

# Load modules:  
module load 2025
module load openmpi/4.1.7
module load cuda
module load gromacs 
module load python/3.11.9
module load py-matplotlib/3.9.2
module load py-numpy/1.26.4

# Print modules 
module list

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
export GMXLIB=../../../force_fields # make sure this is correct
minim=../../mdp_templates/minim.mdp # Minimization mdp file
nvt=../../mdp_templates/nvt.mdp # NVT mdp file
npt=../../mdp_templates/npt.mdp # NPT mdp file
sanity=../../mdp_templates/sanity_check.mdp # Sanity check mdp file
hrex=../../mdp_templates/hrex.mdp # HREX mdp file
tempering=../../mdp_templates/tempering.mdp # Tempering mdp file
pdb=../../inputs/5EKY_fill.BL00440001.pdb # Input PDB file (with correct protonation states)

# mkdir outputs
mkdir -p ./outputs/2_parametrization ./outputs/3_minimization ./outputs/4_equilibration ./outputs/5_sanity_checks ./outputs/6_HREX

### 5 Set up HREX-MD with PLUMED ###
echo "Setting up HREX-MD with GROMACS and PLUMED"
cd ../6_HREX
# Files we need: npt_5.gro as the starting structure, topol.top without pointers to posre.itp (no more restraints)
# 1. Remove constraints from topol.top by commenting out the line that includes posre_5.itp
cp topol_5.top topol_prod.top
sed -i '/#ifdef POSRES/,/#endif/d' topol_prod.top # Remove the protein restraint block (POSRES)
sed -i '/#ifdef POSRES_WATER/,/#endif/d' topol_prod.top # Remove the water restraint block (POSRES_WATER)
echo "Removed position restraints from topol_prod.top."
# 2. Generate a self-contained topology file
gmx_mpi grompp -f $tempering -c npt_5.gro -p topol_prod.top -pp processed.top -o dummy.tpr -r npt_5.gro -maxwarn 2
echo "Generated self-contained processed.top without position restraints."
# 3. Edit the processed.top file to indicate which atoms we want to scale (marked with an _ after the residue name)
python scale_residues.py > processed_scaled.top
echo "Generated processed_scaled.top with selected residues for scaling."
# 4. Scale the Hamiltonian of the selected atoms by the factors 1.00, 0.95, 0.91, 0.87, 0.83, 0.79, 0.76, 0.72, 0.69, 0.66, 0.63, and 0.60
for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
do 
  mkdir -p ./rep${i} # create directory for each replica
  plumed partial_tempering ${i} < processed_scaled.top  > ./rep${i}/scaled_${i}.top
  echo "Generated scaled_${i}.top in ./rep${i} with scaling factor ${i}."
  cd ./rep${i}
  gmx_mpi grompp -f $hrex -c ../npt_5.gro -p scaled_${i}.top -o topol.tpr
  echo "Generated topol.tpr from scaled_${i}.top in ./rep${i}."
  cd ..
done 
echo "Generated all scaled topology files and tpr files for HREX replicas."
# 5. Empty plumed.dat file (no additional CVs, only HREX)
: > plumed.dat # empty plumed file

### 6 Sanity checks ### 
echo "Sanity checks"
cp npt_5.gro ../5_sanity_checks/npt_5.gro
cp processed.top ../5_sanity_checks/processed.top
cp ./rep1.00/scaled_1.00.top ../5_sanity_checks/scaled_1.00.top
cd ../5_sanity_checks
# 1. Sanity check of created topologies: compare energies between original and 1.00 scaled system
# produce tpr + short trajectory (use processed.top)
gmx_mpi grompp -f $sanity_check -c npt_5.gro -p processed.top -o make_traj.tpr -maxwarn 1
gmx_mpi mdrun -deffnm make_traj
cp make_traj.xtc traj.xtc
# produce tpr for original topology
gmx_mpi grompp -f $sanity_check -c npt_5.gro -p processed.top -o topol_pro.tpr -maxwarn 1
# compute energies using original
gmx_mpi mdrun -rerun traj.xtc -s topol_pro.tpr -e ener_pro.edr -g rerun_pro.log
# produce tpr from scaled topology (use same mdp / coords)
gmx_mpi grompp -f $sanity_check -c npt_5.gro -p scaled_1.00.top -o scaled_1.00.tpr -maxwarn 1
# Recompute energies using the scaled topology and the previously made trajectory
gmx_mpi mdrun -rerun traj.xtc -s scaled_1.00.tpr -e ener_scaled_1.00.edr -g rerun_scaled_1.00.log
# compare energies
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | gmx_mpi energy -f ener_pro.edr -o energies_pro.xvg -xvg none
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | gmx_mpi energy -f ener_scaled_1.00.edr -o energies_scaled_1.00.xvg -xvg none
python energy_comparison.py energies_pro.xvg energies_scaled_1.00.xvg > energy_diff_1.00.log
# 2. Sanity check of replica-exchange implementation
# Run a short HREX with two equivalent topology files (topol.tpr and scaled_1.00.tpr)
mkdir -p ./rep0 ./rep1
: > plumed.dat # empty plumed file
cp plumed.dat ./rep0/plumed.dat
cp plumed.dat ./rep1/plumed.dat
# copy mdp files and set different random seeds for both replicas
cp $sanity_check ./rep0/sanity_check_rep0.mdp
cp $sanity_check ./rep1/sanity_check_rep1.mdp
sed -i 's/^gen_seed.*/gen_seed = 12345/' ./rep0/sanity_check_rep0.mdp
sed -i 's/^gen_seed.*/gen_seed = 67890/' ./rep1/sanity_check_rep1.mdp
# Prepare tpr for rep0 (original topology)
cd rep0
gmx_mpi grompp -f sanity_check_rep0.mdp -c ../npt_5.gro -p ../processed.top -o topol.tpr -maxwarn 1
# Prepare tpr for rep1 (scaled topology)
cd ../rep1
gmx_mpi grompp -f ./sanity_check_rep1.mdp -c ../npt_5.gro -p ../scaled_1.00.top -o topol.tpr -maxwarn 1
cd ..
# Run a short HREX simulation with 2 replicas
srun gmx_mpi mdrun -multidir rep0 rep1 -replex 200 -hrex -plumed plumed.dat
# Inspect logs
echo "Replica-exchange acceptance"
grep "Repl" rep0/md.log