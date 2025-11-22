#!/bin/bash  

#SBATCH --job-name="monomeric 5EKY set up"   
#SBATCH --time=05:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu-a100
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=Research-AS-BN
#SBATCH --output=/scratch/blueschmitz/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/5EKYm_setup%j.out
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
export GMXLIB=../../../../force_fields # make sure this is correct
minim=../../../mdp_templates/minim.mdp # Minimization mdp file
nvt=../../../mdp_templates/nvt.mdp # NVT mdp file
npt=../../../mdp_templates/npt.mdp # NPT mdp file
sanity=../../../mdp_templates/sanity_check.mdp # Sanity check mdp file
hrex=../../../mdp_templates/hrex.mdp # HREX mdp file
tempering=../../../mdp_templates/tempering.mdp # Tempering mdp file
pdb=../../inputs/5EKY_fill.BL00440001.pdb # Input PDB file (with correct protonation states)

# mkdir outputs
mkdir -p ./outputs/2_parametrization ./outputs/3_minimization ./outputs/4_equilibration ./outputs/5_sanity_checks ./outputs/6_HREX

### 2 Parametrize ###
echo "============= Parametrization with GROMACS =============" 
cd ./outputs/2_parametrization
# Generate topology and add hydrogens according to the chosen protonation states
  # Amber 99SB*-ILDN force field in combination with TIP3P water model
gmx_mpi pdb2gmx -f $pdb -o processed.gro -p topol.top -ff amber99sb-star-ildnp -water tip3p 
# Define the unit cell as described in paper: 15 A from protein to box edge = 1.5 nm
gmx_mpi editconf -f processed.gro -o boxed.gro -c -d 1.5 -bt cubic
  # -c: center the protein in the box
  # -d 1.5: minimum 15 Å distance from protein to box edge 
  # -bt cubic: cubic box 
# Solvate
gmx_mpi solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
# Add counterions (neutralize system)
gmx_mpi grompp -f $minim -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 # warning ignores net charge 
echo SOL | gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -neutral
# -pname NA -neutral: add Na⁺ to neutralize (paper)
cp solv_ions.gro ../3_minimization/solv_ions.gro
cp topol.top ../3_minimization/topol.top
cp posre.itp ../4_equilibration/posre.itp

### 3 Energy minimization ###
echo "============= Energy minimization with GROMACS ============="
cd ../3_minimization
gmx_mpi grompp -f $minim -c solv_ions.gro -p topol.top -o em.tpr
srun gmx_mpi mdrun -deffnm em
echo 10 0 | gmx_mpi energy -f em.edr -o potential.xvg # choose potential energy (10), 0 terminates input
python ../plot_xvg.py potential.xvg

### 4 Equilibration ###
echo "============= Equilibration with GROMACS ============="
mkdir -p ../4_equilibration
cp em.gro ../4_equilibration/em.gro
cp topol.top ../4_equilibration/topol.top
cd ../4_equilibration
# NVT Equilibration
gmx_mpi grompp -f $nvt -c em.gro -r em.gro -p topol.top -o nvt.tpr
srun gmx_mpi mdrun -deffnm nvt -cpt 15
echo 16 0 | gmx_mpi energy -f nvt.edr -o temperature.xvg # choose Temperature (16), 0 terminates input
python ../plot_xvg.py temperature.xvg
# NPT Equilibration
# Gradually reduce restraints from 1000 to 5 kJ mol−1 nm−2 by running 5 short NPT simulations of 500 ps each (5*500=2.5 ns)
for i in 1000 500 250 100 5;
do
  # Copy posre.itp 5 times and modify the force constant in each file
  cp posre.itp posre_$i.itp
  sed -i "s/\b1000\b/$i/g" posre_$i.itp # \b for whole word match
  echo "Modified posre.itp to posre_$i.itp."
  # Modify topol.top to include the correct posre file for each run
  cp topol.top topol_$i.top
  sed -i "s/posre.itp/posre_$i.itp/g" topol_$i.top
  echo "Modified topol.top to topol_$i.top."
  # Change absolute path to force field to relative path in each topol file
  sed -i "s|/home/blue/master/GROMACS/force_fields|../force_fields|g" topol_$i.top
  echo "Modified path in topol_$i.top."
done

# NPT Equilibration
for i in 1000 500 250 100 5;
do
  echo "Running NPT equilibration with restraints = ${i}"

  gmx_mpi grompp -f $npt \
             -c ${prev:-nvt.gro} \
             -r ${prev:-nvt.gro} \
             -p topol_${i}.top \
             -o npt_${i}.tpr \
              -maxwarn 1
  # ${prev:-nvt.gro} ensures the first run starts from NVT output, then continues from the last .gro
  srun gmx_mpi mdrun -deffnm npt_${i} -cpt 15
  echo 18 0 | gmx_mpi energy -f npt_${i}.edr -o pressure_${i}.xvg # choose Pressure (18), 0 terminates input
  echo 24 0 | gmx_mpi energy -f npt_${i}.edr -o density_${i}.xvg # choose Density (24), 0 terminates input
  prev=npt_${i}.gro
done
python ../plot_xvg.py pressure_*.xvg
python ../plot_xvg.py density_*.xvg
cp npt_5.gro ../6_HREX/npt_5.gro
cp topol_5.top ../6_HREX/topol_5.top

echo "Setup completed successfully."
