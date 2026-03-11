#!/bin/bash  

#SBATCH -J BtDERA_setup  
#SBATCH -t 02:30:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16 
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./BtDERA_setup_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_TASKS=$SLURM_NTASKS

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

Run this script from the projects/project_dir/ directory, it contains relative paths.
This script is written for running inside an apptainer on the Snellius super computer.
'

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=BtDERA
pH=7

# Set paths for mdp_templates, force_fields and pdb file (to change paths quickly)
export GMXLIB=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/inputs/*.pdb

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module list

# Create temporary directory on scratch for this job
tmpdir=$(mktemp -d /gpfs/scratch1/shared/rleveson/Blue/tmp.XXXXXX)
mkdir -p "$tmpdir/MD_simulations/projects/"

# Copy project to scratch 
echo "=== Copying project to scratch ==="
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir "$tmpdir/MD_simulations/projects/"
cd $tmpdir/MD_simulations/projects/$project_dir

# Function to copy back results when error occurs and before the script exits
function copy_back_results {
    set +e +u # Disable exit on error for this function
    echo "=== Copying results back to home at $(date). ==="
    if [[ -d "$tmpdir/MD_simulations/projects/$project_dir/outputs/" ]]; then
        rsync -av --partial --inplace \
          "$tmpdir/MD_simulations/projects/$project_dir/outputs/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/"
        echo "=== Copy complete at $(date) ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
    # Clean up temporary directory
    rm -rf "$tmpdir"
}
trap copy_back_results EXIT

# mkdir outputs directories 
mkdir -p ./outputs/1_protonation ./outputs/2_parametrization ./outputs/3_minimization ./outputs/4_equilibration ./outputs/5_sanity_checks ./outputs/6_HREX

### 1 Protonation ###
echo "============= Protonation states assignment with PDB2PQR and PROPKA 3 ============="
cd ./outputs/1_protonation
# Assign protonation states at the desired pH (and pH7) on the basis of the PROPKA 3 estimate
apptainer exec $PDB2PQR_CONTAINER pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=$pH $pdb ./${project_dir}_pH${pH}.pqr
apptainer exec $PDB2PQR_CONTAINER pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=7 $pdb ./${project_dir}_pH7.pqr
# Convert .pqr to .pdb
cp ${project_dir}_pH${pH}.pqr ${project_dir}_pH${pH}.pdb
cp ${project_dir}_pH7.pqr ${project_dir}_pH7.pdb
# change the names of the residues with differing H numbers at pH 6.8 in comparison to pH 7 so that GROMACS can rebuild them
python $scripts/compare_protonation.py ${project_dir}_pH${pH}.pdb ${project_dir}_pH7.pdb differences.txt $pdb ${project_dir}_pro.pdb
echo "Differences in protonation states (pH ${pH} vs pH 7) written to differences.txt. Modified PDB written to ${project_dir}_pro.pdb."
cp ./${project_dir}_pro.pdb ../2_parametrization/${project_dir}_pro.pdb

### 2 Parametrize ###
echo "============= Parametrization with GROMACS =============" 
cd ../2_parametrization
# Generate topology and initial structure files
  # Amber 99SB*-ILDN force field in combination with TIP3P water model
apptainer exec $GROMACS_CONTAINER gmx_mpi pdb2gmx -f ${project_dir}_pro.pdb -o processed.gro -p topol.top -ff amber99sb-star-ildnp -water tip3p 
# Define the unit cell: 15 A from protein to box edge = 1.5 nm
apptainer exec $GROMACS_CONTAINER gmx_mpi editconf -f processed.gro -o boxed.gro -c -d 1.5 -bt cubic
  # -c: center the protein in the box
  # -d 1.5: minimum 15 Å distance from protein to box edge 
  # -bt cubic: cubic box 
# Solvate
apptainer exec $GROMACS_CONTAINER gmx_mpi solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
# Add counterions (neutralize system)
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/minim.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 # warning ignores net charge 
echo SOL | apptainer exec $GROMACS_CONTAINER gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -neutral
  # -pname NA -neutral: add Na⁺ to neutralize (paper)
cp solv_ions.gro ../3_minimization/solv_ions.gro
cp topol.top ../3_minimization/topol.top
cp posre.itp ../4_equilibration/posre.itp

### 3 Energy minimization ###
echo "============= Energy minimization with GROMACS ============="
cd ../3_minimization
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm em
echo 10 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f em.edr -o potential.xvg # choose potential energy (10), 0 terminates input
python $scripts/plot_xvg.py potential.xvg

### 4 Equilibration ###
echo "============= Equilibration with GROMACS ============="
mkdir -p ../4_equilibration
cp em.gro ../4_equilibration/em.gro
cp topol.top ../4_equilibration/topol.top
cd ../4_equilibration

# NVT Equilibration
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm nvt -cpt 15
echo 16 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f nvt.edr -o temperature.xvg # choose Temperature (16), 0 terminates input
python $scripts/plot_xvg.py temperature.xvg
# Preparation for NPT Equilibration
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
done

# NPT Equilibration
for i in 1000 500 250 100 5;
do
  echo "Running NPT equilibration with restraints = ${i}"

  apptainer exec $GROMACS_CONTAINER gmx_mpi grompp \
             -f $mdp/npt.mdp \
             -c ${prev:-nvt.gro} \
             -r ${prev:-nvt.gro} \
             -p topol_${i}.top \
             -o npt_${i}.tpr \
              -maxwarn 1
  # ${prev:-nvt.gro} ensures the first run starts from NVT output, then continues from the last .gro
  apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm npt_${i} -cpt 15
  echo 18 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_${i}.edr -o pressure_${i}.xvg # choose Pressure (18), 0 terminates input
  echo 24 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_${i}.edr -o density_${i}.xvg # choose Density (24), 0 terminates input
  prev=npt_${i}.gro
done
python $scripts/plot_xvg.py pressure_*.xvg
python $scripts/plot_xvg.py density_*.xvg
cp npt_5.gro ../6_HREX/npt_5.gro
cp topol_5.top ../6_HREX/topol_5.top

echo "============= Setup completed successfully. ============="