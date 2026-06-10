#!/bin/bash

#SBATCH -J EcDERA_dump_frame
#SBATCH -t 10:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task=2
#SBATCH --gpus=0
#SBATCH --output=./EcDERA_dump_frame_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

# Exit immediately on errors
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=EcDERA_sub
output_dir=7_simple_MD
pH=7

export GMXLIB=$TMPDIR/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/inputs/*.pdb

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_TASKS=$SLURM_NTASKS

# Create temporary directory on scratch for this job
tmpdir=$(mktemp -d /gpfs/scratch1/shared/rleveson/Blue/tmp.XXXXXX)
mkdir -p "$tmpdir/MD_simulations/projects/"

### Load modules ###
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module list

### Copy project to scratch ###
echo "=== Copying project to scratch ==="
cp -r /gpfs/work1/0/prjs2080/$project_dir "$tmpdir/MD_simulations/projects/"
cd $tmpdir/MD_simulations/projects/$project_dir

# Function to copy back results when error occurs and before the script exits
function copy_back_results {
    set +e +u # Disable exit on error for this function
    echo "=== Copying results back to project folder at $(date). ==="
    if [[ -d "$tmpdir/MD_simulations/projects/$project_dir/outputs/$output_dir" ]]; then
        rsync -av --partial --inplace \
          "$tmpdir/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
          "/gpfs/work1/0/prjs2080/$project_dir/outputs/$output_dir/"
        echo "=== Copy complete at $(date) ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
    # Clean up temporary directory
    rm -rf "$tmpdir"
    echo "=== Temporary directory cleaned up at $(date). ==="
}
trap copy_back_results EXIT

mkdir -p ./outputs/$output_dir
cd ./outputs/$output_dir/

# Downsample and save trajectory
echo "============= Downsizing and Exporting trajectory ============="
# downsample trajectory 
#echo 1 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s md.tpr -f md.xtc -o md_center_mol.xtc -center -pbc mol -ur compact
# fit trajectory to reference (TIM barrel backbone) to remove overall rotation and translation, keep whole system
#echo 4 22 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s md.tpr -f md_center_mol.xtc -o md_fit.xtc -fit rot+trans -n index.ndx
#rm md_center_mol.xtc
#echo -e "1\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s topol.tpr -n index.ndx -o md_1000.xtc -dt 1000
#echo "Trajectory saved as md_1000.xtc"

# extract closed conformations
#echo -e "0\n0" | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_closed_rep3.xtc -s ./rep3/md3.tpr -o closed_start.gro -n index.ndx -dump 9000

### 3 Energy minimization ###
#apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/minim.mdp -c closed_start.gro -p topol.top -o em.tpr
#apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm em
#echo 10 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f em.edr -o potential.xvg # choose potential energy (10), 0 terminates input
#python $scripts/plot_xvg.py potential.xvg

### 4 Equilibration ###
echo "============= Equilibration with GROMACS ============="
# NVT Equilibration
# Restraint file
#echo -e "q" | apptainer exec $GROMACS_CONTAINER gmx_mpi make_ndx -f em.tpr -o index.ndx << EOF
#1 | 13
#name 22 Protein_KPS

#q
#EOF

#echo 22 | apptainer exec $GROMACS_CONTAINER gmx_mpi genrestr -f em.gro -o posre.itp -n index.ndx

### include posre.itp in the topology file!

#apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
#apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm nvt -cpt 15
#echo 16 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f nvt.edr -o temperature.xvg # choose Temperature (16), 0 terminates input
#python $scripts/plot_xvg.py temperature.xvg
# Preparation for NPT Equilibration
# Gradually reduce restraints from 1000 to 5 kJ mol−1 nm−2 by running 5 short NPT simulations of 500 ps each (5*500=2.5 ns)
#for i in 1000 500 250 100 5;
#do
  # Copy posre.itp 5 times and modify the force constant in each file
#  cp posre.itp posre_$i.itp
#  sed -i "s/\b1000\b/$i/g" posre_$i.itp # \b for whole word match
#  echo "Modified posre.itp to posre_$i.itp."
  # Modify topol.top to include the correct posre file for each run
#  cp topol.top topol_$i.top
#  sed -i "s/posre.itp/posre_$i.itp/g" topol_$i.top
#  echo "Modified topol.top to topol_$i.top."
#done

### fix NPT topol files 

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

echo "============= Setup completed successfully. ============="