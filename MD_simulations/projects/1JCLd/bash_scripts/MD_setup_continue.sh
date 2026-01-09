#!/bin/bash  

#SBATCH -J 1JCLd_setup  
#SBATCH -t 20:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16 
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLd_setup_%j.out
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

Run this script from the projects/5EKY_monomeric/ directory, it contains relative paths
'
# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

# Set paths for mdp_templates, force_fields and pdb file (to change paths quickly)
export GMXLIB=$TMPDIR/MD_simulations/force_fields
project_dir=1JCLd
mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts
pdb=$TMPDIR/MD_simulations/projects/$project_dir/inputs/*.pdb
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif # path to gromacs apptainer container

# Copy input files to scratch
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/ "$TMPDIR"
cd $TMPDIR/MD_simulations/projects/$project_dir

# Function to copy back results when error occurs before the script exits
function copy_back_results {
    set +e # Disable exit on error for this function
    echo "=== Error occured. Copying results back to home. ==="
    if [[ -d "$TMPDIR/MD_simulations/projects/$project_dir/outputs" ]]; then
        rsync -av \
          --exclude 'rleveson.*' \
          --exclude '*.out' \
          "$TMPDIR/MD_simulations/projects/$project_dir/outputs/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/"
        echo "=== Copy complete ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
}

# Trap errors and print line number + command
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
trap copy_back_results EXIT

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
# Print modules 
module list

# mkdir outputs directories 
mkdir -p ./outputs/2_parametrization ./outputs/3_minimization ./outputs/4_equilibration ./outputs/5_sanity_checks ./outputs/6_HREX

# NPT Equilibration
cd ./outputs/4_equilibration
apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -deffnm npt_250 -cpt 15 -cpi npt_250.cpt
echo 18 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_250.edr -o pressure_250.xvg # choose Pressure (18), 0 terminates input
echo 24 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_250.edr -o density_250.xvg # choose Density (24), 0 terminates input

prev=npt_250.gro
for i in 100 5;
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
  apptainer exec $GROMACS_CONTAINER mpirun -np $OMP_NUM_TASKS gmx_mpi mdrun -deffnm npt_${i} -cpt 15 
  echo 18 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_${i}.edr -o pressure_${i}.xvg # choose Pressure (18), 0 terminates input
  echo 24 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_${i}.edr -o density_${i}.xvg # choose Density (24), 0 terminates input
  prev=npt_${i}.gro
done
python $scripts/plot_xvg.py pressure_*.xvg
python $scripts/plot_xvg.py density_*.xvg
cp npt_5.gro ../6_HREX/npt_5.gro
cp topol_5.top ../6_HREX/topol_5.top
cp topol_Protein_chain_A_5.itp ../6_HREX/topol_Protein_chain_A_5.itp
cp topol_Protein_chain_B_5.itp ../6_HREX/topol_Protein_chain_B_5.itp

echo "============= Setup completed successfully. ============="

echo "============= Copying project outputs back to home ============="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude '*.out' \
      $TMPDIR/MD_simulations/projects/$project_dir/outputs/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/
echo "============= Copy complete ============="
