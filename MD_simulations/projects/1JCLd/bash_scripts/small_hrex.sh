#!/bin/bash  

#SBATCH -J 1JCLd_small_HREX
#SBATCH -t 05:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 120
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLd_small_HREX_%j.out
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

Run this script from the projects/5EKY_monomeric/ directory, it contains relative paths.
This script is written for running inside an apptainer on the Snellius super computer.
'

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

### Project-specific settings ###
project_dir=1JCLd
pH=7
# Set paths for mdp_templates, force_fields and pdb file (to change paths quickly)
export GMXLIB=$TMPDIR/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif # path to gromacs apptainer container
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif # path to pdb2pqr apptainer container
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif # path to plumed apptainer container
mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts
pdb=$TMPDIR/MD_simulations/projects/$project_dir/inputs/*.pdb

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
mkdir -p ./outputs/8_small_HREX

echo "============= Small HREX-MD with GROMACS and PLUMED to check acceptance rate ============="
cd ./outputs/8_small_HREX
for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
do 
    cd ./rep${i}
    apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/small_hrex.mdp -c ../npt_5.gro -p scaled_${i}_2loops.top -o topol.tpr -maxwarn 1
    echo "Generated topol.tpr from scaled_${i}_2loops.top in rep${i}."
    cd ..
done
# 500 total exchanges
apptainer exec $GROMACS_CONTAINER mpirun -np 120 \
    gmx_mpi mdrun -multidir rep* \
    -replex 1000 \
    -plumed ../plumed.dat \
    -ntomp 1 \
    -hrex \
    -cpt 15
echo "============= Small HREX-MD complete ============="