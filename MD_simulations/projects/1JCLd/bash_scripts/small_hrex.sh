#!/bin/bash  

#SBATCH -J 1JCLd_small_HREX2
#SBATCH -t 20:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLd_small_HREX2_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --signal=B:USR1@300

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
output_dir=8_small_HREX/16lambdas
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
    echo "=== Error occured. Copying results back to home at $(date). ==="
    if [[ -d "$TMPDIR/MD_simulations/projects/$project_dir/outputs/8_small_HREX" ]]; then
        rsync -av \
          --exclude 'rleveson.*' \
          --exclude '*.out' \
          "$TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/"
        echo "=== Copy complete ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
}

# Trap errors and print line number + command
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
# Trap timeout warning (USR1)
trap 'echo "=== Walltime approaching, copying results ==="; copy_back_results' USR1
trap copy_back_results EXIT

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
# Print modules 
module list

# mkdir outputs directories 
mkdir -p ./outputs/$output_dir

echo "============= Small HREX-MD with GROMACS and PLUMED to check acceptance rate ============="
cd ./outputs/$output_dir
for i in 1.00 0.96 0.93 0.90 0.87 0.84 0.81 0.78 0.75 0.72 0.70 0.68 0.66 0.64 0.62 0.60;
do 
    mkdir -p ./rep${i} # create directory for each replica
    apptainer exec $PLUMED_CONTAINER plumed partial_tempering ${i} < ./processed_scaled_2loops.top  > ./rep${i}/scaled_${i}_2loops.top
    echo "Generated scaled_${i}.top in 2loops/rep${i} with scaling factor ${i}."
    cd ./rep${i}
    apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/small_hrex.mdp -c ../npt_5.gro -p scaled_${i}_2loops.top -o topol.tpr -maxwarn 1
    echo "Generated topol.tpr from scaled_${i}_2loops.top in rep${i}."
    cd ..
done
# 500 total exchanges
apptainer exec $GROMACS_CONTAINER mpirun -np 128 \
    gmx_mpi mdrun -multidir rep* \
    -replex 1000 \
    -plumed ../plumed.dat \
    -ntomp 1 \
    -hrex \
    -cpt 15 \
    -cpi state.cpt
echo "============= Small HREX-MD complete ============="

echo "============= Copying project outputs back to home ============="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude '*.out' \
      $TMPDIR/MD_simulations/projects/$project_dir/outputs/$output_dir/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/
echo "============= Copy complete ============="