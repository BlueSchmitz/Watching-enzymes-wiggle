#!/bin/bash  

#SBATCH -J 1JCLm_minimal_lambda
#SBATCH -t 04:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLm_minimal_lambda_%j.out
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
project_dir=1JCLm
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
mkdir -p ./outputs/9_minimal_scaling

echo "============= Find lowest lambda required for HREX ============="
lambdas="1.0 0.9 0.8 0.7 0.6 0.5"
cp ./outputs/6_HREX/processed_scaled.top ./outputs/9_minimal_scaling/processed_scaled.top
cp ./outputs/6_HREX/npt_5.gro ./outputs/9_minimal_scaling/npt_5.gro
cd ./outputs/9_minimal_scaling
for l in $lambdas; do
    mkdir -p l$l
    cd l$l || exit 1

    echo "Running lambda $l"
    apptainer exec $PLUMED_CONTAINER plumed partial_tempering ${l} < ../processed_scaled.top  > scaled_${l}.top
    apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/small_hrex.mdp -c ../npt_5.gro -p scaled_${l}.top -o scaling_${l}.tpr -maxwarn 1
    apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm scaling_${l} -s scaling_${l}.tpr -cpt 15

    # generate plumed.dat
    echo "MOLINFO STRUCTURE=../../1_protonation/${project_dir}_pro.pdb" > plumed.dat
    echo "" >> plumed.dat

    args=""
    for r in $(seq 250 259); do
        echo "phi$r: TORSION ATOMS=@phi-$r" >> plumed.dat
        echo "psi$r: TORSION ATOMS=@psi-$r" >> plumed.dat
        args="$args,phi$r,psi$r"
    done

    args=${args#,}  # remove leading comma
    echo "" >> plumed.dat
    echo "PRINT ARG=$args FILE=COLVAR" >> plumed.dat

    # run plumed
    apptainer exec $PLUMED_CONTAINER plumed driver --ixtc scaling_${l}.xtc --plumed plumed.dat

    cd ../../
done

python $scripts/plot_torsions.py