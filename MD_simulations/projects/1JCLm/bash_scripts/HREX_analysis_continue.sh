#!/bin/bash  

#SBATCH -J EcDERA_HREX_analysis
#SBATCH -t 00:30:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=hrex_analysis_%j.out

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

### Project-specific settings ###
project_dir=1JCLm
output_dir=6_HREX/analysis
pH=7

export GMXLIB=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/force_fields
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif
export PDB2PQR_CONTAINER=$HOME/Blue/software/apptainer_pdb2pqr/pdb2pqr.sif
export PLUMED_CONTAINER=$HOME/Blue/software/apptainer_plumed/plumed.sif

mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/inputs/*.pdb

# Create temporary directory on scratch for this job
tmpdir=$(mktemp -d /gpfs/scratch1/shared/rleveson/Blue/tmp.XXXXXX)
mkdir -p "$tmpdir/MD_simulations/projects/"

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a
module load Python/3.11.3-GCCcore-12.3.0
module list

### Copy project to scratch ###
echo "=== Copying project to scratch ==="
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir "$tmpdir/MD_simulations/projects/"
cd $tmpdir/MD_simulations/projects/$project_dir

# Function to copy back results when error occurs and before the script exits
function copy_back_results {
    set +e +u # Disable exit on error for this function
    echo "=== Copying results back to home at $(date). ==="
    if [[ -d "$tmpdir/MD_simulations/projects/$project_dir/outputs/$output_dir" ]]; then
        rsync -av --partial --inplace \
          "$tmpdir/MD_simulations/projects/$project_dir/outputs/$output_dir/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/$output_dir/"
        echo "=== Copy complete at $(date) ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
    # Clean up temporary directory
    rm -rf "$tmpdir"
    echo "Temporary directory cleaned up."
}
trap copy_back_results EXIT

mkdir -p ./outputs/$output_dir
cd ./outputs/$output_dir/

### Analysis ###
echo "============= Analysis of trajectory ============="

# quick access
tpr="../rep1.00/topol.tpr"
xtc="../rep1.00/traj_comp.xtc"

# Continue 
# h-bonds and hydrophobic contacts analysis with MDAnalysis
#python $scripts/contact_matrices.py $tpr md_closed.xtc
#python $scripts/plot_RMSF_red.py rmsf_Ca.xvg

# PCA
# Compute covariance matrix
#echo 19 3 | apptainer exec $GROMACS_CONTAINER gmx_mpi covar -s $tpr -f md_fit.xtc -n index.ndx -b 20000 -o eigenvalues.xvg -v eigenvectors.trr
    # fit to CA TIM, covariance of whole protein CA (so side-chains do not contribute to covariance)
    # eigenvalues.xvg contains the eigenvalues (variance along each PC as mean square fluctuation captured by that PC in nm^2), eigenvectors.trr contains the eigenvectors (PCs) as a trajectory
# Project trajectory onto PCs
#echo 19 3 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f md_fit.xtc -s $tpr -n index.ndx -proj proj.xvg -first 1 -last 2
#echo 19 3 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f md_fit.xtc -s $tpr -n index.ndx -proj proj_20_pcs.xvg -first 1 -last 20
# Extract extreme projections along PC1 and PC2
#echo 19 3 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f md_fit.xtc -s $tpr -n index.ndx -extr extremes.pdb -first 1 -last 2
# Eigenvector components per atom (which residues dominate the motion)
#echo 3 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f md_fit.xtc -s $tpr -n index.ndx -rmsf PC_rmsf_per_atom.xvg -first 1 -last 2

python $scripts/PCA_try.py proj.xvg eigenvalues.xvg lys167_tyr259_distance.xvg proj_20_pcs.xvg

# Extract extreme projections (from python script)
#min_pc1=$(awk '/min_pc1/ {print $2}' pc_extreme_frames.dat)
#echo 1 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s $tpr -dump $min_pc1 -o min_pc1.pdb
#max_pc1=$(awk '/max_pc1/ {print $2}' pc_extreme_frames.dat)
#echo 1 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s $tpr -dump $max_pc1 -o max_pc1.pdb
#min_pc2=$(awk '/min_pc2/ {print $2}' pc_extreme_frames.dat)
#echo 1 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s $tpr -dump $min_pc2 -o min_pc2.pdb
#max_pc2=$(awk '/max_pc2/ {print $2}' pc_extreme_frames.dat)
#echo 1 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s $tpr -dump $max_pc2 -o max_pc2.pdb

echo "Analysis complete. Results will be copied back to home directory."