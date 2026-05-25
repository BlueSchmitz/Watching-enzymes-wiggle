#!/bin/bash  

#SBATCH -J MlDERA_HREX_analysis
#SBATCH -t 10:00:00
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
project_dir=MlDERA
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
module load Python/3.11.3-GCCcore-12.3.0
module load matplotlib/3.7.2-gfbf-2023a
module load scikit-learn/1.3.1-gfbf-2023a
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
}
trap copy_back_results EXIT

mkdir -p ./outputs/$output_dir
cd ./outputs/$output_dir/

### Analysis ###
echo "============= Analysis of trajectory ============="

# quick access
tpr="../rep1.00/topol.tpr"
xtc="../rep1.00/traj_comp.xtc"

# make index file with default + custom groups
apptainer exec $GROMACS_CONTAINER gmx_mpi make_ndx -f $tpr -o index.ndx << EOF
r 1-219
name 18 TIM_barrel
r 1-219 & a CA
name 19 CA_TIM
4 & 18 
name 20 TIM_barrel_backbone
r 220-229
name 21 tail
r 220-229 & a CA
name 22 CA_tail
4 & 21 
name 23 tail_backbone
r 166 & a NZ
name 24 Lys166_NZ
r 229 & a OD1
name 25 Asp229_OH

q
EOF

# center the trajectory on the whole protein and remove PBC artifacts, output in compact format
#echo 1 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s $tpr -f $xtc -o md_center_mol.xtc -center -pbc mol -ur compact
# fit trajectory to reference (TIM barrel backbone) to remove overall rotation and translation, keep whole system
#echo 20 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s $tpr -f md_center_mol.xtc -o md_fit.xtc -fit rot+trans -n index.ndx
#rm md_center_mol.xtc
# calculate RMSD of TIM barrel backbone over time, output in xvg format
#echo 20 20 | apptainer exec $GROMACS_CONTAINER gmx_mpi rms -s $tpr -f md_fit.xtc -n index.ndx -o rmsd_tim_barrel_backbone.xvg -tu ns
#python $scripts/plotxvg.py rmsd_tim_barrel_backbone.xvg 
# calculate RMSD of tail backbone over time, output in xvg format
#echo 23 23 | apptainer exec $GROMACS_CONTAINER gmx_mpi rms -s $tpr -f md_fit.xtc -n index.ndx -o rmsd_tail_backbone.xvg -tu ns
#python $scripts/plotxvg.py rmsd_tail_backbone.xvg 
# calculate distance between Lys166 NZ and Asp229 OH over time, output in xvg format
#apptainer exec $GROMACS_CONTAINER gmx_mpi distance -s $tpr -f md_fit.xtc -n index.ndx -oall lys166_asp229_distance.xvg -tu ns -select 'com of group 24 plus com of group 25'
#python $scripts/plotxvg.py lys166_asp229_distance.xvg
# histogram of the distance between Lys166 NZ and Asp229 OH with bin width of 0.2 nm, output in xvg format
#apptainer exec $GROMACS_CONTAINER gmx_mpi analyze -f lys166_asp229_distance.xvg -dist lys166_asp229_hist.xvg -bw 0.2
#python $scripts/plotxvg_hist.py lys166_asp229_hist.xvg
# RSMF of CA atoms of the whole protein, output in xvg format
#echo 3 | apptainer exec $GROMACS_CONTAINER gmx_mpi rmsf -f md_fit.xtc -s $tpr -o rmsf_Ca.xvg -n index.ndx -b 20000 -res  # start at 20 ns (time in ps)
#python $scripts/plot_RMSF_red.py rmsf_Ca.xvg
### define frames with distance between Lys166 NZ and Asp229 OH < 0.6 nm as "closed" and >= 0.6 nm as "open", output in xvg format
# Compute distance time series (ps)
#apptainer exec $GROMACS_CONTAINER gmx_mpi distance -f md_fit.xtc -s $tpr -n index.ndx -select 'com of group 24 plus com of group 25' -oall dist_k167_y259_ps.xvg
# Create new trajectory with selected frames where distance < 0.6 nm
#echo 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit.xtc -s $tpr -o md_closed.xtc -drop dist_k167_y259_ps.xvg -dropover 0.6

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

#python $scripts/PCA.py proj.xvg eigenvalues.xvg lys166_asp229_distance.xvg proj_20_pcs.xvg

# How many frames in closed trajectory vs full?
#apptainer exec $GROMACS_CONTAINER gmx_mpi check -f md_fit.xtc
#apptainer exec $GROMACS_CONTAINER gmx_mpi check -f md_closed.xtc

# Clustering 
#python $scripts/run_clustering_Ml.py $tpr md_fit.xtc
#python $scripts/plot_clustering.py
# Extract medoids
python $scripts/extract_medoids.py

# h-bonds and hydrophobic contacts analysis with MDAnalysis
#python $scripts/contact_matrices_Ml.py $tpr md_closed.xtc
#python $scripts/contact_distributions_Ml.py $tpr md_closed.xtc

echo "Analysis complete. Results will be copied back to home directory."