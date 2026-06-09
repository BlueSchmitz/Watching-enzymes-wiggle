#!/bin/bash  

#SBATCH -J BbDERA_sub_MD_analysis
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=md_analysis_%j.out

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
project_dir=BbDERA_sub
output_dir=7_simple_MD
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

### Analysis ###
echo "============= Analysis of trajectory ============="

# make index file with default + custom groups
apptainer exec $GROMACS_CONTAINER gmx_mpi make_ndx -f ./rep1/md1.tpr -o index.ndx << EOF
1 | 13
name 22 Protein_KPS
r 1-212
name 23 TIM_barrel
r 1-212 & a CA
name 24 CA_TIM
4 & 23
name 25 TIM_barrel_backbone
r 213-221
name 26 tail
r 213-221 & a CA
name 27 CA_tail
4 & 26 
name 28 tail_backbone
r 151 & a N2
name 29 Lys151_NZ
r 221 & a OH 
name 30 Tyr221_OH
22 & a CA
name 31 Protein_CA
r 151 & a C5
name 32 KPS_CA
31 | 32
name 33 Protein_KPS_CA


q
EOF

# center the trajectory on the whole protein and remove PBC artifacts, output in compact format
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./rep${i}/md${i}.xtc"
#    echo 1 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s $tpr -f $xtc -o md_center_mol_rep${i}.xtc -center -pbc mol -ur compact
#done
# fit trajectory to reference (TIM barrel backbone) to remove overall rotation and translation, keep whole system
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_center_mol_rep${i}.xtc"
#    echo 25 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -s $tpr -f $xtc -o md_fit_rep${i}.xtc -fit rot+trans -n index.ndx 
#    rm md_center_mol_rep${i}.xtc
#done
# calculate RMSD of TIM barrel backbone over time, output in xvg format
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_fit_rep${i}.xtc"
#    echo 25 25 | apptainer exec $GROMACS_CONTAINER gmx_mpi rms -s $tpr -f $xtc -n index.ndx -o rmsd_tim_barrel_backbone_rep${i}.xvg -tu ns
#done
#python $scripts/plotxvg_reps.py rmsd_tim_barrel_backbone_rep1.xvg rmsd_tim_barrel_backbone_rep2.xvg rmsd_tim_barrel_backbone_rep3.xvg
# calculate RMSD of tail backbone over time, output in xvg format
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_fit_rep${i}.xtc"
#    echo 28 28 | apptainer exec $GROMACS_CONTAINER gmx_mpi rms -s $tpr -f $xtc -n index.ndx -o rmsd_tail_backbone_rep${i}.xvg -tu ns
#done
#python $scripts/plotxvg_reps.py rmsd_tail_backbone_rep1.xvg rmsd_tail_backbone_rep2.xvg rmsd_tail_backbone_rep3.xvg
# calculate distance between Lys151 NZ and Tyr221 OH over time, output in xvg format
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_fit_rep${i}.xtc"
#    apptainer exec $GROMACS_CONTAINER gmx_mpi distance -s $tpr -f $xtc -n index.ndx -oall dist_k151_y221_rep${i}.xvg -tu ns -select 'com of group 29 plus com of group 30'
#done
#python $scripts/plotxvg_reps.py dist_k151_y221_rep1.xvg dist_k151_y221_rep2.xvg dist_k151_y221_rep3.xvg
# histogram of the distance between Lys151 NZ and Tyr221 OH with bin width of 0.2 nm, output in xvg format
#for i in 1 2 3; do
#    apptainer exec $GROMACS_CONTAINER gmx_mpi analyze -f dist_k151_y221_rep${i}.xvg -dist dist_k151_y221_hist_rep${i}.xvg -bw 0.2
#done
#python $scripts/plotxvg_hist_reps.py dist_k151_y221_hist_rep1.xvg dist_k151_y221_hist_rep2.xvg dist_k151_y221_hist_rep3.xvg
# RSMF of CA atoms of the whole protein, output in xvg format
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_fit_rep${i}.xtc"
#    echo 33 | apptainer exec $GROMACS_CONTAINER gmx_mpi rmsf -f md_fit_rep${i}.xtc -s $tpr -o rmsf_Ca_rep${i}.xvg -n index.ndx -b 20000 -res  # start at 20 ns (time in ps)
#done
#python $scripts/plot_RMSF_red_Bb_reps.py rmsf_Ca_rep1.xvg rmsf_Ca_rep2.xvg rmsf_Ca_rep3.xvg
### define frames with distance between Lys151 NZ and Tyr221 OH < 0.6 nm as "closed" and >= 0.6 nm as "open", output in xvg format
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_fit_rep${i}.xtc"
#    # Compute distance time series (ps)
#    apptainer exec $GROMACS_CONTAINER gmx_mpi distance -s $tpr -f $xtc -n index.ndx -oall dist_k151_y221_rep${i}.xvg -tu ps -select 'com of group 29 plus com of group 30'
    # create new trajectory with selected frames where distance < 0.6 nm
#    echo 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi trjconv -f md_fit_rep${i}.xtc -s $tpr -o md_closed_rep${i}.xtc -drop dist_k151_y221_ps.xvg -dropover 0.6
#done

# h-bonds and hydrophobic contacts analysis with MDAnalysis
#python $scripts/contact_matrices_reps.py ./rep1/md1.tpr ./md_closed_rep1.xtc ./md_closed_rep2.xtc ./md_closed_rep3.xtc

# Find proton relay systems
#for i in 1 2 3; do
#    tpr="./rep${i}/md${i}.tpr"
#    xtc="./md_fit_rep${i}.xtc"
#    echo "Finding proton relay systems in replicate ${i}..."
#    python $scripts/scan_proton_relay_MD_speedup.py $tpr $xtc
#    echo "Proton relay analysis for replicate ${i} complete. Copying results..."
#    if [[ -f relay_candidates.csv ]]; then
#        cp relay_candidates.csv relay_candidates.csv_rep${i}.csv
#        rm relay_candidates.csv
#    fi
#    if [[ -f direct_candidates.csv ]]; then
#        cp direct_candidates.csv direct_candidates.csv_rep${i}.csv
#        rm direct_candidates.csv
#    fi
#    echo "Results for replicate ${i} copied."
#done

# Find proton relay systems
tpr="./rep3/md3.tpr"
xtc="./md_fit_rep3.xtc"
echo "Finding proton relay systems in replicate 3..."
python $scripts/scan_proton_relay_MD_speedup.py $tpr $xtc
echo "Proton relay analysis for replicate 3 complete. Copying results..."
if [[ -f relay_candidates.csv ]]; then
    cp relay_candidates.csv relay_candidates.csv_rep3.csv
    rm relay_candidates.csv
fi
if [[ -f direct_candidates.csv ]]; then
    cp direct_candidates.csv direct_candidates.csv_rep3.csv
    rm direct_candidates.csv
fi
echo "Results for replicate 3 copied."

# PCA
# Combine trajectories of all 3 replicates for PCA
#apptainer exec $GROMACS_CONTAINER gmx_mpi trjcat -f md_fit_rep1.xtc md_fit_rep2.xtc md_fit_rep3.xtc -o concat.xtc
#tpr="./rep1/md1.tpr"
# Compute covariance matrix
#echo 24 33 | apptainer exec $GROMACS_CONTAINER gmx_mpi covar -s $tpr -f concat.xtc -n index.ndx -b 20000 -o eigenvalues.xvg -v eigenvectors.trr
    # fit to CA TIM, covariance of whole protein CA (so side-chains do not contribute to covariance)
    # eigenvalues.xvg contains the eigenvalues (variance along each PC as mean square fluctuation captured by that PC in nm^2), eigenvectors.trr contains the eigenvectors (PCs) as a trajectory
# Project trajectory onto PCs
#echo 24 33 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f concat.xtc -s $tpr -n index.ndx -proj proj.xvg -first 1 -last 2
#echo 24 33 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f concat.xtc -s $tpr -n index.ndx -proj proj_20_pcs.xvg -first 1 -last 20
# Extract extreme projections along PC1 and PC2
#echo 24 33 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f concat.xtc -s $tpr -n index.ndx -extr extremes.pdb -first 1 -last 2
# Eigenvector components per atom (which residues dominate the motion)
#echo 33 | apptainer exec $GROMACS_CONTAINER gmx_mpi anaeig -v eigenvectors.trr -f concat.xtc -s $tpr -n index.ndx -rmsf PC_rmsf_per_atom.xvg -first 1 -last 2
#calculate distance between Lys151 NZ and Tyr221 OH for each frame of concatenated traj
#apptainer exec $GROMACS_CONTAINER gmx_mpi distance -s $tpr -f concat.xtc -n index.ndx -oall lys151_y221_distance.xvg -tu ns -select 'com of group 29 plus com of group 30'
#python $scripts/PCA_Bb.py proj.xvg eigenvalues.xvg lys151_y221_distance.xvg proj_20_pcs.xvg

# How many frames in closed trajectory vs full?
#for i in 1 2 3; do
#    apptainer exec $GROMACS_CONTAINER gmx_mpi check -f md_fit_rep${i}.xtc
#    apptainer exec $GROMACS_CONTAINER gmx_mpi check -f md_closed_rep${i}.xtc
#done

# Clustering 
#python $scripts/run_clustering_Ec.py $tpr md_fit.xtc
#python $scripts/plot_clustering.py
# Extract medoids
#python $scripts/extract_medoids.py

echo "Analysis complete. Results will be copied back to home directory."