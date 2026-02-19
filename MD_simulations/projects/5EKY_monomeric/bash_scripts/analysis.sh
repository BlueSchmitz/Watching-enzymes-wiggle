#!/bin/bash  

#SBATCH --job-name="analyse 5EKYm trajectory"   
#SBATCH --time=01:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu-a100
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=Research-AS-BN
#SBATCH --output=/scratch/blueschmitz/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/5EKYm_analysis_%j.out
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
mdp=../../../../mdp_templates
scripts=../../../../scripts
pdb=../../inputs/5EKY_fill.BL00440001.pdb # Input PDB file (with correct protonation states)

# cd to outputs
cd ./outputs/7_simple_MD

### Analysis ###
echo "============= Analysis of trajectory ============="
# make index file with default + custom groups
gmx_mpi make_ndx -f md.tpr -o index.ndx << EOF
r 1-248
name 18 TIM_barrel
r 1-248 & a CA
name 19 CA_TIM
4 & 18 
name 20 TIM_barrel_backbone
r 249-259
name 21 tail
r 249-259 & a CA
name 22 CA_tail
4 & 21 
name 23 tail_backbone
r 167 & a NZ
name 24 Lys167_NZ
r 259 & a OH
name 25 Tyr259_OH

q
EOF

# center the trajectory on the whole protein and remove PBC artifacts, output in compact format
echo 1 0 | gmx_mpi trjconv -s md.tpr -f md.xtc -o md_center_mol.xtc -center -pbc mol -ur compact
# fit trajectory to reference (TIM barrel backbone) to remove overall rotation and translation, keep whole system
echo 20 0 | gmx_mpi trjconv -s md.tpr -f md_center_mol.xtc -o md_fit.xtc -fit rot+trans -n index.ndx
# calculate RMSD of TIM barrel backbone over time, output in xvg format
echo 20 20 | gmx_mpi rms -s md.tpr -f md_fit.xtc -n index.ndx -o rmsd_tim_barrel_backbone.xvg -tu ns
python $scripts/plotxvg.py -i rmsd_tim_barrel_backbone.xvg 
# calculate RMSD of tail backbone over time, output in xvg format
echo 23 23 | gmx_mpi rms -s md.tpr -f md_fit.xtc -n index.ndx -o rmsd_tail_backbone.xvg -tu ns
python $scripts/plotxvg.py -i rmsd_tail_backbone.xvg 
# calculate distance between Lys167 NZ and Tyr259 OH over time, output in xvg format
gmx_mpi distance -s md.tpr -f md_fit.xtc -n index.ndx -oall lys167_tyr259_distance.xvg -tu ns -select 'com of group 24 plus com of group 25'
python $scripts/plotxvg.py -i lys167_tyr259_distance.xvg
# histogram of the distance between Lys167 NZ and Tyr259 OH with bin width of 0.2 nm, output in xvg format
gmx_mpi analyze -f lys167_tyr259_distance.xvg -dist lys167_tyr259_hist.xvg -bw 0.2
# RSMF of CA atoms of the whole protein, output in xvg format
echo 3 | gmx_mpi rmsf -f md_fit.xtc -s md.tpr -o rmsf_Ca.xvg -n index.ndx -b 20000 -res  # start at 20 ns (time in ps)
python $scripts/plot_RMSF.py rmsf_Ca_10000.xvg
### define frames with distance between Lys167 NZ and Tyr259 OH < 0.6 nm as "closed" and >= 0.6 nm as "open", output in xvg format
# Compute distance time series (ps)
gmx_mpi distance -f md_fit.xtc -s md.tpr -n index.ndx -select 'com of group 24 plus com of group 25' -oall dist_k167_y259_ps.xvg
# Create new trajectory with selected frames where distance < 0.6 nm
echo 0 | gmx_mpi trjconv -f md_fit.xtc -s md.tpr -o md_closed.xtc -drop dist_k167_y259_ps.xvg -dropover 0.6
# number of h-bonds over time between protein and tail (DHA angle > 120 degrees --> HDA angle > 60 degrees, distance < 0.35 nm), output in xvg format
echo 2 21 | gmx_mpi hbond -s md.tpr -f md_closed.xtc -n index.ndx -num hbond_time_closed.xvg -tu ns -b 20000 -a 60 -r 0.25
# 
