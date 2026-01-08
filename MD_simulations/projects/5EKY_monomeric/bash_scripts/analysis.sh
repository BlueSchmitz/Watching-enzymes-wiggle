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
name 17 TIM_barrel
r 1-248 & a CA
name 18 CA_TIM
r 249-259
name 19 tail
r 249-259 & a CA
name 20 CA_loop
r 167 & a NZ
name 21 Lys167_NZ
r 259 & a OH
name 22 Tyr259_OH
q
EOF
# center protein and remove jumps, keep whole protein 
gmx_mpi trjconv -f md.xtc -s md.tpr -o md_nojump.xtc -pbc nojump -center
# Fit trajectory to reference (usually backbone or Cα)
gmx_mpi trjconv -s md.tpr -f md_nojump.xtc -o md_fit.xtc -fit rot+trans
