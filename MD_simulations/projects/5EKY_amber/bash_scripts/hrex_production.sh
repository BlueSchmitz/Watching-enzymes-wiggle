#!/bin/bash  

#SBATCH --job-name="hrex production run monomer 5EKY"   
#SBATCH --time=48:00:00
#SBATCH --ntasks=12 
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=0
#SBATCH --partition=compute-p2
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=Research-AS-BN
#SBATCH --output=/scratch/blueschmitz/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/5EKYm_hrex_%j.out
#SBATCH --mail-type=ALL ##you can also set BEGIN/END
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

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
module load 2024r1
module load openmpi
# Add your Spack MPI module path
module use /home/blueschmitz/spack/share/spack/lmod/linux-rhel8-x86_64/openmpi/4.1.6-w6w5qi5/Core
# Load the modules
module load gromacs/2020.6
module load plumed/2.8.0


# Print modules 
module list

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
export GMXLIB=../../../../force_fields # make sure this is correct
mdp=../../../../mdp_templates
scripts=../../../../scripts

echo "Production HREX-MD with GROMACS and PLUMED"
srun gmx_mpi mdrun -multidir rep* -replex 2000 -f $mdp/hrex.mdp -plumed ../plumed.dat -cpt 15 # 4 ps/0.002 ps = exchange every 2000 steps