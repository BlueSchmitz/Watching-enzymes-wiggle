#!/bin/bash  

#SBATCH -J setup_Ec5EKY_amber_apptainer  
#SBATCH -t 12:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 2 
#SBATCH --cpus-per-task 8
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./outputs/setup_Ec5EKY_amber_apptainer_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl

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
set -o errtrace

# path to gromacs apptainer container
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif

# Load modules:  
module load 2023
module load matplotlib/3.7.2-gfbf-2023a


# Print modules 
module list

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
#export GMXLIB=../../../../force_fields # make sure this is correct
mdp=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/mdp_templates
pdb=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_amber/inputs # Input PDB file (with correct protonation states)
scripts=$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/scripts

# NPT Equilibration
cd ./outputs/4_equilibration
for i in 1000 500 250 100 5;
do
  echo "Running NPT equilibration with restraints = ${i}"

  apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/npt.mdp \
             -c ${prev:-nvt.gro} \
             -r ${prev:-nvt.gro} \
             -p topol_${i}.top \
             -o npt_${i}.tpr \
              -maxwarn 1
  # ${prev:-nvt.gro} ensures the first run starts from NVT output, then continues from the last .gro
  apptainer exec $GROMACS_CONTAINER mpirun -np 2 gmx_mpi mdrun -deffnm npt_${i} -cpt 15 -cpi npt_${i}.cpt -append
  echo 18 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_${i}.edr -o pressure_${i}.xvg # choose Pressure (18), 0 terminates input
  echo 24 0 | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f npt_${i}.edr -o density_${i}.xvg # choose Density (24), 0 terminates input
  prev=npt_${i}.gro
done
python $scripts/plot_xvg.py pressure_*.xvg
python $scripts/plot_xvg.py density_*.xvg
cp npt_5.gro ../6_HREX/npt_5.gro
cp topol_5.top ../6_HREX/topol_5.top

echo "Setup completed successfully."
