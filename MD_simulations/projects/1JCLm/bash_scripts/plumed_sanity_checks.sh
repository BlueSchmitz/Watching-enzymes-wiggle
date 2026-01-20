#!/bin/bash  

#SBATCH -J 1JCLm_plumed_sanity_checks
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 16 
#SBATCH --cpus-per-task 1
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --output=./1JCLm_plumed_%j.out
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

Run this script from the projects/$project_dir/ directory, it contains relative paths
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

# mkdir outputs
mkdir -p ./outputs/5_sanity_checks ./outputs/6_HREX

### 6 Set up HREX-MD with PLUMED ###
echo "============= Scale topologies with PLUMED =============" 
cd ./outputs/6_HREX
# Files we need: npt_5.gro as the starting structure, topol.top without pointers to posre.itp (no more restraints)
# 1. Remove constraints from topol.top by commenting out the line that includes posre_5.itp
cp topol_5.top topol_prod.top
if [ -f topol_Protein_chain_A_5.itp ]; then # check for dimer
    sed -i '/#ifdef POSRES_WATER/,/#endif/d' topol_prod.top
    sed -i '/#ifdef POSRES/,/#endif/d' topol_Protein_chain_A_5.itp
    sed -i '/#ifdef POSRES/,/#endif/d' topol_Protein_chain_B_5.itp
    echo "Removed position restraints from topol_prod.top and both chains."
else 
    sed -i '/#ifdef POSRES/,/#endif/d' topol_prod.top # Remove the protein restraint block (POSRES)
    sed -i '/#ifdef POSRES_WATER/,/#endif/d' topol_prod.top # Remove the water restraint block (POSRES_WATER)
    echo "Removed position restraints from topol_prod.top."
fi 
# 2. Generate a self-contained topology file
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/tempering.mdp -c npt_5.gro -p topol_prod.top -pp processed.top -o dummy.tpr -r npt_5.gro -maxwarn 2
echo "Generated self-contained processed.top without position restraints."
# 3. Edit the processed.top file to indicate which atoms we want to scale (marked with an _ after the residue name)
if [ -f topol_Protein_chain_A_5.itp ]; then # check for dimer
    python $scripts/scale_residues.py > processed_scaled_2loops.top
    python $scripts/scale_residues.py Protein_chain_A > processed_scaled_1loop.top
    echo "Generated processed_scaled_1loop.top and processed_scaled_2loops.top with selected residues for scaling."
else
    python $scripts/scale_residues.py > processed_scaled.top
    echo "Generated processed_scaled.top with selected residues for scaling."
fi
# 4. Scale the Hamiltonian of the selected atoms by the factors 1.00, 0.95, 0.91, 0.87, 0.83, 0.79, 0.76, 0.72, 0.69, 0.66, 0.63, and 0.60
if [ -f topol_Protein_chain_A_5.itp ]; then # check for dimer
    mkdir -p ./1loop ./2loops # scale either one or both loops
    for chains in 1loop 2loops;
    do 
        cd ./$chains
        : > plumed.dat # empty plumed file
        for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
        do 
          mkdir -p ./rep${i} # create directory for each replica
          apptainer exec $PLUMED_CONTAINER plumed partial_tempering ${i} < processed_scaled.top  > ./rep${i}/scaled_${i}.top
          echo "Generated scaled_${i}.top in ${chains}/rep${i} with scaling factor ${i}."
          cd ./rep${i}
          apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/hrex.mdp -c ../npt_5.gro -p scaled_${i}.top -o topol.tpr -maxwarn 1
          echo "Generated topol.tpr from scaled_${i}.top in ${chains}/rep${i}."
          cd ..
        done 
        cd ..
    done
    echo "Generated all scaled topology files and tpr files for HREX replicas in 1loop and 2loops."
else
    for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
    do 
      mkdir -p ./rep${i} # create directory for each replica
      apptainer exec $PLUMED_CONTAINER plumed partial_tempering ${i} < processed_scaled.top  > ./rep${i}/scaled_${i}.top
      echo "Generated scaled_${i}.top in ./rep${i} with scaling factor ${i}."
      cd ./rep${i}
      apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/hrex.mdp -c ../npt_5.gro -p scaled_${i}.top -o topol.tpr -maxwarn 1
      echo "Generated topol.tpr from scaled_${i}.top in ./rep${i}."
      cd ..
    done 
    : > plumed.dat # empty plumed file
    echo "Generated all scaled topology files and tpr files for HREX replicas."
fi

### 5 Sanity checks ### 
echo "Sanity checks"
cp npt_5.gro ../5_sanity_checks/npt_5.gro
cp processed.top ../5_sanity_checks/processed.top
if [ -f topol_Protein_chain_A_5.itp ]; then # check for dimer
    cp ./1loop/rep1.00/scaled_1.00.top ../5_sanity_checks/scaled_1loop_1.00.top
    cp ./2loops/rep1.00/scaled_1.00.top ../5_sanity_checks/scaled_2loops_1.00.top
else
    cp ./rep1.00/scaled_1.00.top ../5_sanity_checks/scaled_1.00.top
fi
cd ../5_sanity_checks
# 1. Sanity check of created topologies: compare energies between original and 1.00 scaled system
# produce tpr + short trajectory (use processed.top)
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p processed.top -o make_traj.tpr -maxwarn 1
apptainer exec $GROMACS_CONTAINER gmx_mpi mdrun -deffnm make_traj
cp make_traj.xtc traj.xtc
# produce tpr for original topology
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p processed.top -o topol_pro.tpr -maxwarn 1
# compute energies using original
apptainer exec $GROMACS_CONTAINER gmx_mpi mdrun -rerun traj.xtc -s topol_pro.tpr -e ener_pro.edr -g rerun_pro.log
# produce tpr from scaled topology (use same mdp / coords)
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c npt_5.gro -p scaled_1.00.top -o scaled_1.00.tpr -maxwarn 1
# Recompute energies using the scaled topology and the previously made trajectory
apptainer exec $GROMACS_CONTAINER gmx_mpi mdrun -rerun traj.xtc -s scaled_1.00.tpr -e ener_scaled_1.00.edr -g rerun_scaled_1.00.log
# compare energies
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_pro.edr -o energies_pro.xvg -xvg none
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_scaled_1.00.edr -o energies_scaled_1.00.xvg -xvg none
python $scripts/energy_comparison.py energies_pro.xvg energies_scaled_1.00.xvg > energy_diff_1.00.log
# 2. Sanity check of scaled topologies: compare energies between original and 0.5 scaled system
python $scripts/scale_all_residues.py 0.5 > processed_scaled_0.5_all.top
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f $mdp/sanity_check.mdp -c ./npt_5.gro -p ./scaled_0.5_all.top -o scaled_0.5_all.tpr -maxwarn 2
apptainer exec $GROMACS_CONTAINER gmx_mpi mdrun -rerun traj.xtc -s scaled_0.5_all.tpr -e ener_scaled_0.5_all.edr -g rerun_scaled_0.5_all.log
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n10" | apptainer exec $GROMACS_CONTAINER gmx_mpi energy -f ener_scaled_0.5_all.edr -o energies_scaled_0.5_all.xvg -xvg none
python $scripts/energy_comparison.py energies_scaled_1.00.xvg energies_scaled_0.5_all.xvg > energy_diff_0.5_all.log
# 3. Sanity check of replica-exchange implementation
# Run a short HREX with two equivalent topology files (topol.tpr and scaled_1.00.tpr)
mkdir -p ./rep0 ./rep1
: > plumed.dat # empty plumed file
cp plumed.dat ./rep0/plumed.dat
cp plumed.dat ./rep1/plumed.dat
# copy mdp files and set different random seeds for both replicas
cp $mdp/sanity_check.mdp ./rep0/sanity_check_rep0.mdp
cp $mdp/sanity_check.mdp ./rep1/sanity_check_rep1.mdp
sed -i 's/^gen_seed.*/gen_seed = 12345/' ./rep0/sanity_check_rep0.mdp
sed -i 's/^gen_seed.*/gen_seed = 67890/' ./rep1/sanity_check_rep1.mdp
# Prepare tpr for rep0 (original topology)
cd rep0
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f sanity_check_rep0.mdp -c ../npt_5.gro -p ../processed.top -o topol.tpr -maxwarn 1
# Prepare tpr for rep1 (scaled topology)
cd ../rep1
apptainer exec $GROMACS_CONTAINER gmx_mpi grompp -f ./sanity_check_rep1.mdp -c ../npt_5.gro -p ../scaled_1.00.top -o topol.tpr -maxwarn 1
cd ..
# Run a short HREX simulation with 2 replicas
apptainer exec $GROMACS_CONTAINER mpirun -np $SLURM_NTASKS gmx_mpi mdrun -multidir rep0 rep1 -replex 200 -hrex -plumed plumed.dat
# Inspect logs
echo "Replica-exchange acceptance"
grep "Repl" rep0/md.log

echo "============= Plumed sanity checks completed successfully. ============="

echo "============= Copying project outputs back to home ============="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude '*.out' \
      $TMPDIR/MD_simulations/projects/$project_dir/outputs/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/
echo "============= Copy complete ============="