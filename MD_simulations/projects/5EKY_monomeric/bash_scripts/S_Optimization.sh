#!/bin/bash  

#SBATCH -J optimization_setup_Ec5EKYm  
#SBATCH -t 00:30:00
#SBATCH -p genoa
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --cpus-per-task 10
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./outputs/optimization%j.out

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

function copy_back_results {
    echo "=== Copying results back to home ==="
    rsync -av \
          --exclude 'rleveson.*' \
          --exclude 'sanity_checks_Ec5EKYm*.out' \
        $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/ \
        $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/
    echo "=== Copy complete ==="
}

# Trap errors and print line number + command
trap 'echo "Error in ${BASH_SOURCE[0]} at line ${LINENO}: ${BASH_COMMAND}"' ERR
trap copy_back_results EXIT

# Load modules:  
module load 2024
module load mpicopy/4.2-gompi-2024a
module load Miniconda3/24.7.1-0

# Copy input files to scratch
mpicopy $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/
#echo "Contents of $TMPDIR:"
#tree -a -L 10 $TMPDIR
cd $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs
#cd $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs

# activate environment with GROMACS installed
unset CONDA_SHLVL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gromacs2020
# Force use of Conda MPI
export PATH=$CONDA_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v "$CONDA_PREFIX/lib" | paste -sd:)

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
export GMXLIB=$TMPDIR/MD_simulations/force_fields # make sure this is correct
project=$TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs
mdp=$TMPDIR/MD_simulations/mdp_templates
pdb=$TMPDIR/MD_simulations/projects/5EKY_monomeric/inputs # Input PDB file (with correct protonation states)
pyscripts=$TMPDIR/MD_simulations/scripts

# Candidate MPI ranks per replica (np) = max. -n
np_list=(1 2 4 8 10)

# Candidate threads per rank (ntomp) = max. --cpus-per-task
ntomp_list=(1 2 4 8 10)

# Maximum cores per replica
max_cores_per_replica=16

# Temporary file to store results
cd ./5_sanity_checks
RESULTS_FILE="tune_summary.csv"
echo "np,ntomp,total_cores,step_per_sec,ns_per_day,imbalance_percent,opt_npme" > $RESULTS_FILE

# Loop over candidate combinations
for np in "${np_list[@]}"; do
    for nt in "${ntomp_list[@]}"; do
        total_cores=$((np * nt))
        if [ $total_cores -le $max_cores_per_replica ]; then
            echo "=== Testing np=$np ntomp=$nt (total cores=$total_cores) ==="
            export OMP_NUM_THREADS=$nt

            # Run tune_pme or mdrun -tune (short)
            # Use -deffnm temporary output to avoid overwriting
            TMP_PREFIX="tune_np${np}_nt${nt}"
            mpirun -np $np gmx_mpi mdrun -s scaled_1.00.tpr -dlb no -tune pme -deffnm $TMP_PREFIX 2>&1 | tee ${TMP_PREFIX}.log

            # Extract relevant metrics from log
            step_per_sec=$(grep "Performance:" ${TMP_PREFIX}.log | awk '{print $3}')
            ns_per_day=$(grep "Performance:" ${TMP_PREFIX}.log | awk '{print $7}')
            imbalance=$(grep "Average load imbalance:" ${TMP_PREFIX}.log | awk '{print $5}')
            opt_npme=$(grep "npme:" ${TMP_PREFIX}.log | awk '{print $2}')

            # Append to results CSV
            echo "$np,$nt,$total_cores,$step_per_sec,$ns_per_day,$imbalance,$opt_npme" >> $RESULTS_FILE
        fi
    done
done

echo "=== Tuning complete! Results saved in $RESULTS_FILE ==="
cat $RESULTS_FILE

#echo "=== Copying project outputs back to home ==="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude 'sanity_checks_Ec5EKYm*.out' \
      $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/
#echo "=== Copy complete ==="