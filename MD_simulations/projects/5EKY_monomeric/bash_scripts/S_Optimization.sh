#!/bin/bash  

#SBATCH -J optimization_setup_Ec5EKYm_apptainer  
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --cpus-per-task 10
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./outputs/optimization_%j.out

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

function copy_back_results {
    set +e # Disable exit on error for this function
    echo "=== Copying results back to home ==="
    if [[ -d "$TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs" ]]; then
        rsync -av \
          --exclude 'rleveson.*' \
          --exclude 'sanity_checks_Ec5EKYm*.out' \
          "$TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/" \
          "$HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/"
        echo "=== Copy complete ==="
    else
        echo "Nothing to copy back (outputs directory not found)"
    fi
}

# Trap errors and print line number + command
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
trap copy_back_results EXIT

# path to gromacs apptainer container
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif

# Load modules:  
module load 2023

# Copy input files to scratch
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/ "$TMPDIR"
cd $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs

# Set paths for mdp_templates, force_fields and pdb file (to change quickly)
export GMXLIB=$TMPDIR/MD_simulations/force_fields # make sure this is correct
mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts

### Tuning parallelization for HREX MD runs ###
# one node has 128 cores (rome), we have 12 replicas --> max. 10 cores per replica
# possible combinations of MPI ranks and OpenMP threads per replica
# np = MPI ranks  ntomp = OpenMP threads per rank
# 1                       10
# 2                        5
# 3                        3
# 4                        2
# 5                        2
# 6                        1
# 7                        1
# 8                        1
# 9                        1
# 10                       1
# Possible MPI ranks per replica (np) = max. -n
np_list=(1 2 3 4 5 6 7 8 9 10)

# Candidate threads per rank (ntomp) = max. --cpus-per-task
ntomp_list=(1 2 3 4 5 6 7 8 9 10)

# Maximum cores per replica
max_cores_per_replica=10

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
            echo "Running tuning with deffnm=$TMP_PREFIX"
            apptainer exec $GROMACS_CONTAINER mpirun -np $np \
                gmx_mpi mdrun \
                -s scaled_1.00.tpr \
                -ntomp $nt \
                -dlb no \
                -tunepme \
                -deffnm $TMP_PREFIX 2>&1 | tee ${TMP_PREFIX}.log || true
            #apptainer exec $GROMACS_CONTAINER gmx_mpi tune_pme -ntmpi $np -ntomp $nt -s scaled_1.00.tpr -mdrun "gmx_mpi mdrun" -dlb no -deffnm $TMP_PREFIX 2>&1 | tee ${TMP_PREFIX}.log || true

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

echo "=== Copying project outputs back to home ==="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude 'sanity_checks_Ec5EKYm*.out' \
      $TMPDIR/MD_simulations/projects/5EKY_monomeric/outputs/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/
echo "=== Copy complete ==="