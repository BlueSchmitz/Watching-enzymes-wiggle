#!/bin/bash  

#SBATCH -J optimization_setup_Ec5EKYm_apptainer  
#SBATCH -t 02:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --cpus-per-task 10
#SBATCH --gpus=0
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=./optimization_%j.out

# Exit immediately on errors, undefined vars, or failed pipes
set -euo pipefail
set -o errtrace

# Set paths for mdp_templates, force_fields, scripts, and project directory
export GMXLIB=$TMPDIR/MD_simulations/force_fields
mdp=$TMPDIR/MD_simulations/mdp_templates
scripts=$TMPDIR/MD_simulations/scripts
project_dir=5EKY_monomeric
export GROMACS_CONTAINER=$HOME/Blue/software/apptainer_2021/gromacs_plumed.sif # path to gromacs apptainer container

# Copy input files to scratch
cp -r $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/ "$TMPDIR"
cd $TMPDIR/MD_simulations/projects/$project_dir/outputs

# Function to copy back results on exit
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

### Tuning parallelization for HREX MD runs ###
# one node has 128 cores (rome), we have 12 replicas --> max. 10 cores per replica
# possible combinations of MPI ranks and OpenMP threads per replica
# np = MPI ranks  ntomp = OpenMP threads per rank
# Possible MPI ranks per replica (np) = max. -n
np_list=(1 2 3 4 5 6 7 8 9 10)

# Candidate threads per rank (ntomp) = max. --cpus-per-task
ntomp_list=(1 2 3 4 5 6 7 8 9 10)

# Maximum cores per replica
max_cores_per_replica=10

mkdir -p ./6_HREX
cd ./6_HREX
# generate scaled topology files and tpr files for all replicas
for i in 1.00 0.95 0.91 0.87 0.83 0.79 0.76 0.72 0.69 0.66 0.63 0.60;
do 
  cd ./rep${i}
  apptainer exec $GROMACS_CONTAINER mpirun gmx_mpi grompp -f $mdp/optimization.mdp -c ../npt_5.gro -p scaled_${i}.top -o topol.tpr -maxwarn 1
  echo "Generated topol.tpr from scaled_${i}.top in ./rep${i}."
  cd ..
done 
echo "Generated all scaled topology files and tpr files for HREX replicas."
# 5. Empty plumed.dat file (no additional CVs, only HREX)
: > plumed.dat # empty plumed file

# Temporary file to store results
RESULTS_FILE="tune_summary.csv"
echo "MPI ranks (np),OpenMP threads per rank (ntomp),Total cores,Performance (ns/day),Performance (h/ns),Speed per core (ns/day/core)" > $RESULTS_FILE

# Loop over candidate combinations
for np in "${np_list[@]}"; do
    for nt in "${ntomp_list[@]}"; do
        total_cores=$((np * nt))
        if [ $total_cores -le $max_cores_per_replica ]; then
            echo "=== Testing np=$np ntomp=$nt (total cores=$total_cores) ==="
            totalranks=$((np * 12)) # 12 replicas
            export OMP_NUM_THREADS=$nt

            # Use -deffnm temporary output to avoid overwriting
            TMP_PREFIX="tune_np${np}_nt${nt}"
            echo "Running tuning with deffnm=$TMP_PREFIX"
            # Run HREX MD for a short duration to measure performance
            apptainer exec $GROMACS_CONTAINER mpirun -np $totalranks \
                gmx_mpi mdrun multidir rep* \
                -replex 1000 \
                -ntomp $nt \
                -f $mdp/optimization.mdp \
                -plumed ../plumed.dat \
                -cpt 15 \
                -deffnm $TMP_PREFIX 2>&1 | tee ${TMP_PREFIX}.log || true

            # Extract relevant metrics from log
            ns_per_day=$(grep "Performance:" ${TMP_PREFIX}.log | awk '{print $2}' || echo "NA")
            hours_per_ns=$(grep "Performance:" ${TMP_PREFIX}.log | awk '{print $3}' || echo "NA")
            speed_per_core=$(echo "$ns_per_day / $total_cores" | bc -l || echo "NA")

            # Append to results CSV
            echo "$np,$nt,$total_cores,$ns_per_day,$hours_per_ns,$speed_per_core" >> $RESULTS_FILE
        fi
    done
done

# Visualize results
python $scripts/analyse_optimization.py ./tune_summary.csv

echo "=== Tuning complete! Results saved in $RESULTS_FILE ==="
cat $RESULTS_FILE

echo "=== Copying project outputs back to home ==="
rsync -av \
      --exclude 'rleveson.*' \
      --exclude '*.out' \
      $TMPDIR/MD_simulations/projects/$project_dir/outputs/ \
      $HOME/Blue/Watching-enzymes-wiggle/MD_simulations/projects/$project_dir/outputs/
echo "=== Copy complete ==="