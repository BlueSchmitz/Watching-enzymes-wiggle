#!/bin/bash  

#SBATCH -J copy_to_TU
#SBATCH -t 00:20:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=copy_to_TU_%j.out

SRC="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/5EKY_monomeric/outputs/"
DEST="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/EcDERA/Old_5EKY"

echo "HOSTNAME: $(hostname)"
echo "PATH: $PATH"
which rclone || echo "RCLONE NOT FOUND"
rclone version || echo "RCLONE NOT RUNNING"

rclone copy "$SRC" "$DEST" \
    -P \
    -vv \
    --transfers=8 \
    --checkers=16 \
    --timeout=1h \
    --retries=10 \
    --retries-sleep=30s \
    --log-file=rclone_${SLURM_JOB_ID}.log \