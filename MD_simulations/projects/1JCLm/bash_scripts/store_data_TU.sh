#!/bin/bash  

#SBATCH -J copy_to_TU
#SBATCH -t 10:00:00
#SBATCH -p rome
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blueschmitz@tudelft.nl
#SBATCH --output=copy_to_TU_%j.out

SRC="./outputs/*"
DEST="tudelft_sftp:/staff-umbrella/biocat dera"

rclone copy "$SRC" "$DEST" \
    -P \
    --transfers=8 \
    --checkers=16 \
    --timeout=1h \
    --retries=10 \
    --retries-sleep=30s \
    --log-file=rclone_${SLURM_JOB_ID}.log \
    --log-level=INFO