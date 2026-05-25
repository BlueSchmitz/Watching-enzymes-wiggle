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

#SRC="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/BtDERA"
#SRC="/gpfs/work1/0/prjs2080/BtDERA/outputs/6_HREX/rep1.00"
#DEST="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/BtDERA"
#DEST="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/BtDERA/outputs/6_HREX/rep1.00"

#echo "HOSTNAME: $(hostname)"
#echo "PATH: $PATH"
#which rclone || echo "RCLONE NOT FOUND"
#rclone version || echo "RCLONE NOT RUNNING"

#rclone copy "$SRC" "$DEST" \
#    -P \
#    -vv \
#    --transfers=8 \
#    --checkers=16 \
#    --timeout=1h \
#    --retries=10 \
#    --retries-sleep=30s \
#    --log-file=rclone_${SLURM_JOB_ID}.log \

SRC2="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/BtDERA"
DEST2="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/BtDERA"

echo "HOSTNAME: $(hostname)"
echo "PATH: $PATH"
which rclone || echo "RCLONE NOT FOUND"
rclone version || echo "RCLONE NOT RUNNING"

rclone copy "$SRC2" "$DEST2" \
    -P \
    -vv \
    --transfers=8 \
    --checkers=16 \
    --timeout=1h \
    --retries=10 \
    --retries-sleep=30s \
    --log-file=rclone_${SLURM_JOB_ID}.log \


#SRC3="/gpfs/work1/0/prjs2080/BtDERA/outputs/6_HREX"
#DEST3="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/BtDERA/outputs/6_HREX"

#rclone copy "$SRC3" "$DEST3" \
#    -P \
#    -vv \
#    --transfers=8 \
#    --checkers=16 \
#    --timeout=1h \
#    --retries=10 \
#    --retries-sleep=30s \
#    --log-file=rclone_${SLURM_JOB_ID}.log \