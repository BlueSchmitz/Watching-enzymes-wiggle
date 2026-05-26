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

SRC="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/BtDERA/outputs/6_HREX/analysis"
#SRC="/gpfs/work1/0/prjs2080/BtDERA/outputs/6_HREX/rep1.00"
DEST="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/BtDERA/outputs/6_HREX/analysis"
#DEST="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/BtDERA/outputs/6_HREX/rep1.00"

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

SRC2="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/BbDERA/outputs/6_HREX/analysis"
DEST2="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/BbDERA/outputs/6_HREX/analysis"

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


SRC3="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/CbDERA/outputs/6_HREX/analysis"
DEST3="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/CbDERA/outputs/6_HREX/analysis"

echo "HOSTNAME: $(hostname)"
echo "PATH: $PATH"
which rclone || echo "RCLONE NOT FOUND"
rclone version || echo "RCLONE NOT RUNNING"

rclone copy "$SRC3" "$DEST3" \
    -P \
    -vv \
    --transfers=8 \
    --checkers=16 \
    --timeout=1h \
    --retries=10 \
    --retries-sleep=30s \
    --log-file=rclone_${SLURM_JOB_ID}.log \

SRC4="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/MlDERA/outputs/6_HREX/analysis"
DEST4="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/MlDERA/outputs/6_HREX/analysis"

echo "HOSTNAME: $(hostname)"
echo "PATH: $PATH"
which rclone || echo "RCLONE NOT FOUND"
rclone version || echo "RCLONE NOT RUNNING"

rclone copy "$SRC4" "$DEST4" \
    -P \
    -vv \
    --transfers=8 \
    --checkers=16 \
    --timeout=1h \
    --retries=10 \
    --retries-sleep=30s \
    --log-file=rclone_${SLURM_JOB_ID}.log \

SRC5="/home/rleveson/Blue/Watching-enzymes-wiggle/MD_simulations/projects/1JCLm/outputs/6_HREX/analysis"
DEST5="tudelft_sftp:/staff-umbrella/biocat dera/MD_simulations/1JCLm/outputs/6_HREX/analysis"

echo "HOSTNAME: $(hostname)"
echo "PATH: $PATH"
which rclone || echo "RCLONE NOT FOUND"
rclone version || echo "RCLONE NOT RUNNING"

rclone copy "$SRC5" "$DEST5" \
    -P \
    -vv \
    --transfers=8 \
    --checkers=16 \
    --timeout=1h \
    --retries=10 \
    --retries-sleep=30s \
    --log-file=rclone_${SLURM_JOB_ID}.log \