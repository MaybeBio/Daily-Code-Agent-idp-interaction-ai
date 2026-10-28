#!/bin/bash
#SBATCH --job-name=Q3LI60_1_44_0a4
#SBATCH --account=OPEN-32-28
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=48:00:00
#SBATCH --partition=qgpu
#SBATCH --output=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/Q3LI60_1_44_0a470d24/Q3LI60_1_44_0a470d24.out
#SBATCH --error=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/Q3LI60_1_44_0a470d24/Q3LI60_1_44_0a470d24.err

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /home/shakespeare2/.conda/envs/calvados

cd /mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/Q3LI60_1_44_0a470d24
python run.py
