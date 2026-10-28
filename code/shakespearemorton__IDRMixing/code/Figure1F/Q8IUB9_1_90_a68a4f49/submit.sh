#!/bin/bash
#SBATCH --job-name=Q8IUB9_1_90_a68
#SBATCH --account=OPEN-32-28
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=48:00:00
#SBATCH --partition=qgpu
#SBATCH --output=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/Q8IUB9_1_90_a68a4f49/Q8IUB9_1_90_a68a4f49.out
#SBATCH --error=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/Q8IUB9_1_90_a68a4f49/Q8IUB9_1_90_a68a4f49.err

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /home/shakespeare2/.conda/envs/calvados

cd /mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/Q8IUB9_1_90_a68a4f49
python run.py
