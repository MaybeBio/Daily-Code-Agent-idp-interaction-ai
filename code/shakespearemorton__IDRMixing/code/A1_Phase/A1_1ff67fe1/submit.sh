#!/bin/bash
#SBATCH --job-name=A1_1ff67fe1
#SBATCH --account=OPEN-32-28
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=48:00:00
#SBATCH --partition=qgpu
#SBATCH --output=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/A1_1ff67fe1/A1_1ff67fe1.out
#SBATCH --error=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/A1_1ff67fe1/A1_1ff67fe1.err

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /home/shakespeare2/.conda/envs/calvados

cd /mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/A1_1ff67fe1
python run.py
