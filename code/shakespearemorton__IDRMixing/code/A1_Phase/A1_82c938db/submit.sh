#!/bin/bash
#SBATCH --job-name=A1_82c938db
#SBATCH --account=OPEN-32-28
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=48:00:00
#SBATCH --partition=qgpu
#SBATCH --output=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/A1_82c938db/A1_82c938db.out
#SBATCH --error=/mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/A1_82c938db/A1_82c938db.err

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /home/shakespeare2/.conda/envs/calvados

cd /mnt/storage-brno14-ceitec/shared/softmatter/shakespearem/Extra_Jobs/A1_82c938db
python run.py
