#!/bin/bash
#SBATCH --job-name=A1_da2d9e14
#SBATCH --account=OPEN-32-28
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=48:00:00
#SBATCH --partition=qgpu
#SBATCH --output=/scratch/project/open-32-28/Extra_Jobs/A1_da2d9e14/A1_da2d9e14.out
#SBATCH --error=/scratch/project/open-32-28/Extra_Jobs/A1_da2d9e14/A1_da2d9e14.err

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate /home/shakespeare2/.conda/envs/calvados

cd /scratch/project/open-32-28/Extra_Jobs/A1_da2d9e14
python run.py
