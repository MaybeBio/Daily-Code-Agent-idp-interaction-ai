#!/bin/bash
#SBATCH --job-name="test"
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --qos=job_gpu
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=8G


module load CUDA/11.8.0
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate openmm
python -m openmm.testInstallation
