#!/bin/bash
#SBATCH --job-name="test"
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --qos=job_gpu
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=8G
#SBATCH --exclude=gnode23

module load CUDA
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate openmm
bash utils/run.sh '0,1'
