#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --partition=PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

# ── adjust module name to your cluster ──────────────────────────────────────
module load gromacs/2024.1-gpu

# ── paths ───────────────────────────────────────────────────────────────────
TPR="runs/prod/prod.tpr"
DEFFNM="runs/prod/prod"

mkdir -p runs/prod

# ── single-GPU mdrun ─────────────────────────────────────────────────────────
gmx mdrun \
    -v \
    -s "${TPR}" \
    -deffnm "${DEFFNM}" \
    -ntmpi 1 \
    -ntomp "${SLURM_CPUS_PER_TASK}" \
    -pme gpu \
    -bonded gpu

echo "Done. Check ${DEFFNM}.log for performance summary."
