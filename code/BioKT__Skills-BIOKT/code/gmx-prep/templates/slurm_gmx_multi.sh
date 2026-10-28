#!/bin/bash
#SBATCH --job-name=HREMD
#SBATCH --partition=PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=NREPLICAS        # e.g. 8 for 8 replicas
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:NREPLICAS
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.out

module load gromacs/2024.1-mpi-gpu

NREP=NREPLICAS
PREFIX="runs/hremd/hremd_rep"

mpirun -np ${NREP} gmx_mpi mdrun \
    -v \
    -s "${PREFIX}" \
    -deffnm "${PREFIX}" \
    -multi ${NREP} \
    -replex 1000 \
    -ntomp "${SLURM_CPUS_PER_TASK}"
