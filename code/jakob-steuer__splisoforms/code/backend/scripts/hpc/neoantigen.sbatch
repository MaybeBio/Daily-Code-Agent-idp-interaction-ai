#!/bin/bash
#SBATCH --job-name=neoantigens
#SBATCH --partition=GPU_long
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --array=1-20%3
#SBATCH --output=./logs/neoantigens_%A_%a.out
#SBATCH --error=./logs/neoantigens_%A_%a.err

# ── Paths ────────────────────────────────────────────────────────────────────
HPC_ROOT="/path/to/neoantigen_predictions"
ENV_PREFIX="${HPC_ROOT}/envs/neo"
TOOLS_DIR="${HPC_ROOT}/tools"
WORK_DIR="${HPC_ROOT}"

mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/output"

# ── Env setup (use python binary directly, no conda activate needed) ───────
PYTHON="${ENV_PREFIX}/bin/python"
export PATH="${ENV_PREFIX}/bin:${PATH}"
export NETMHCPAN_BIN="${TOOLS_DIR}/netmhcpan-4.2/netMHCpan"
export BIGMHC_DIR="${TOOLS_DIR}/bigmhc"

# ── Run ─────────────────────────────────────────────────────────────────────
cd "${WORK_DIR}"

echo "═══════════════════════════════════════════════════════════"
echo "  SLURM Array Job:  ${SLURM_ARRAY_JOB_ID}"
echo "  Task:             ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_COUNT}"
echo "  Node:             $(hostname)"
echo "  GPUs:             ${CUDA_VISIBLE_DEVICES:-none}"
echo "  CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "═══════════════════════════════════════════════════════════"

nvidia-smi 2>/dev/null || echo "(no GPU info)"

# Each array task processes its slice of the input and writes a separate parquet.
# --workers 1: one process per GPU avoids BigMHC/MHCflurry GPU contention.
# --batch-size 100: BigMHC loads once per batch (not per isoform).
# --checkpoint-every 500: write parquet every 500 isoforms (survives wall-time kills).
# The script auto-reads SLURM_ARRAY_TASK_ID and SLURM_ARRAY_TASK_COUNT.
"${PYTHON}" hpc_predict_neoantigens.py \
    --input  hpc_input.tsv \
    --output "output/predictions_${SLURM_ARRAY_TASK_ID}.parquet" \
    --workers 1 \
    --batch-size 100 \
    --checkpoint-every 500

echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
