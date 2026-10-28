#!/bin/bash
#SBATCH --job-name=neo_rerun
#SBATCH --partition=GPU_long
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --array=1-8%4
#SBATCH --output=./logs/neo_rerun_%A_%a.out
#SBATCH --error=./logs/neo_rerun_%A_%a.err

# Targeted re-run of 8,231 isoforms (see rerun_isoform_ids.txt):
#   - 3,733 never submitted to the original run (clean AA, have a structure)
#   - 4,512 whose canonical basis moved in the 2026-07 gene re-attribution, so their
#     original "novel regions" were diffed against the WRONG canonical protein.
#
# Reads hpc_input_rerun.tsv and writes to output_rerun/ — a SEPARATE input and output dir.
# This is deliberate and load-bearing:
#   * output/ still holds the original 20 parquets + 20 .done files. Writing there would
#     overwrite them.
#   * hpc_predict_neoantigens.py SKIPS any isoform_id already present in a sibling .done
#     file. The 4,512 basis-changed isoforms appear in the ORIGINAL .done/parquets, so
#     reusing output/ would silently skip exactly the isoforms we are trying to recompute.
# Keep output_rerun/ free of stale .done files before submitting.

HPC_ROOT="/path/to/neoantigen_predictions"
ENV_PREFIX="${HPC_ROOT}/envs/neo"
TOOLS_DIR="${HPC_ROOT}/tools"
WORK_DIR="${HPC_ROOT}"

mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/output_rerun"

PYTHON="${ENV_PREFIX}/bin/python"
export PATH="${ENV_PREFIX}/bin:${PATH}"
export NETMHCPAN_BIN="${TOOLS_DIR}/netmhcpan-4.2/netMHCpan"
export BIGMHC_DIR="${TOOLS_DIR}/bigmhc"

cd "${WORK_DIR}"

echo "═══════════════════════════════════════════════════════════"
echo "  SLURM Array Job:  ${SLURM_ARRAY_JOB_ID}"
echo "  Task:             ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_COUNT}"
echo "  Node:             $(hostname)"
echo "  GPUs:             ${CUDA_VISIBLE_DEVICES:-none}"
echo "═══════════════════════════════════════════════════════════"

nvidia-smi 2>/dev/null || echo "(no GPU info)"

"${PYTHON}" hpc_predict_neoantigens.py \
    --input  hpc_input_rerun.tsv \
    --output "output_rerun/predictions_rerun_${SLURM_ARRAY_TASK_ID}.parquet" \
    --workers 1 \
    --batch-size 100 \
    --checkpoint-every 500

echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
