#!/bin/bash
#SBATCH --job-name=neo_v2
#SBATCH --partition=GPU_long
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --array=1-10%4
#SBATCH --output=./logs/neo_v2_%A_%a.out
#SBATCH --error=./logs/neo_v2_%A_%a.err

# v2 restage (2026-10-23): 27,593 isoforms = NetMHCpan-less shards of the 2026-07 run
# + isoforms whose canonical changed in the splisoforms_v2 rebuild. Input ids are
# splisoforms_v2 isoform ids → load with load_neoantigen_results.py --replace (no id map).
HPC_ROOT="/path/to/neoantigen_predictions"
ENV_PREFIX="${HPC_ROOT}/envs/neo"
TOOLS_DIR="${HPC_ROOT}/tools"
WORK_DIR="${HPC_ROOT}"
mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/output_v2"
PYTHON="${ENV_PREFIX}/bin/python"
export PATH="${ENV_PREFIX}/bin:${PATH}"          # also provides tcsh for the netMHCpan wrapper
export NETMHCPAN_BIN="${TOOLS_DIR}/netmhcpan-4.2/netMHCpan"
export BIGMHC_DIR="${TOOLS_DIR}/bigmhc"
export TMPDIR="${TMPDIR:-/tmp}"
cd "${WORK_DIR}"
echo "Array ${SLURM_ARRAY_JOB_ID} task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} node $(hostname) gpus ${CUDA_VISIBLE_DEVICES:-none}"
"${NETMHCPAN_BIN}" -version | head -1 || { echo "netMHCpan does not run on $(hostname)"; exit 2; }
"${PYTHON}" hpc_predict_neoantigens.py \
    --input  hpc_input_v2_restage_20261023.tsv \
    --output "output_v2/predictions_${SLURM_ARRAY_TASK_ID}.parquet" \
    --num-chunks 10 --workers 1 --batch-size 100 --checkpoint-every 500
echo "Task ${SLURM_ARRAY_TASK_ID} exit=$? finished at $(date)"
