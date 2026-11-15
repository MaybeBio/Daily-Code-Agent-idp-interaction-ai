#!/usr/bin/env bash
# Ablation of the IDP changes: 6 arms x 16 crops, 50 trajectories each, binder 130 aa.
# jobs.txt is tab separated: <settings>\t<advanced>\t<filters>, so each run carries its
# own arm config. Lines are grouped 16 per arm, so array task N == arm N.
# Packed 16 crops per node (4 per GPU, ~89 trajectories/node-hour).
#SBATCH --job-name=up500_full
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --account=ab011
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=4
#SBATCH --mem=700G
#SBATCH --array=1-7
#SBATCH --signal=B:USR1@120
#SBATCH --output=/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/upscale_500/logs/%A_%a.log

set -uo pipefail

ROOT=/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/upscale_500
BINDCRAFT_DIR=/capstor/store/cscs/swissai/ab011/jkuhn/Bindome-private/BindCraft

CROPS_PER_NODE=16
GPUS_PER_NODE=4
MEM_FRACTION=0.22
MAX_LEGS=4

mkdir -p "${ROOT}/logs"
THIS=$(realpath "$0")
LEG_FILE="${ROOT}/logs/.leg_up500_full_task${SLURM_ARRAY_TASK_ID}"
LEG=$(cat "${LEG_FILE}" 2>/dev/null || echo 1)
RESUBMITTED=0

resubmit() {
  [[ "${RESUBMITTED}" -eq 1 ]] && return 0
  RESUBMITTED=1
  if [[ "${LEG}" -ge "${MAX_LEGS}" ]]; then echo "[$(date -Is)] leg ${LEG}/${MAX_LEGS}, kein resubmit ($1)"; return 0; fi
  echo "$((LEG+1))" > "${LEG_FILE}"
  echo "[$(date -Is)] resubmit task ${SLURM_ARRAY_TASK_ID} leg $((LEG+1)) ($1)"
  sbatch --dependency=afterany:"${SLURM_JOB_ID}" --array="${SLURM_ARRAY_TASK_ID}" "${THIS}"
}
trap 'echo "[$(date -Is)] USR1"; resubmit walltime' USR1

FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * CROPS_PER_NODE + 1 ))
ARM=$(sed -n "${FIRST}p" "${ROOT}/jobs_full.txt" | cut -f2 | xargs basename | sed 's/^advanced_//; s/\.json$//')
echo "[$(date -Is)] task=${SLURM_ARRAY_TASK_ID} ARM=${ARM} leg=${LEG} crops ${FIRST}-$((FIRST+CROPS_PER_NODE-1)) host=$(hostname)"

srun --ntasks="${CROPS_PER_NODE}" --cpus-per-task=16 --kill-on-bad-exit=0 \
     --environment=jax_env \
     --output="${ROOT}/logs/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_task%t.log" \
     bash -c "
  set -e
  export CUDA_VISIBLE_DEVICES=\$(( \${SLURM_PROCID} % ${GPUS_PER_NODE} ))
  export XLA_PYTHON_CLIENT_PREALLOCATE=false
  export XLA_PYTHON_CLIENT_MEM_FRACTION=${MEM_FRACTION}
  export OMP_NUM_THREADS=4
  LINE=\$(( ${FIRST} + \${SLURM_PROCID} ))
  ROW=\$(sed -n \"\${LINE}p\" '${ROOT}/jobs_full.txt')
  if [[ -z \"\${ROW}\" ]]; then echo \"[task \${SLURM_PROCID}] keine Zeile \${LINE}\"; exit 0; fi
  SETTINGS=\$(printf '%s' \"\${ROW}\" | cut -f1)
  ADVANCED=\$(printf '%s' \"\${ROW}\" | cut -f2)
  FILTERS=\$(printf '%s' \"\${ROW}\" | cut -f3)
  echo \"[task \${SLURM_PROCID}] gpu=\${CUDA_VISIBLE_DEVICES} \$(basename \"\${SETTINGS}\" .json)\"
  echo \"    advanced=\$(basename \"\${ADVANCED}\")  filters=\$(basename \"\${FILTERS}\")\"
  source /miniforge3/etc/profile.d/conda.sh
  conda activate BindCraft
  cd '${BINDCRAFT_DIR}'
  python -u bindcraft.py --settings \"\${SETTINGS}\" --filters \"\${FILTERS}\" --advanced \"\${ADVANCED}\"
" &
pid=$!
rc=0; wait "$pid" || rc=$?
[[ "${rc}" -eq 124 || "${rc}" -eq 137 || "${rc}" -eq 143 ]] && { resubmit "walltime rc=${rc}"; exit 0; }
echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID} (${ARM}) fertig (rc=${rc})"
