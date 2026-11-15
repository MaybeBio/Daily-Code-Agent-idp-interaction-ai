#!/usr/bin/env bash
#SBATCH --job-name=boltz_cmp
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --account=ab011
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --mem=100G
#SBATCH --output=/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/boltz2_env/actif_vs_idp/boltz_%A.log

set -uo pipefail
REPO=/capstor/store/cscs/swissai/ab011/jkuhn/.uenv-images
VENV=/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/boltz2_env/venv
B2=/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/boltz2_env/actif_vs_idp
CACHE=/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/boltz2_env/boltz_cache
srun --ntasks=1 --gpus-per-task=1 --cpus-per-task=16 uenv --repo "${REPO}" run --view=default pytorch/v2.9.1:v2 -- bash -c "
  source ${VENV}/bin/activate
  boltz predict ${B2}/inputs --out_dir ${B2}/out --cache ${CACHE} --output_format pdb --no_kernels 2>&1 | tail -100
"
echo "[$(date -Is)] done"
