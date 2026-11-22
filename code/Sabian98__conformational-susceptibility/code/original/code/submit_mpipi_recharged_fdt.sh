#!/bin/bash
# =============================================================================
# submit_mpipi_recharged_fdt.sh
#
# FDT susceptibility of the Mpipi-Recharged (OpenMpipi) runs at one temperature on GAIVI
# (mpipi_recharged_fdt.py run): one Slurm job, CPU platform, 16 workers, 32 GB, 24 h. Resumable: a resubmission
# skips proteins already in the output table.
#
# Usage (on the HPC, in the environment that runs run_openmpipi_idp.py; start it with bash: it submits the job
# itself and never changes directory):
#   bash /path/to/submit_mpipi_recharged_fdt.sh RUNS TEMP [IDS]
#     RUNS  run_openmpipi_idp.py --outdir root, holding <TEMP>K/<ID>/<ID>.dcd + top.pdb (e.g. /general/taseefr)
#     TEMP  temperature in K
#     IDS   optional file with one ID per line
# First try one protein interactively:
#   python /path/to/mpipi_recharged_fdt.py check --runs RUNS --temp TEMP
# Needs mpipi_recharged_fdt.py, TODO_T2_fdt_computation.py, TODO_T2_trajectory_loader.py, bender_common.py and
# graph_features.py in this script's folder.
# Output: $OUT_ROOT/recharged_fdt_<TEMP>K.csv (+ .provenance.json) and $OUT_ROOT/npz_<TEMP>K/; copy them to
# BENDER-MT/Traj Level Analysis/Analysis_FDT_RECHARGED/. Nothing is written into the script folder.
# Slurm logs go to $OUT_ROOT/logs; if that cannot be created, everything goes under ~/mpipi_recharged_fdt.
#   OUT_ROOT  default /data/taseefr/mpipi_recharged_fdt
# =============================================================================
set -euo pipefail
RUNS="${1:?usage: bash submit_mpipi_recharged_fdt.sh RUNS TEMP [IDS]}"
TEMP="${2:?usage: bash submit_mpipi_recharged_fdt.sh RUNS TEMP [IDS]}"
IDS="${3:-}"
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
RUNS="$(realpath "$RUNS")"
OUT_ROOT="${OUT_ROOT:-/data/taseefr/mpipi_recharged_fdt}"
if ! { mkdir -p "$OUT_ROOT/logs" 2>/dev/null && [ -w "$OUT_ROOT/logs" ]; }; then
    OUT_ROOT="$HOME/mpipi_recharged_fdt"
    mkdir -p "$OUT_ROOT/logs"
fi
LOG_DIR="$OUT_ROOT/logs"
for f in mpipi_recharged_fdt.py TODO_T2_fdt_computation.py TODO_T2_trajectory_loader.py bender_common.py graph_features.py; do
    [ -f "$SCRIPT_DIR/$f" ] || { echo "missing $SCRIPT_DIR/$f" >&2; exit 1; }
done
IDS_OPT=""
[ -n "$IDS" ] && IDS_OPT="--ids \"$(realpath "$IDS")\""
sbatch --job-name="recharged_fdt_${TEMP}K" --ntasks=1 --cpus-per-task=16 --mem=32G --time=24:00:00 \
    --output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" \
    --wrap="BENDER_REANALYSIS_DIR=\"$SCRIPT_DIR\" python \"$SCRIPT_DIR/mpipi_recharged_fdt.py\" run --runs \"$RUNS\" --temp $TEMP --out \"$OUT_ROOT\" --workers 16 --platform CPU $IDS_OPT"
echo "submitted: output $OUT_ROOT/recharged_fdt_${TEMP}K.csv, logs $LOG_DIR"
