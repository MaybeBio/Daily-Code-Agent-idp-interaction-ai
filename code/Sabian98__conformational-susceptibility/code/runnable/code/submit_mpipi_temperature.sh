#!/bin/bash
# =============================================================================
# submit_mpipi_temperature.sh
#
# Analyse the MPIPI-GG runs of the 992 paired sequences at one temperature on GAIVI
# (mpipi_temperature_hpc.py): one Slurm job, 16 CPUs, 32 GB, 12 h.
#
# Usage (on the HPC; start it with bash: it submits the job itself and never changes directory):
#   bash /path/to/submit_mpipi_temperature.sh RUNS TEMP
#     RUNS  folder holding one sub-folder per sequence (seq.dat, traj.xtc, run.in)
#     TEMP  temperature in K (checked against each run.in)
# Needs mpipi_temperature_hpc.py, bender_common.py and graph_features.py in this script's folder
# (the last two from Scripts/bender_reanalysis) and a Python environment with numpy, scipy, pandas,
# mdtraj and networkx.
# Output: $OUT_ROOT/mpipi_<TEMP>K.csv (+ .provenance.json); copy both back to
# BENDER-BIO/MPIPI/MPIPI Results/. Nothing is written into the script folder, which may be read-only.
# Slurm logs go to $OUT_ROOT/logs; if that cannot be created, output and logs go under ~/mpipi_temperature.
#   OUT_ROOT  default /data/taseefr/mpipi_temperature
# =============================================================================
set -euo pipefail
RUNS="${1:?usage: bash submit_mpipi_temperature.sh RUNS TEMP}"
TEMP="${2:?usage: bash submit_mpipi_temperature.sh RUNS TEMP}"
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
RUNS="$(realpath "$RUNS")"
OUT_ROOT="${OUT_ROOT:-/data/taseefr/mpipi_temperature}"
if ! { mkdir -p "$OUT_ROOT/logs" 2>/dev/null && [ -w "$OUT_ROOT/logs" ]; }; then
    OUT_ROOT="$HOME/mpipi_temperature"
    mkdir -p "$OUT_ROOT/logs"
fi
LOG_DIR="$OUT_ROOT/logs"
for f in mpipi_temperature_hpc.py bender_common.py graph_features.py; do
    [ -f "$SCRIPT_DIR/$f" ] || { echo "missing $SCRIPT_DIR/$f" >&2; exit 1; }
done
sbatch --job-name="mpipi_${TEMP}K" --ntasks=1 --cpus-per-task=16 --mem=32G --time=12:00:00 \
    --output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" \
    --wrap="python \"$SCRIPT_DIR/mpipi_temperature_hpc.py\" --runs \"$RUNS\" --temp $TEMP --out \"$OUT_ROOT/mpipi_${TEMP}K.csv\" --workers 16"
echo "submitted: output $OUT_ROOT/mpipi_${TEMP}K.csv, logs $LOG_DIR"
