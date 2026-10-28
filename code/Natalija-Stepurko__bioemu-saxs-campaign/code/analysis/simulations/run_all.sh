#!/usr/bin/env bash
# Rerun the simulation analysis (docs/SIMULATION_DESIGN.md).
#   SIM_WORK   downloads, tools, trajectories and computed curves (default analysis/simulations/work, not tracked)
#   SIM_PY     a Python with requirements.txt (the simulation environment)
#   BSC_ROOT   where `bsc fetch` put data/PeptoneDB-SAXS (the measured curves)
#   WORKERS    parallel single-thread processes (default 6)
# Every step resumes: finished proteins are skipped, interrupted runs continue from their last checkpoint.
set -euo pipefail
cd "$(dirname "$0")"
PY=${SIM_PY:-python3}
W=${WORKERS:-6}
uv run python 01_prepare_inputs.py
$PY 02_fetch_sources.py
$PY 03_circularity.py
# pilot (design §8): two proteins in full, the first 20 frames of the longest; then the plan and the rest
$PY 05_run.py --labels SASDJH6 SASDJT7 --workers 2 &
$PY 05_run.py --labels SASDLK9 --frames 20 --workers 1
wait
$PY 04_plan.py
$PY 05_run.py --workers "$W"
$PY 06_saxs.py --workers "$W"
uv run python 07_analyse.py
