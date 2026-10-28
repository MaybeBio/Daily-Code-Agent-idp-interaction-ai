#!/usr/bin/env bash
# The v1.2 re-sample analysis (docs/RESAMPLE_DESIGN.md), in order. Every step can be rerun; 04 resumes.
#   BIOEMU_PY        python of the conversion environment (requirements.txt)
#   BSC_ROOT         where `bsc fetch` put data/;  ARCHIVE  Predictions.tar.gz, for ensembles not yet extracted
#   RESAMPLE_SAMPLES the notebook's output folder;  RESAMPLE_WORK  a folder with ~2 GB free
#   SIM_WORK         the simulation analysis's work folder (PULCHRA and Pepsi-SAXS under tools/)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-uv run python}
$PY 01_inventory.py ${ARCHIVE:+--archive "$ARCHIVE"}
"$BIOEMU_PY" 02_convert.py
$PY 03_rg.py
$PY 04_saxs.py
$PY 05_chi2.py
$PY 06_analyse.py
