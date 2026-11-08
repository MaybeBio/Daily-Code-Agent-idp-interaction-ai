#!/bin/bash
# recompute_v3.sh — all labels with the code-review fixes of 2026-10-03 (unbiased tie-breaking at the hub cutoff,
# Jaccard null from the true hub-set sizes, split-half of the position-adjusted z, one frame rule, truncated-DCD check).
# Writes new folders only (old outputs are kept). Run from redo_evidence/.
# Steps 1-3 read the NAS. When it stalls, the label scripts exit with code 75 and are rerun here; they keep what is
# already labelled, so a rerun continues where the stall happened. To start a label folder from scratch after a code
# change, delete it first.
set -e
cd "$(dirname "$0")"
PY=~/.venvs/bender/bin/python
resumable() {
  local k rc
  for k in $(seq 1 40); do
    rc=0
    "$@" || rc=$?
    [ $rc -eq 0 ] && return 0
    [ $rc -eq 75 ] || return $rc
    echo "[$(date +%H:%M)] NAS stall (attempt $k), continuing in 3 min"
    sleep 180
  done
  return 75
}
export WORKERS=8
echo "[$(date +%H:%M)] 1/6 CALVADOS-2, all 10,282 (hub_frames.py)"
resumable python3 hub_frames.py --ids ids_all_10282.csv --out hub_frames_v3_cal --resume
echo "[$(date +%H:%M)] 2/6 Mpipi-GG 908 (cluster code path)"
resumable python3 gaivi/hub_labels_gaivi.py --layout lammps --root "$BSS_DATA/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k " --fasta labels_input_mpipi.fasta --out hub_frames_v3_mpipi --workers 8
echo "[$(date +%H:%M)] 3/6 Mpipi-Recharged pilot 289 (cluster code path)"
resumable python3 gaivi/hub_labels_gaivi.py --layout openmm --root "$BSS_DATA/MPIPI Recharged/300K" --fasta labels_input_recharged.fasta --out hub_frames_v3_recharged --workers 8 --no-require-done
echo "[$(date +%H:%M)] 4/6 PED labels (18 development, 8 blind) and the five-pool null"
mkdir -p v2_backup && cp -n hub_predictor/experimental_occupancy.csv hub_predictor/fresh_ped_occupancy.csv fresh_ped_test_summary.txt hub_predictor/null_occupancy_multi.csv v2_backup/ 2>/dev/null || true
python3 exp_occupancy.py
$PY fresh_ped_test.py > fresh_ped_test_v3.log 2>&1
python3 null_multi.py 5
echo "[$(date +%H:%M)] 5/6 dataset"
mkdir -p hub_predictor/data_v3
HUB_LABELS=v3 python3 hub_predictor/hub_dataset.py hub_predictor/data_v3
cp hub_predictor/experimental_occupancy.csv hub_predictor/data_v3/
echo "[$(date +%H:%M)] 6/6 cross-validation on the new labels"
cd hub_predictor && $PY cv_full.py --data data_v3 > cv_full_v3_run.log 2>&1 && cd ..
echo "[$(date +%H:%M)] done"
