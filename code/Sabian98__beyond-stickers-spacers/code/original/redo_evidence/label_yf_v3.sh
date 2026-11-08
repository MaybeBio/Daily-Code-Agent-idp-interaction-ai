#!/bin/bash
# label_shuffles_v3.sh — v3 hub labels for the Y-to-F mutants of the earlier BENDER-DYN version
# (800 proteins with every tyrosine replaced by phenylalanine, CALVADOS-2, OpenMM DCD + top.pdb, same protocol as BENDER). Same code path and settings
# as the Mpipi-Recharged labels; incremental and resumable.
cd "$(dirname "$0")"
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
echo "[$(date +%H:%M)] Y-to-F labels, incremental"
rc=0
resumable python3 gaivi/hub_labels_gaivi.py --layout openmm --root "/Volumes/User Homes/Nu Project/BENDER-DYN/Mutation Experiment/Trajectories" \
  --fasta labels_input_yf.fasta --out hub_frames_v3_yf --workers ${WORKERS:-3} --no-require-done || rc=$?
echo "[$(date +%H:%M)] exit $rc"
exit $rc
