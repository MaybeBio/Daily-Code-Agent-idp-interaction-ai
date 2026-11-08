#!/bin/bash
# label_recharged_v3.sh — add newly copied Mpipi-Recharged runs to hub_frames_v3_recharged (same code path and settings
# as recompute_v3.sh step 3; the label script skips IDs already labelled and appends the rest). Resumable after a NAS
# stall (exit code 75). Input: labels_input_recharged_present.fasta, the IDs with a folder on the NAS (2026-10-04: 5,421; 2026-10-09: all 10,282).
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
echo "[$(date +%H:%M)] Mpipi-Recharged labels, incremental"
rc=0
resumable python3 gaivi/hub_labels_gaivi.py --layout openmm --root "/Volumes/User Homes/Nu Project/MPIPI Recharged/300K" \
  --fasta labels_input_recharged_present.fasta --out hub_frames_v3_recharged --workers 8 --no-require-done || rc=$?
echo "[$(date +%H:%M)] exit $rc"
exit $rc
