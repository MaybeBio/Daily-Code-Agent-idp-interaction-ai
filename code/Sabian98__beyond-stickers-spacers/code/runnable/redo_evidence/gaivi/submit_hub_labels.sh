#!/bin/bash
# =============================================================================
# submit_hub_labels.sh — per-frame hub labels on GAIVI (one CPU job; standard-library Python only)
#
# Run from the folder that holds hub_labels_gaivi.py, hub_frames.py, hub_z_v2_stdlib.py and xtc_stdlib.py:
#   sbatch submit_hub_labels.sh <layout> <root> <fasta> <outdir>
#
#   Mpipi-Recharged, full BENDER (OpenMPIPI/submit_openmpipi.sh output):
#     sbatch submit_hub_labels.sh openmm /general/taseefr/300K ~/OpenMPIPI/bender_10282.fasta labels_recharged_all
#   Design batch:
#     sbatch submit_hub_labels.sh openmm /data/taseefr/hub_design/300K    hub_design.fasta labels_design_calvados
#     sbatch submit_hub_labels.sh openmm /general/taseefr/hub_design/300K hub_design.fasta labels_design_recharged
#     sbatch submit_hub_labels.sh lammps /data/taseefr/mpipi_hubdesign    hub_design.fasta labels_design_mpipi
#
# Rerun the same command whenever more simulations have finished: IDs already labelled are skipped and new ones
# appended; unfinished runs are listed in <outdir>/failed.csv and picked up next time.
# =============================================================================
#SBATCH --job-name=hub_labels
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=hub_labels_%j.out

[ $# -eq 4 ] || { echo "Usage: sbatch $0 <openmm|lammps> <root> <fasta> <outdir>"; exit 1; }
cd "${SLURM_SUBMIT_DIR:-.}"

# any Python >= 3.6 works; use the CALVADOS env when it is there
source /apps/anaconda3/etc/profile.d/conda.sh 2>/dev/null && conda activate calvados 2>/dev/null || true

python3 hub_labels_gaivi.py --layout "$1" --root "$2" --fasta "$3" --out "$4" --workers "${SLURM_CPUS_PER_TASK:-8}"
