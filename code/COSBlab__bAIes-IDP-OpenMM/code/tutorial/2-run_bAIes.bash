#!/usr/bin/bash
set -euo pipefail

# Target number of frames:
Nsteps=500 #00000
# Simulation timestep in fs
tstep=0.002
# Save interval in number of steps.
save_interval=50 #000

# GROMACS executable
GROMACS=gmx

usage() {
    echo "Usage: $0 -pdb <file.pdb> -baies <file.pkl> -out <output_directory>" >&2
    exit 1
}

PDB="" BAIES="" OUT="bAIes_simulations"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -pdb|--pdb) PDB="${2:-}"; shift 2 ;;
        -baies|--baies) BAIES="${2:-}"; shift 2 ;;
        -out|--out) OUT="${2:-}"; shift 2 ;;
        -h|--help)  usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

# Check that all arguments were given
[[ -z "$PDB" || -z "$OUT" ]] && usage

# Check that the inputs exist
[[ -f "$PDB" ]] || { echo "PDB file not found: $PDB" >&2; exit 1; }

if [[ -n "$BAIES" ]]; then
    [[ -f "$BAIES" ]] || { echo "BAIES file not found: $BAIES" >&2; exit 1; }
fi

mkdir -p $OUT #"$(dirname "$OUT")"

echo "PDB: $PDB"
echo "BAIES: $BAIES"
echo "OUT: $OUT"


if [[ -n "$BAIES" ]]; then
echo "Running bAIes..."
python ./scripts/baies_openMM.py \
                  -pdb ${PDB} \
                  -xtc ${OUT}/bAIes_trajectory.xtc \
                  -ff ./force_field/coil_202602.xml \
                  -cmap ./force_field/cmap_20240524.cmap \
                  -temp 298.15 \
                  -friction 0.01 \
                  -tstep 0.002 \
                  -Nsteps ${Nsteps} \
                  -save_xtc ${save_interval} \
                  -save_stdout ${save_interval} \
                  -baies ${BAIES}
                  
else
echo "Running coil..."
python ./scripts/baies_openMM.py \
                  -pdb ${PDB} \
                  -xtc ${OUT}/bAIes_trajectory.xtc \
                  -ff ./force_field/coil_202602.xml \
                  -cmap ./force_field/cmap_20240524.cmap \
                  -temp 298.15 \
                  -friction 0.01 \
                  -tstep 0.002 \
                  -Nsteps ${Nsteps} \
                  -save_xtc ${save_interval} \
                  -save_stdout ${save_interval}
fi





