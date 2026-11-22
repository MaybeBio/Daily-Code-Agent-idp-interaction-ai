#!/usr/bin/bash
set -euo pipefail

# GROMACS executable
GROMACS=gmx

usage() {
    echo "Usage: $0 -pdb <file.pdb> -pkl <file.pkl> -out <output_directory>" >&2
    exit 1
}

PDB="" PKL="" OUT="bAIes_simulations"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -pdb|--pdb) PDB="${2:-}"; shift 2 ;;
        -pkl|--pkl) PKL="${2:-}"; shift 2 ;;
        -out|--out) OUT="${2:-}"; shift 2 ;;
        -h|--help)  usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

# Check that all three were given
[[ -z "$PDB" || -z "$PKL" || -z "$OUT" ]] && usage

# Check that the inputs exist
[[ -f "$PDB" ]] || { echo "PDB file not found: $PDB" >&2; exit 1; }
[[ -f "$PKL" ]] || { echo "PKL file not found: $PKL" >&2; exit 1; }

mkdir -p $OUT #"$(dirname "$OUT")"

echo "PDB: $PDB"
echo "PKL: $PKL"
echo "OUT: $OUT"


# PDB model reading
BAIES_PDB=${OUT}/`basename ${PDB}`
echo 6 | ${GROMACS} pdb2gmx -f ${PDB} -o ${BAIES_PDB} -ignh -water none -p ${BAIES_PDB/.pdb}.top -i ${BAIES_PDB/.pdb}.itp

# Preprocessing
python ./scripts/preprocess_bAIes.py -pdb ${BAIES_PDB} -pkl ${PKL} -out ${OUT}/bAIes_params.dat


