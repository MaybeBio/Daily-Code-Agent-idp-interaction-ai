# source this file before running anything:  export BSS_DATA=/path/to/data; source runnable/env.sh
export BSS_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export BSS_DATA="${BSS_DATA:?set BSS_DATA to the folder that holds the unpacked data (see README)}"
export BENDER_REANALYSIS_DIR="$BSS_REPO/bender_reanalysis"
export PYTHONPATH="$BSS_REPO/redo_evidence:$BSS_REPO/bender_reanalysis:$BSS_REPO/methods_paper/figures:$PYTHONPATH"
