# source this file before running anything:  export CS_DATA=/path/to/data; source runnable/env.sh
export CS_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export CS_DATA="${CS_DATA:?set CS_DATA to the folder that holds the data (see README)}"
export BENDER_REANALYSIS_DIR="$CS_REPO/bender_reanalysis"
export PYTHONPATH="$CS_REPO/code:$CS_REPO/bender_reanalysis:$PYTHONPATH"
