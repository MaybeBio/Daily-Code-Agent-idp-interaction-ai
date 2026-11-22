#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"

# get the project main directory
SCRIPT_DIR="$SCRIPT_PATH"

for i in $(seq 1 4);
do
    SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done

# change directory to project main directory
cd $SCRIPT_DIR

# run script
python -m src.scripts.train --config src/experiments/esm2/train_config.yml
