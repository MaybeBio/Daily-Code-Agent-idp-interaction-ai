#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"

# script's parent directory
SCRIPT_DIR="$SCRIPT_PATH"

# get the project main directory
for i in $(seq 1 4);
do
    SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done

# change directory to project main directory
cd $SCRIPT_DIR

# run script
python -m src.scripts.train_ensemble --config src/experiments/ensemble_partner/train_config.yml
