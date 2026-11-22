#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
EXPERIMENT_DIR="src/experiments/critical_region"

# script's parent directory
SCRIPT_DIR="$SCRIPT_PATH"

# get the project main directory
for i in $(seq 1 4);
do
    SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done

# change directory to project main directory
cd $SCRIPT_DIR

# create logging directory
LOGGING_DIR="${EXPERIMENT_DIR}/logging"
if [ ! -d $LOGGING_DIR ]; then
    mkdir $LOGGING_DIR
fi

# create score directory
SCORE_DIR="${EXPERIMENT_DIR}/output"
if [ ! -d $SCORE_DIR ]; then
    mkdir $SCORE_DIR
fi

# run script
python -m src.scripts.predict_critical_region \
    --config $EXPERIMENT_DIR/phasepro_critical_region_config.yml \
    --data_path data/llps_driving_region/PhaSePro_LLPS.csv \
    --score_path $EXPERIMENT_DIR/output/PhaSePro_Critical_Region.csv
