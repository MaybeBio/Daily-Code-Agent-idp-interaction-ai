#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
EXPERIMENT_DIR="src/experiments/intermapcnn"

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
python -m src.scripts.test_intermap \
    --config $EXPERIMENT_DIR/val_config.yml \
    --data_path data/validation/llps_vs_non_llps_val.csv \
    --score_path $EXPERIMENT_DIR/output/intermapcnn_baseline_h4_l1_val_score.csv \
    > $EXPERIMENT_DIR/logging/intermapcnn_baseline_h4_l1_val.log
python -m src.scripts.test_intermap \
    --config $EXPERIMENT_DIR/val_config.yml \
    --data_path data/validation/driver_vs_non_llps_val.csv \
    --score_path $EXPERIMENT_DIR/output/intermapcnn_baseline_h4_l1_val_driver_score.csv \
    > $EXPERIMENT_DIR/logging/intermapcnn_baseline_h4_l1_val_driver.log
python -m src.scripts.test_intermap \
    --config $EXPERIMENT_DIR/val_config.yml \
    --data_path data/validation/partner_vs_non_llps_val.csv \
    --score_path $EXPERIMENT_DIR/output/intermapcnn_baseline_h4_l1_val_partner_score.csv \
    > $EXPERIMENT_DIR/logging/intermapcnn_baseline_h4_l1_val_partner.log
