#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
EXPERIMENT_DIR="src/experiments/esm2"

# get the project main directory
SCRIPT_DIR="$SCRIPT_PATH"

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
python -m src.scripts.test \
    --config $EXPERIMENT_DIR/test_config.yml \
    --data_path data/test/llps_vs_non_llps_disjoint_test.csv \
    --score_path $EXPERIMENT_DIR/output/esm2_650m_baseline_h4_l1_test_score.csv \
    > $EXPERIMENT_DIR/logging/esm2_650m_baseline_h4_l1_test.log
