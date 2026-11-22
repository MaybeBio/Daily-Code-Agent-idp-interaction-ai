#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
EXPERIMENT_DIR="src/experiments/denovo_db"

# get the project main directory
SCRIPT_DIR="$SCRIPT_PATH"

for i in $(seq 1 4);
do
    SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done

# change directory to project main directory
cd $SCRIPT_DIR

# create output and logging directories if not exist
OUTPUT_DIR="$EXPERIMENT_DIR/output"
LOGGING_DIR="$EXPERIMENT_DIR/logging"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

if [ ! -d "$LOGGING_DIR" ]; then
    mkdir -p "$LOGGING_DIR"
fi

# run script
python -m src.scripts.test_ensemble \
    --config $EXPERIMENT_DIR/denovo-db_1.8_ssc2_score_config.yml \
    --data_path data/denovo_db/denovo_db_ssc_mutations_full.csv \
    --score_path $EXPERIMENT_DIR/output/denovo_db_ssc_mutations_full_score.csv \
    > $EXPERIMENT_DIR/logging/denovo_db_ssc_mutations_full_score.log
