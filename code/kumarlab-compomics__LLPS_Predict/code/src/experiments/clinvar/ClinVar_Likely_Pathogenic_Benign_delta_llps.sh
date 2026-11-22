#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
EXPERIMENT_DIR="src/experiments/clinvar"

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
python -m src.scripts.predict_delta_llps \
    --config $EXPERIMENT_DIR/ClinVar_Likely_Pathogenic_Benign_delta_llps_config.yml \
    --data_path data/clinvar/ClinVar_mutations_full.csv \
    --score_path $EXPERIMENT_DIR/output/ClinVar_mutations_full_delta_llps.csv \
    > $EXPERIMENT_DIR/logging/ClinVar_mutations_full_delta_llps.log
