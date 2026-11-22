#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
EXPERIMENT_DIR="src/experiments/mutagenesis"

# get the project main directory
SCRIPT_DIR="$SCRIPT_PATH"

for i in $(seq 1 4);
do
    SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done

# change directory to project main directory
cd $SCRIPT_DIR

# create output and logging directories if not exist
OUTPUT_DIR="$EXPERIMENT_DIR/output/denovo_db_ssc_top_genes"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# run script
python -m src.scripts.predict_mutagenesis \
    --config $EXPERIMENT_DIR/mutagenesis_denovo_db_ssc_top_genes.yml \
    --data_path data/denovo_db/denovo-db_1.8_ssc_Top_Genes_All.csv \
    --score_path $EXPERIMENT_DIR/output/denovo_db_ssc_top_genes/{i}_Top_Genes_Unique_{gene_symbol}_mutagenesis.csv
