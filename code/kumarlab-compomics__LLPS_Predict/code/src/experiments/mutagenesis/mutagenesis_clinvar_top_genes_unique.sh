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
OUTPUT_DIR="$EXPERIMENT_DIR/output/clinvar_top_genes_unique"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# run script
python -m src.scripts.predict_mutagenesis \
    --config $EXPERIMENT_DIR/mutagenesis_clinvar_top_genes_unique.yml \
    --data_path data/clinvar/ClinVar_Likely_Pathogenic_Benign_Top_Genes_Unique.csv \
    --score_path $EXPERIMENT_DIR/output/clinvar_top_genes_unique/{i}_Top_Genes_Unique_{gene_symbol}_mutagenesis.csv
