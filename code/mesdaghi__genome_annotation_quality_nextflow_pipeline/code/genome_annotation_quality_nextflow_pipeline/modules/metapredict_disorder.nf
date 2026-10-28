#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process METAPREDICT_DISORDER {
    tag { dataset_name }

    label 'cpu'

    publishDir "results/metapredict", mode: 'copy'

    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("${dataset_name}_metapredict.csv")

    script: 
    // # This will need to change a lot as it's no longer sourcing miniconda
    """
    echo "Running on: \$(hostname)"
    echo "Start time: \$(date)"

    # module load cuda/13.0.2
    # source ~/miniconda3/etc/profile.d/conda.sh
    # conda activate metapredict


    metapredict-predict-disorder \
        -o ${dataset_name}_metapredict.csv \
        ${fasta_file}

    echo "End time: \$(date)"
    """
}