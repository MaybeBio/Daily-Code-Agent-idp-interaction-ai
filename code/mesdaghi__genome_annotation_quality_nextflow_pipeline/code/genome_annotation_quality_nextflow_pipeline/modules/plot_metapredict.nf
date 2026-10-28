#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_METAPREDICT {
    tag { dataset_name }

    label 'cpu'

    publishDir "results/metapredict_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(csv_file)

    output:
    tuple val(dataset_name), path("*.png")

    script:
    """ 
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/plot_metapredict.py ${csv_file}

    # Where is the output file going. Need to update this to be more flexible. 
    """
}