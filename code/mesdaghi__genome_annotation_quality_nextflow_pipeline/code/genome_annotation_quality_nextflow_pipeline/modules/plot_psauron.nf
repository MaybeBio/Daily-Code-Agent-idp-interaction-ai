#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_PSAURON {
    tag { dataset_name }
    
    label 'cpu'

    publishDir "results/psauron_plots", mode: 'copy' // having dirs of the same name sometimes causes issues with singularity.

    input:
    tuple val(dataset_name), path(psauron_csv)

    output:
    tuple val(dataset_name), path("psauron_categorical_distribution_${dataset_name}.png")

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    # This is an issue as the HPC doesn't like linking to files on the home directory. Need to update this to be more flexible.
    # The output is named after the dataset so that several inputs don't
    # overwrite each other in results/psauron_plots.
    python3 ${projectDir}/bin/plot_psauron_distribution.py \
        ${projectDir}/reference/combined_psauron_results.csv \
        ${psauron_csv} \
        --out psauron_categorical_distribution_${dataset_name}.png
    """
}