#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process GENERATE_REPORT {
    tag { dataset_name }

    label 'cpu'
    
    publishDir "results", mode: 'copy'

    input:
    tuple val(dataset_name), path(all_pngs)

    output:
    path "${dataset_name}_report.pdf"

    script:
    """

    python3 ${projectDir}/bin/generate_report.py \
        ${dataset_name} \
        ${dataset_name}_report.pdf \
        ${all_pngs}
    """
}