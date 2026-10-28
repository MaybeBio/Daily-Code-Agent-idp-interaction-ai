#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SPLIT_JSON {
    tag { dataset_name }

    label 'cpu'

    input:
    tuple path(json_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("json_chunks/*.json")

    script:
    """
    mkdir -p json_chunks
    python3 ${projectDir}/bin/split_json.py ${json_file} ${params.chunk_size} json_chunks
    """
}