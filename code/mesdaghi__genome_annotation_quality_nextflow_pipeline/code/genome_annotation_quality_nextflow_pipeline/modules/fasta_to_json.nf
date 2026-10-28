#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process FASTA_TO_JSON {
    
    tag { dataset_name }

    label 'cpu'

    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple path("out.json"), val(dataset_name)

    script:
    """
    python3 ${projectDir}/bin/fasta_to_protenix_json_all.py ${fasta_file} out.json
    """
}