#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process COLLECT_CHUNKS {
    tag { dataset_name }

    label 'cpu'
   
    input:
    tuple val(dataset_name), path(pred_dirs, stageAs: "chunk*/*")

    output:
    tuple val(dataset_name), path("${dataset_name}_all_predictions")

    script:
    """
    mkdir -p ${dataset_name}_all_predictions
    printf '%s\\n' ${pred_dirs} | xargs -P 8 -I{} cp -r {}/. ${dataset_name}_all_predictions/
    """

}
