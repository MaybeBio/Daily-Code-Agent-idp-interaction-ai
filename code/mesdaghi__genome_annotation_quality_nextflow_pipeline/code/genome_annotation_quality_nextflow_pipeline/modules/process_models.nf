#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PROCESS_MODELS {
    tag { dataset_name }

    label 'cpu'
    
    input:
    tuple val(dataset_name), path(pred_dir)

    output:
    tuple val(dataset_name), path("plddt_all_values_${dataset_name}_all_one.pkl")

    script:
    """
    python3 ${projectDir}/bin/process_models.py \
      ${pred_dir} \
      plddt_all_values_${dataset_name}_all_one.pkl
    """
}