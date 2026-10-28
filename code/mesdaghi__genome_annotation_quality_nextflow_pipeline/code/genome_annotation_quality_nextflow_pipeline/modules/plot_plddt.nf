#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_PLDDT {
    tag { dataset_name }

    label 'cpu'

    publishDir "results/plots", mode: 'copy' // having dirs of the same name sometimes causes issues with singularity.

    input:
    tuple val(dataset_name), path(pkl_file), path(metapredict_csv)

    output:
    tuple val(dataset_name), path("plddt_*_${dataset_name}.png")

    script:
    // metapredict_csv is a NO_FILE_* placeholder when the metapredict
    // branch is disabled; plot_plddt.py then skips the disorder scatter.
    def disorder_arg = metapredict_csv.name.startsWith('NO_FILE') ? '' : metapredict_csv
    """
    export MPLCONFIGDIR=\$(mktemp -d)
    
    python3 ${projectDir}/bin/plot_plddt.py ${dataset_name} ${disorder_arg}
    """
}
