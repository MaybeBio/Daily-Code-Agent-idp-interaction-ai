#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process INTERPROSCAN {
    tag { "${dataset_name}_${chunk_file.simpleName}" }

    label 'cpu'

    input:
    tuple val(dataset_name), path(chunk_file)

    output:
    tuple val(dataset_name), path("*.tsv"), path("*.xml")

    script:
    """
    # module purge
    # module load adoptopenjdk/11.0.12+7
    # module load interproscan/5.66-98.0

    export TMPDIR=\$(mktemp -d)


    /opt/interproscan/interproscan.sh \
        -i ${chunk_file} \
        -f TSV,XML \
        --cpu ${task.cpus} \
        --tempdir \$TMPDIR
    """
}
