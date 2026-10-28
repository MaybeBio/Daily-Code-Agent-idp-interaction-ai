#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * BUILD_SUMMARY_TSV — writes one TSV per dataset with the headline
 * metrics shown on the pipeline's plots.
 *
 * Every analysis branch is optional. Nextflow's `path` input has no clean
 * "optional" mode in DSL2, so main.nf passes a NO_FILE_* placeholder from
 * assets/ for any branch that did not run, and the matching flag is simply
 * left off the command line. build_summary_tsv.py then leaves that column
 * blank so the TSV schema stays the same across runs.
 */


process BUILD_SUMMARY_TSV {
    tag { dataset_name }

    label 'cpu'

    publishDir "results", mode: 'copy'

    input:
    tuple val(dataset_name),
          path(fasta_file),
          path(plddt_pkl),
          path(metapredict_csv),
          path(psauron_csv),
          path(ipr_xml)

    output:
    tuple val(dataset_name), path("${dataset_name}_summary.tsv")

    script:
    def provided = { f -> !f.name.startsWith('NO_FILE') }
    def args = []
    if ( provided(plddt_pkl) ) {
        args << "--plddt-pkl ${plddt_pkl}"
        args << "--plddt-reference-csv ${projectDir}/reference/plddt_model_organisms.csv"
    }
    if ( provided(metapredict_csv) ) {
        args << "--metapredict-csv ${metapredict_csv}"
    }
    if ( provided(psauron_csv) ) {
        args << "--psauron-csv ${psauron_csv}"
    }
    if ( provided(ipr_xml) ) {
        args << "--ipr-fasta ${fasta_file}"
        args << "--ipr-xml ${ipr_xml}"
        args << "--plot-interpro-dir ${projectDir}/bin"
    }
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/build_summary_tsv.py \\
        --name "${dataset_name}" \\
        --fasta ${fasta_file} \\
        ${args.join(' \\\n        ')} \\
        --out ${dataset_name}_summary.tsv
    """
}
