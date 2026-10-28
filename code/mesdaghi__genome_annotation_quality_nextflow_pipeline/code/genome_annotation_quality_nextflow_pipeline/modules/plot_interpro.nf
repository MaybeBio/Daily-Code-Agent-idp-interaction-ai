#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PLOT_INTERPRO {
    tag { dataset_name }

    label 'cpu'

    publishDir "results/interpro_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(fasta_file), path(xml_file)

    output:
    tuple val(dataset_name),
          path("interpro_domain_coverage_${dataset_name}.png"),
          path("interpro_merged_coverage_distribution_${dataset_name}.png"),
          path("interpro_summary_${dataset_name}.tsv")

    script:
    """
    python3 ${projectDir}/bin/plot_interpro.py \
        ${fasta_file} \
        ${xml_file} \
        --json ${params.ipr_reference_json} \
        --name "${dataset_name}" \
        --outdir .

    # Rename to stable filenames so the report step can find them. The
    # dataset name keeps several inputs apart in results/interpro_plots.
    mv query_ipr_coverage.png                       interpro_domain_coverage_${dataset_name}.png
    mv query_ipr_merged_coverage_distribution.png   interpro_merged_coverage_distribution_${dataset_name}.png
    mv query_ipr_summary.tsv                        interpro_summary_${dataset_name}.tsv
    """
}