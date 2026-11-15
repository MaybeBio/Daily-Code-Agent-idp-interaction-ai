#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * UniPPI-OWL v1.8 — Nextflow DSL2 High-Throughput Structure-to-Affinity Pipeline
 * Integrates unippi_cli.py as an automated batch node for AlphaFold3 / PDB screening.
 */

params.input_dir      = "${projectDir}/demo_structures"
params.output_dir     = "${projectDir}/demo_results"
params.temp           = 298.15
params.lambda_pae     = 8.0
params.ionic_strength = 150.0
params.ph             = 7.2
params.ptm_state      = "native"
params.mutate         = "A:LEU5PRO,A:PRO5LEU,A:ASP2ALA,A:TYR5ALA"

process UNIPPI_BATCH_SCREEN {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path input_folder

    output:
    path "unippi_out/master_interaction_ranking.csv", emit: master_csv
    path "unippi_out/unippi_full_report.json",        emit: full_json
    path "unippi_out/converted_structures/*",         emit: sanitized_structures

    script:
    def mut_arg = params.mutate ? "--mutate ${params.mutate}" : ""
    """
    python3 ${projectDir}/unippi_cli.py \
        -i ${input_folder} \
        -o unippi_out \
        --temp ${params.temp} \
        --lambda-pae ${params.lambda_pae} \
        --ionic-strength ${params.ionic_strength} \
        --ph ${params.ph} \
        --ptm-state ${params.ptm_state} \
        ${mut_arg}
    """
}

workflow {
    in_ch = Channel.fromPath(params.input_dir, checkIfExists: true)
    UNIPPI_BATCH_SCREEN(in_ch)
}
