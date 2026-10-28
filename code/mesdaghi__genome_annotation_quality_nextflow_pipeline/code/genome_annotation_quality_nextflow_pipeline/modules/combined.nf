#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Processes for --report combined.
 *
 * Each plot process receives every input at once, as parallel lists that
 * main.nf has sorted by dataset name, and draws all of them on the same
 * figure alongside the reference organisms. Files are staged into numbered
 * sub-folders (stageAs) so that inputs which share a file name, such as the
 * merged InterProScan XML, never collide.
 */


// A single staged file arrives as a Path rather than a one-element list
def asList(x) {
    return x instanceof List ? x : [x]
}

// File name without any stageAs sub-folder (a staged path's .name is
// e.g. 'metapredict1/NO_FILE_METAPREDICT', not just the file name)
def baseFileName(f) {
    return f.toString().tokenize('/').last()
}

// Quote a dataset name for the shell
def q(String s) {
    return "'" + s.replace("'", "'\\''") + "'"
}


process PLOT_PLDDT_COMBINED {
    tag { prefix }

    label 'cpu'

    publishDir "results/combined_plots", mode: 'copy'

    input:
    tuple val(names), path(pkls, stageAs: 'pkl*/*'), path(metapredict_csvs, stageAs: 'metapredict*/*')
    val prefix

    output:
    path "${prefix}_plddt_*.png"

    script:
    // metapredict_csvs holds NO_FILE_METAPREDICT placeholders when the
    // metapredict branch is off; the disorder scatter is then skipped.
    def query_args = [asList(names), asList(pkls), asList(metapredict_csvs)]
        .transpose()
        .collect { name, pkl, csv ->
            "--query ${q(name)} ${pkl}" + (baseFileName(csv).startsWith('NO_FILE') ? '' : " ${csv}")
        }
        .join(' \\\n        ')
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/plot_plddt.py --combined ${prefix} \\
        ${query_args}
    """
}


process PLOT_METAPREDICT_COMBINED {
    tag { prefix }

    label 'cpu'

    publishDir "results/combined_plots", mode: 'copy'

    input:
    tuple val(names), path(csvs, stageAs: 'metapredict*/*')
    val prefix

    output:
    path "${prefix}_metapredict_mean_disorder_density_scipy.png"

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/plot_metapredict.py \\
        ${asList(csvs).join(' ')} \\
        --names ${asList(names).collect { q(it) }.join(' ')} \\
        --out-prefix ${prefix}
    """
}


process PLOT_PSAURON_COMBINED {
    tag { prefix }

    label 'cpu'

    publishDir "results/combined_plots", mode: 'copy'

    input:
    tuple val(names), path(csvs, stageAs: 'psauron*/*')
    val prefix

    output:
    path "${prefix}_psauron_categorical_distribution.png"

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/plot_psauron_distribution.py \\
        ${projectDir}/reference/combined_psauron_results.csv \\
        ${asList(csvs).join(' ')} \\
        --names ${asList(names).collect { q(it) }.join(' ')} \\
        --out ${prefix}_psauron_categorical_distribution.png
    """
}


process PLOT_INTERPRO_COMBINED {
    tag { prefix }

    label 'cpu'

    publishDir "results/combined_plots", mode: 'copy'

    input:
    tuple val(names), path(fastas, stageAs: 'fasta*/*'), path(xmls, stageAs: 'xml*/*')
    val prefix

    output:
    tuple path("${prefix}_interpro_domain_coverage.png"),
          path("${prefix}_interpro_merged_coverage_distribution.png"),
          path("${prefix}_interpro_summary.tsv")

    script:
    def query_args = [asList(names), asList(fastas), asList(xmls)]
        .transpose()
        .collect { name, fasta, xml -> "--query ${fasta} ${xml} ${q(name)}" }
        .join(' \\\n        ')
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/plot_interpro.py \\
        ${query_args} \\
        --json ${params.ipr_reference_json} \\
        --prefix ${prefix} \\
        --outdir .

    # Stable names so the report step can label them
    mv ${prefix}_ipr_coverage.png                      ${prefix}_interpro_domain_coverage.png
    mv ${prefix}_ipr_merged_coverage_distribution.png  ${prefix}_interpro_merged_coverage_distribution.png
    mv ${prefix}_ipr_summary.tsv                       ${prefix}_interpro_summary.tsv
    """
}


process GENERATE_COMBINED_REPORT {
    tag { prefix }

    label 'cpu'

    publishDir "results", mode: 'copy'

    input:
    path all_pngs
    val prefix
    val subtitle

    output:
    path "${prefix}_report.pdf"

    script:
    """
    python3 ${projectDir}/bin/generate_report.py \\
        "Combined report" \\
        ${prefix}_report.pdf \\
        ${all_pngs} \\
        --subtitle ${q(subtitle)}
    """
}


process COMBINE_SUMMARY_TSV {
    tag { prefix }

    label 'cpu'

    publishDir "results", mode: 'copy'

    input:
    path tsvs, stageAs: 'summary*/*'
    val prefix

    output:
    path "${prefix}_summary.tsv"

    script:
    // Every per-input TSV has the same header; keep the first one only.
    """
    first=1
    for f in ${asList(tsvs).join(' ')}; do
        if [ "\$first" = 1 ]; then
            cat "\$f"
            first=0
        else
            tail -n +2 "\$f"
        fi
    done > ${prefix}_summary.tsv
    """
}
