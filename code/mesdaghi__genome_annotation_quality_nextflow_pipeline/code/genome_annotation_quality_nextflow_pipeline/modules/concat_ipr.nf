process CONCAT_IPR {
    tag { dataset_name }

    label 'cpu'

    publishDir "results/interpro", mode: 'copy'

    input:
    tuple val(dataset_name), path(tsv_files), path(xml_files)

    output:
    tuple val(dataset_name), path("${dataset_name}_interproscan.tsv"), path("${dataset_name}_interproscan.xml")

    script:
    // Outputs are named after the dataset so that several inputs don't
    // overwrite each other in results/interpro.
    """
    python3 ${projectDir}/bin/concat_xml.py ${xml_files} -o ${dataset_name}_interproscan.xml
    python3 ${projectDir}/bin/concat_tsv.py ${tsv_files}
    mv concatenated.tsv ${dataset_name}_interproscan.tsv
    """
}