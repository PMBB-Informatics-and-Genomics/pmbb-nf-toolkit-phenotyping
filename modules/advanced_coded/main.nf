process ADVANCED_CODED {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy',
        saveAs: { it.endsWith('.manifest.tsv') ? null : it }

    input:
    tuple path(config_json), path(long_tsvs)

    output:
    path "*.diag.tsv",     emit: result
    path "*.manifest.tsv", emit: manifest

    script:
    def name = config_json.baseName.replace('.json', '')
    """
    advanced_coded.py \\
        --input    ${long_tsvs} \\
        --config   ${config_json} \\
        --output   ${name}.diag.tsv

    phenotype_manifest.py \\
        --result   ${name}.diag.tsv \\
        --config   ${config_json} \\
        --output   ${name}.manifest.tsv
    """

    stub:
    def name = config_json.baseName.replace('.json', '')
    """
    touch ${name}.diag.tsv
    printf 'column_name\\tphenotype\\tdata_type\\tcode\\n' > ${name}.manifest.tsv
    """
}
