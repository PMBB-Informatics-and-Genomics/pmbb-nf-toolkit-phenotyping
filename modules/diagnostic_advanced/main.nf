process DIAGNOSTIC_ADVANCED {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy'

    input:
    tuple path(config_json), path(long_tsvs)

    output:
    path "*.diag.tsv", emit: result

    script:
    def name = config_json.baseName.replace('.json', '')
    """
    diagnostic_advanced.py \\
        --input    ${long_tsvs} \\
        --config   ${config_json} \\
        --output   ${name}.diag.tsv
    """

    stub:
    def name = config_json.baseName.replace('.json', '')
    """
    touch ${name}.diag.tsv
    """
}
