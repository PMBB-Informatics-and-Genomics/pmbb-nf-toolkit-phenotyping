process DIAGNOSTIC_BASIC {
    tag "${long_tsv.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy'

    input:
    tuple path(config_json), path(long_tsv)

    output:
    path "*.diag.tsv", emit: result

    script:
    def code = long_tsv.baseName.replace('.long', '')
    """
    diagnostic_basic.py \\
        --input    ${long_tsv} \\
        --config   ${config_json} \\
        --output   ${code}.diag.tsv
    """

    stub:
    def code = long_tsv.baseName.replace('.long', '')
    """
    touch ${code}.diag.tsv
    """
}
