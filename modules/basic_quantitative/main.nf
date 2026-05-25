process BASIC_QUANTITATIVE {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy'

    input:
    tuple path(config_json), path(long_tsv)

    output:
    path "*.quant.tsv",             emit: result
    path "quantitative_report.txt", emit: report, optional: true

    script:
    def name = config_json.baseName
    """
    basic_quantitative.py \\
        --input   ${long_tsv} \\
        --config  ${config_json} \\
        --output  ${name}.quant.tsv \\
        --report  quantitative_report.txt
    """

    stub:
    def name = config_json.baseName
    """
    touch ${name}.quant.tsv
    """
}
