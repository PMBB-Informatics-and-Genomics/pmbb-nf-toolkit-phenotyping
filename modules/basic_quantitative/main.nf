process BASIC_QUANTITATIVE {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy',
        saveAs: { it.endsWith('.manifest.tsv') ? null : it }

    input:
    tuple path(config_json), path(long_tsv)

    output:
    path "*.quant.tsv",             emit: result
    path "*.manifest.tsv",          emit: manifest
    path "quantitative_report.txt", emit: report, optional: true

    script:
    def name = config_json.baseName
    """
    basic_quantitative.py \\
        --input   ${long_tsv} \\
        --config  ${config_json} \\
        --output  ${name}.quant.tsv \\
        --report  quantitative_report.txt

    phenotype_manifest.py \\
        --result  ${name}.quant.tsv \\
        --config  ${config_json} \\
        --output  ${name}.manifest.tsv
    """

    stub:
    def name = config_json.baseName
    """
    touch ${name}.quant.tsv
    printf 'column_name\\tphenotype\\tdata_type\\tcode\\n' > ${name}.manifest.tsv
    """
}
