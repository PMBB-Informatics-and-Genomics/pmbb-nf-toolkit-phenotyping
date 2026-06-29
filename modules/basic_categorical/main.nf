process BASIC_CATEGORICAL {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy'

    input:
    tuple path(config_json), path(long_tsv)

    output:
    path "*.cat.tsv",               emit: result
    path "categorical_report.txt",  emit: report, optional: true

    script:
    def name = config_json.baseName
    """
    basic_categorical.py \\
        --input   ${long_tsv} \\
        --config  ${config_json} \\
        --output  ${name}.cat.tsv \\
        --report  categorical_report.txt
    """

    stub:
    def name = config_json.baseName
    """
    touch ${name}.cat.tsv
    """
}
