process PREPROCESS {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/long_format", mode: 'copy'

    input:
    path config_json
    path sample_list

    output:
    tuple path(config_json), path("*.long.tsv"), emit: with_config
    path "preprocess_report_${config_json.baseName}.txt", emit: report, optional: true

    script:
    def sl_arg = sample_list ? "--sample_list ${sample_list}" : ''
    """
    preprocess.py \\
        --config ${config_json} \\
        --output_dir . \\
        ${sl_arg}
    """

    stub:
    """
    touch ${config_json.baseName}.long.tsv
    """
}
