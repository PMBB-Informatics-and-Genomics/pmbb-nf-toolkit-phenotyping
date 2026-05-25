process RESOLVE_CONFIG {
    tag "resolve_config"
    publishDir "${params.output_dir}/configs", mode: 'copy'

    input:
    path biobank_config

    output:
    path "quantitative/*.json",             emit: quant_configs,         optional: true
    path "categorical/*.json",              emit: cat_configs,           optional: true
    path "diagnostic/*.json",              emit: diag_configs,          optional: true
    path "advanced_diagnostic/*.json",      emit: advanced_diag_configs, optional: true
    path "risk_classification_config.json", emit: risk_config,           optional: true

    script:
    """
    resolve_config.py \\
        --config ${biobank_config} \\
        --output_dir . \\
        --diagnostic_chunk_size ${params.diagnostic_chunk_size ?: 0}
    """

    stub:
    """
    mkdir -p quantitative categorical diagnostic advanced_diagnostic
    echo '{}' > quantitative/stub.json
    echo '{}' > categorical/stub.json
    echo '{}' > diagnostic/stub.json
    echo '{}' > advanced_diagnostic/stub.json
    echo '{"rules":[]}' > risk_classification_config.json
    """
}
