process RESOLVE_CONFIG {
    tag "resolve_config"
    publishDir "${params.output_dir}/configs", mode: 'copy'

    input:
    path biobank_config

    output:
    path "quantitative/*.json",             emit: quant_configs,         optional: true
    path "categorical/*.json",              emit: cat_configs,           optional: true
    path "coded/*.json",              emit: coded_configs,          optional: true
    path "advanced_coded/*.json",      emit: advanced_coded_configs, optional: true
    path "risk_classification_config.json", emit: risk_config,           optional: true

    script:
    """
    resolve_config.py \\
        --config ${biobank_config} \\
        --output_dir . \\
        --coded_chunk_size ${params.coded_chunk_size ?: 0}
    """

    stub:
    """
    mkdir -p quantitative categorical coded advanced_coded
    echo '{}' > quantitative/stub.json
    echo '{}' > categorical/stub.json
    echo '{}' > coded/stub.json
    echo '{}' > advanced_coded/stub.json
    echo '{"rules":[]}' > risk_classification_config.json
    """
}
