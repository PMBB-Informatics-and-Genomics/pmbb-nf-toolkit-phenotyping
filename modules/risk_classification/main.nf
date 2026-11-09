process RISK_CLASSIFICATION {
    tag "risk_classification"
    publishDir "${params.output_dir}/risk_classification", mode: 'copy'

    input:
    path wide_tsv
    path risk_config

    output:
    path "risk_matrix.tsv",  emit: matrix
    path "risk_report.txt",  emit: report

    script:
    """
    risk_classification.py \\
        --wide_tsv ${wide_tsv} \\
        --risk_config ${risk_config} \\
        --output risk_matrix.tsv \\
        --report risk_report.txt
    """

    stub:
    """
    touch risk_matrix.tsv risk_report.txt
    """
}
