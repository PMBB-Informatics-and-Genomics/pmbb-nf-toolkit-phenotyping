process BASIC_CATEGORICAL {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy',
        saveAs: { it.endsWith('.manifest.tsv') ? null : it }

    input:
    tuple path(config_json), path(long_tsv)

    output:
    path "*.cat.tsv",               emit: result
    path "*.manifest.tsv",          emit: manifest
    path "*.codes.tsv",             emit: codes, optional: true
    path "categorical_report.txt",  emit: report, optional: true

    script:
    def name = config_json.baseName
    """
    basic_categorical.py \\
        --input   ${long_tsv} \\
        --config  ${config_json} \\
        --output  ${name}.cat.tsv \\
        --codes   ${name}.codes.tsv \\
        --report  categorical_report.txt

    phenotype_manifest.py \\
        --result  ${name}.cat.tsv \\
        --config  ${config_json} \\
        --output  ${name}.manifest.tsv
    """

    stub:
    def name = config_json.baseName
    """
    touch ${name}.cat.tsv
    printf 'column_name\\tphenotype\\tdata_type\\tcode\\n' > ${name}.manifest.tsv
    """
}
