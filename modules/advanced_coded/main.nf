process ADVANCED_CODED {
    tag "${config_json.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy',
        saveAs: { it.endsWith('.manifest.tsv') ? null : it }

    input:
    tuple path(config_json), path(long_tsvs)
    path sample_list

    output:
    path "*.diag.tsv",     emit: result
    path "*.manifest.tsv", emit: manifest

    // sample_list is the roster: without it, samples with no source records stay NA
    script:
    def name = config_json.baseName.replace('.json', '')
    def sl_arg = sample_list ? "--all_samples ${sample_list}" : ''
    """
    advanced_coded.py \\
        --input    ${long_tsvs} \\
        --config   ${config_json} \\
        --output   ${name}.diag.tsv \\
        ${sl_arg}

    phenotype_manifest.py \\
        --result   ${name}.diag.tsv \\
        --config   ${config_json} \\
        --output   ${name}.manifest.tsv
    """

    stub:
    def name = config_json.baseName.replace('.json', '')
    """
    touch ${name}.diag.tsv
    printf 'column_name\\tphenotype\\tdata_type\\tcode\\n' > ${name}.manifest.tsv
    """
}
