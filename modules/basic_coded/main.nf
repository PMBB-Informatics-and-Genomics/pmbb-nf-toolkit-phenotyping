process BASIC_CODED {
    tag "${long_tsv.baseName}"
    publishDir "${params.output_dir}/phenotypes", mode: 'copy',
        saveAs: { it.endsWith('.manifest.tsv') ? null : it }

    input:
    tuple path(config_json), path(long_tsv)
    path sample_list

    output:
    path "*.diag.tsv",     emit: result
    path "*.manifest.tsv", emit: manifest

    // Output filename is {column_prefix}{code}.diag.tsv, derived by basic_coded.py
    // from the config's column_prefix. The input long.tsv keeps the raw code.
    // sample_list is the roster: without it, missing_as_control cannot write 0s.
    script:
    def sl_arg = sample_list ? "--all_samples ${sample_list}" : ''
    """
    basic_coded.py \\
        --input      ${long_tsv} \\
        --config     ${config_json} \\
        --output_dir . \\
        ${sl_arg}

    for f in *.diag.tsv; do
        phenotype_manifest.py \\
            --result "\$f" \\
            --config ${config_json} \\
            --output "\${f%.diag.tsv}.manifest.tsv"
    done
    """

    stub:
    def code = long_tsv.baseName.replace('.long', '')
    """
    touch ${code}.diag.tsv
    printf 'column_name\\tphenotype\\tdata_type\\tcode\\n' > ${code}.manifest.tsv
    """
}
