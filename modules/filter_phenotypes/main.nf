process FILTER_PHENOTYPES {
    tag "filter"
    publishDir "${params.output_dir}/filtered", mode: 'copy'

    input:
    path wide_table
    path filter_config, stageAs: 'input/filter.yaml'
    path covariates
    path summary
    path manifest

    output:
    path "phenotype_table.filtered.wide.tsv", emit: wide
    path "filter_log.tsv",                    emit: log
    path "phenotype_manifest.filtered.tsv",   emit: manifest
    path "filter_config_used.yaml",           emit: config_used

    // Writes new files only — phenotype_table.wide.tsv is never modified.
    // The filter config is republished next to the outputs for provenance.
    script:
    """
    cp input/filter.yaml filter_config_used.yaml
    filter_phenotypes.py \\
        --config     input/filter.yaml \\
        --input      ${wide_table} \\
        --covariates ${covariates} \\
        --summary    ${summary} \\
        --manifest   ${manifest} \\
        --output_dir .
    """

    stub:
    """
    touch phenotype_table.filtered.wide.tsv filter_log.tsv phenotype_manifest.filtered.tsv \\
          filter_config_used.yaml
    """
}
