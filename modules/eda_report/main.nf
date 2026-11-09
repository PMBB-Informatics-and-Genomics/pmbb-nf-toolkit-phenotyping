process EDA_REPORT {
    tag "eda_report"
    publishDir "${params.output_dir}/eda", mode: 'copy',
        saveAs: { it.startsWith('manifest/') ? null : it }
    // With --run_eda the manifest at the output root is this flagged copy
    publishDir "${params.output_dir}", mode: 'copy', pattern: 'manifest/*',
        saveAs: { it.replaceFirst('manifest/', '') }

    input:
    path eda_stats
    path manifest, stageAs: 'input/phenotype_manifest.tsv'
    path eda_config
    // Every EDA_STATS batch emits eda_cases.batch.tsv; index them to avoid collisions
    path case_tsvs, stageAs: 'cases/cases_*.tsv'
    path covariates, stageAs: 'input/eda_covariates.tsv'

    output:
    path "eda_summary.tsv",                 emit: summary
    path "eda_flags.tsv",                   emit: flags
    path "eda_report.html",                 emit: html
    path "profile_selection.tsv",           emit: profile_selection
    path "restriction_suggestions.tsv",     emit: restrictions
    path "filter_suggested.yaml",           emit: filter_suggested
    path "eda_covariates.tsv",              emit: covariates
    path "manifest/phenotype_manifest.tsv", emit: manifest
    path "pairwise_overlap.tsv",            emit: pairwise, optional: true

    script:
    def cases_arg = case_tsvs ? "--cases cases/*" : ''
    """
    mkdir -p manifest
    cp input/eda_covariates.tsv eda_covariates.tsv
    eda_report.py \\
        --stats           ${eda_stats} \\
        --manifest        input/phenotype_manifest.tsv \\
        --eda_config      ${eda_config} \\
        --output_dir      . \\
        --manifest_output manifest/phenotype_manifest.tsv \\
        ${cases_arg}
    """

    stub:
    """
    mkdir -p manifest
    touch eda_summary.tsv eda_flags.tsv eda_report.html restriction_suggestions.tsv \\
          filter_suggested.yaml eda_covariates.tsv manifest/phenotype_manifest.tsv
    printf 'column_name\\tsource_file\\tphenotype\\tgroup\\treason\\n' > profile_selection.tsv
    """
}
