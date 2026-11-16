process EDA_STATS {
    tag "batch_${task.index}"

    input:
    path result_tsvs, stageAs: 'results/*'
    path covariate_tsvs, stageAs: 'covariates/*'
    path eda_config
    path sample_list

    output:
    path "eda_stats.batch.tsv", emit: stats
    path "eda_cases.batch.tsv", emit: cases
    path "eda_covariates.tsv",  emit: covariates

    // Covariates are staged separately: a covariate is also a result, so the
    // same file may arrive on both inputs.
    script:
    def cov_arg = covariate_tsvs ? "--covariates covariates/*" : ''
    def sl_arg  = sample_list ? "--universe ${sample_list}" : ''
    """
    eda_stats.py \\
        --inputs       results/* \\
        --eda_config   ${eda_config} \\
        --output       eda_stats.batch.tsv \\
        --cases_output eda_cases.batch.tsv \\
        --covariates_output eda_covariates.tsv \\
        ${cov_arg} \\
        ${sl_arg}
    """

    stub:
    """
    printf 'column_name\\tsource_file\\n' > eda_stats.batch.tsv
    printf 'column_name\\tsource_file\\tcase_ids\\n' > eda_cases.batch.tsv
    printf 'sample_ID\\n' > eda_covariates.tsv
    """
}
