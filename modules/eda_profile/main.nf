process EDA_PROFILE {
    tag "${group}"
    publishDir "${params.output_dir}/eda/profiles", mode: 'copy'

    input:
    tuple val(group), path(result_tsvs, stageAs: 'results/*')
    path profile_selection
    path eda_config

    output:
    path "${group}.html", emit: html, optional: true

    // Writes nothing (and warns) if ydata-profiling is not installed.
    script:
    """
    eda_profile.py \\
        --group      '${group}' \\
        --inputs     results/* \\
        --selection  ${profile_selection} \\
        --eda_config ${eda_config} \\
        --output     '${group}.html'
    """

    stub:
    """
    touch '${group}.html'
    """
}
