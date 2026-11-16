process EDA_PROFILE {
    tag "${group}"
    publishDir "${params.output_dir}/eda/profiles", mode: 'copy'
    // Exit 3 = profiler not installed: don't fail the run, but don't cache it
    // as a success either, so -resume retries once it is installed
    errorStrategy { task.exitStatus == 3 ? 'ignore' : 'terminate' }

    input:
    tuple val(group), path(result_tsvs, stageAs: 'results/*')
    path profile_selection
    path eda_config

    output:
    path "${group}.html", emit: html, optional: true

    // Optional: a group with no staged data writes nothing (and warns).
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
