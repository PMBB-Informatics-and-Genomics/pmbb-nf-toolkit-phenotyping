process EVENT_MASK {
    tag "${source_name}"
    publishDir "${params.output_dir}/event_mask", mode: 'copy'

    input:
    tuple val(source_name), path(source_tsvs)
    path event_mask_config

    output:
    path "${source_name}.mask_intervals.tsv", emit: intervals

    script:
    def tsv_list = (source_tsvs instanceof List ? source_tsvs : [source_tsvs]).join(' ')
    """
    event_mask.py \\
        --config   ${event_mask_config} \\
        --source   ${source_name} \\
        --input    ${tsv_list} \\
        --output   ${source_name}.mask_intervals.tsv
    """

    stub:
    """
    touch ${source_name}.mask_intervals.tsv
    """
}
