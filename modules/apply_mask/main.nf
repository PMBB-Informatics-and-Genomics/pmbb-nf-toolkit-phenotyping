process APPLY_MASK {
    tag "${config_json.baseName}"

    input:
    tuple path(config_json), path(long_tsv)
    path mask_intervals

    output:
    tuple path(config_json), path("${long_tsv.baseName}.masked.long.tsv"), emit: with_config

    script:
    """
    apply_mask.py \\
        --input          ${long_tsv} \\
        --config         ${config_json} \\
        --mask_intervals ${mask_intervals} \\
        --output         ${long_tsv.baseName}.masked.long.tsv
    """

    stub:
    """
    cp ${long_tsv} ${long_tsv.baseName}.masked.long.tsv
    """
}
