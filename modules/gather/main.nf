process GATHER {
    tag "gather"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path phenotype_tsvs

    output:
    path "phenotype_table.wide.tsv", emit: wide,        optional: true
    path "phenotype_table.long.tsv", emit: long_format, optional: true

    script:
    def fmt = params.gather_format ?: 'wide'
    """
    gather.py \\
        --inputs ${phenotype_tsvs} \\
        --output_prefix phenotype_table \\
        --format ${fmt}
    """

    stub:
    """
    touch phenotype_table.wide.tsv phenotype_table.long.tsv
    """
}
