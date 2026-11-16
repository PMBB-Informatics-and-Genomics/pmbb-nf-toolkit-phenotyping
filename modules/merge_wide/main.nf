// MERGE_WIDE — row-concatenate wide-format TSV chunk files into a single output.
//
// Each chunk file has one row per sample; this process stacks them vertically
// and sorts by sample_ID.  Columns that appear in some chunks but not others
// are filled with NaN in the merged output.
//
// This process is intended to be included with an alias so it can be reused
// for different module outputs without name collisions:
//
//   include { MERGE_WIDE as MERGE_QUANTITATIVE } from './modules/merge_wide/main.nf'
//   include { MERGE_WIDE as MERGE_CATEGORICAL  } from './modules/merge_wide/main.nf'
//
// Inputs:
//   chunk_results — all chunk TSV files (pass via .collect())
//   label         — output file stem, e.g. "quantitative_summary"
//   subdir        — publishDir subdirectory, e.g. "quantitative"
//
// Outputs:
//   {label}.tsv   — merged wide-format TSV, one row per sample

process MERGE_WIDE {
    tag "merge_${label}"
    publishDir "${params.output_dir}/${subdir}", mode: 'copy'

    input:
    path chunk_results   // collected chunk TSV files (staged flat in work dir)
    val  label           // output filename stem
    val  subdir          // publishDir sub-path

    output:
    path "${label}.tsv", emit: merged

    script:
    """
    merge_wide.py \\
        --inputs *.tsv \\
        --output ${label}.tsv
    """

    stub:
    """
    touch ${label}.tsv
    """
}
