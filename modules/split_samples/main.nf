// SPLIT_SAMPLES — partition a long-format TSV by sample_ID into N chunk files.
//
// All rows for a given sample_ID are guaranteed to land in the same chunk,
// so downstream processes can treat each chunk as a fully independent sub-cohort.
//
// Inputs:
//   long_format  — canonical long-format TSV
//   chunk_size   — maximum samples per chunk (val, not path)
//
// Outputs:
//   chunks/chunk_*.tsv  — emitted as a list; use .flatten() to process one-per-job

process SPLIT_SAMPLES {
    tag "split_samples(chunk_size=${chunk_size})"

    input:
    path long_format
    val  chunk_size

    output:
    path "chunks/chunk_*.tsv", emit: chunks

    script:
    """
    mkdir -p chunks
    split_samples.py \\
        --input      ${long_format} \\
        --chunk_size ${chunk_size} \\
        --output_dir chunks/
    """

    stub:
    """
    mkdir -p chunks
    touch chunks/chunk_0.tsv
    """
}
