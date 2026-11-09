#!/usr/bin/env python3
"""
Split a long-format TSV into N chunk files partitioned by sample_ID.

All rows for a given sample_ID are guaranteed to stay in the same chunk.
Chunks are written as chunk_001.tsv, chunk_002.tsv, … into --output_dir.

Usage:
  split_samples.py --input long_format.tsv --chunk_size 10000 --output_dir chunks/
"""

import argparse
import os
import sys

import pandas as pd

from utils import resolve_sep


def parse_args():
    p = argparse.ArgumentParser(
        description='Split a long-format TSV into sample-partitioned chunks.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--input',       required=True,  help='Long-format TSV input file.')
    p.add_argument('--chunk_size',  type=int, required=True,
                   help='Maximum number of unique samples per chunk.')
    p.add_argument('--output_dir',  default='.',    help='Directory to write chunk files.')
    p.add_argument('--sample_col',  default='sample_ID',
                   help='Name of the sample identifier column.')
    return p.parse_args()


def split_samples(input_path, chunk_size, output_dir, sample_col='sample_ID'):
    """
    Read input_path, split by sample_col into chunks of up to chunk_size samples,
    write chunk files to output_dir.  Returns list of written file paths.
    """
    if chunk_size < 1:
        raise ValueError(f'chunk_size must be ≥ 1, got {chunk_size}')

    sep, engine = resolve_sep(filepath=input_path)
    df  = pd.read_csv(input_path, sep=sep, engine=engine, dtype=str)

    if sample_col not in df.columns:
        raise ValueError(f"Column '{sample_col}' not found in {input_path}. "
                         f"Available columns: {list(df.columns)}")

    samples   = sorted(df[sample_col].dropna().unique().tolist())
    n_samples = len(samples)
    n_chunks  = max(1, (n_samples + chunk_size - 1) // chunk_size)
    n_digits  = len(str(n_chunks))

    os.makedirs(output_dir, exist_ok=True)
    written = []

    for i in range(n_chunks):
        chunk_samples = set(samples[i * chunk_size : (i + 1) * chunk_size])
        chunk_df      = df[df[sample_col].isin(chunk_samples)].copy()
        chunk_id      = str(i + 1).zfill(n_digits)
        out_path      = os.path.join(output_dir, f'chunk_{chunk_id}.tsv')

        chunk_df.to_csv(out_path, sep='\t', index=False)
        written.append(out_path)
        print(f'  chunk_{chunk_id}: {len(chunk_samples)} samples, '
              f'{len(chunk_df)} rows → {out_path}', file=sys.stderr)

    print(f'Split {n_samples} samples into {n_chunks} chunks '
          f'(chunk_size={chunk_size}) → {output_dir}', file=sys.stderr)
    return written


def main():
    args = parse_args()
    split_samples(
        input_path  = args.input,
        chunk_size  = args.chunk_size,
        output_dir  = args.output_dir,
        sample_col  = args.sample_col,
    )


if __name__ == '__main__':
    main()
