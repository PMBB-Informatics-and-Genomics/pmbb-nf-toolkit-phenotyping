#!/usr/bin/env python3
"""
Merge multiple wide-format TSV chunk files into a single output file.

Each chunk file has one row per sample and the same set of columns
(sample_ID plus phenotype columns).  Files are row-concatenated and
sorted by sample_ID.  Missing columns in any chunk are filled with NaN.

Usage:
  merge_wide.py --inputs chunk_001.tsv chunk_002.tsv … --output merged.tsv

  # Shell glob expansion also works:
  merge_wide.py --inputs result_chunk_*.tsv --output merged.tsv
"""

import argparse
import glob as _glob
import os
import sys

import pandas as pd

from utils import resolve_sep


def parse_args():
    p = argparse.ArgumentParser(
        description='Row-concatenate wide-format TSV chunk files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--inputs', nargs='+', required=True,
                   help='Input TSV files (space-separated or shell-expanded glob).')
    p.add_argument('--output', required=True, help='Output merged TSV file.')
    p.add_argument('--sort_by', default='sample_ID',
                   help='Column to sort the output by.')
    p.add_argument('--no_sort', action='store_true',
                   help='Skip sorting (preserves concatenation order).')
    return p.parse_args()


def merge_wide(input_files, output_path, sort_by='sample_ID', no_sort=False):
    """
    Concatenate wide-format TSVs row-wise.  Returns the merged DataFrame.
    """
    # Resolve any glob patterns that the shell didn't expand (e.g. when
    # called directly from Python without a shell).
    resolved = []
    for f in input_files:
        expanded = _glob.glob(f)
        resolved.extend(expanded if expanded else [f])
    resolved = sorted(set(resolved))

    if not resolved:
        raise FileNotFoundError(f'No input files found from: {input_files}')

    missing = [f for f in resolved if not os.path.isfile(f)]
    if missing:
        raise FileNotFoundError(f'Input files not found: {missing}')

    chunks = []
    for f in resolved:
        sep, engine = resolve_sep(filepath=f)
        df  = pd.read_csv(f, sep=sep, engine=engine, dtype=str)
        chunks.append(df)
        print(f'  loaded {f}: {len(df)} rows', file=sys.stderr)

    merged = pd.concat(chunks, axis=0, ignore_index=True)

    if not no_sort and sort_by in merged.columns:
        merged = merged.sort_values(sort_by).reset_index(drop=True)

    merged.to_csv(output_path, sep='\t', index=False)
    print(f'Merged {len(resolved)} files → {output_path} '
          f'({len(merged)} rows × {len(merged.columns)} cols)', file=sys.stderr)
    return merged


def main():
    args = parse_args()
    merge_wide(
        input_files = args.inputs,
        output_path = args.output,
        sort_by     = args.sort_by,
        no_sort     = args.no_sort,
    )


if __name__ == '__main__':
    main()
