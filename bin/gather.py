#!/usr/bin/env python3
"""
gather.py — Outer-join per-phenotype TSVs into a final wide or long table.

Usage:
    gather.py --inputs *.quant.tsv *.cat.tsv *.diag.tsv \
              --output_prefix phenotype_table \
              --format wide
"""

import argparse
import os
import sys
import pandas as pd


def run_gather(input_paths, output_prefix, fmt='wide'):
    """
    Outer-join all per-phenotype TSVs on sample_ID.

    Parameters
    ----------
    input_paths   : list of TSV file paths; each must have a 'sample_ID' column
    output_prefix : output file prefix (without extension)
    fmt           : 'wide', 'long', or 'both'
    """
    if not input_paths:
        print('WARNING: no input files provided to gather', file=sys.stderr)
        return

    dfs = []
    for path in input_paths:
        df = pd.read_csv(path, sep='\t')
        if 'sample_ID' not in df.columns:
            print(f'WARNING: {path} has no sample_ID column — skipping', file=sys.stderr)
            continue
        dfs.append(df)

    if not dfs:
        raise ValueError('No valid input files — all inputs missing sample_ID column')

    # Option 3: long-only — stack-and-melt each file individually; avoids
    # building the wide table and never produces NaN cross-product rows.
    if fmt == 'long':
        long_parts = [
            df.melt(id_vars='sample_ID', var_name='phenotype', value_name='value')
            for df in dfs
        ]
        long = pd.concat(long_parts, axis=0, ignore_index=True)
        long = long.sort_values(['sample_ID', 'phenotype']).reset_index(drop=True)
        out_path = f'{output_prefix}.long.tsv'
        long.to_csv(out_path, sep='\t', index=False)
        print(f'Wrote: {out_path}', file=sys.stderr)
        return

    # Option 2: outer-join via a single pd.concat on indexed DataFrames;
    # avoids O(N) intermediate allocations from sequential merges.
    dfs_indexed = [df.set_index('sample_ID') for df in dfs]
    wide = pd.concat(dfs_indexed, axis=1, join='outer').reset_index()
    wide = wide.sort_values('sample_ID').reset_index(drop=True)

    if fmt in ('wide', 'both'):
        out_path = f'{output_prefix}.wide.tsv'
        wide.to_csv(out_path, sep='\t', index=False)
        print(f'Wrote: {out_path}', file=sys.stderr)

    if fmt == 'both':
        # Option 3 for 'both': generate long independently via stack-and-melt
        # so the long file is consistent with fmt='long' behavior.
        long_parts = [
            df.melt(id_vars='sample_ID', var_name='phenotype', value_name='value')
            for df in dfs
        ]
        long = pd.concat(long_parts, axis=0, ignore_index=True)
        long = long.sort_values(['sample_ID', 'phenotype']).reset_index(drop=True)
        out_path = f'{output_prefix}.long.tsv'
        long.to_csv(out_path, sep='\t', index=False)
        print(f'Wrote: {out_path}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inputs', nargs='+', required=True, help='Per-phenotype TSV files')
    p.add_argument('--output_prefix', required=True)
    p.add_argument('--format', default='wide', choices=['wide', 'long', 'both'])
    args = p.parse_args()
    run_gather(args.inputs, args.output_prefix, fmt=args.format)


if __name__ == '__main__':
    main()
