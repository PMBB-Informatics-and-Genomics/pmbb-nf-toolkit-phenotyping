#!/usr/bin/env python3
"""
basic_coded.py — Binary presence/absence for a single coded code.

Usage:
    basic_coded.py --input N80.0.long.tsv --config diag.json \
                        --output N80.0.diag.tsv [--all_samples samples.txt]
    basic_coded.py --input 401.1.long.tsv --config phecode.json \
                        --output_dir .              # → {column_prefix}{code}.diag.tsv

If the config sets column_prefix (e.g. "phe_"), the output column and the
--output_dir filename become {column_prefix}{code}. The input long.tsv is never
modified — its 'value' column keeps the raw code for advanced_coded/event_mask.
"""

import argparse
import json
import os
import sys
import pandas as pd


def run_basic_coded(long_tsv, config_path, output_path=None, all_samples=None,
                    output_dir=None):
    """
    Compute binary presence for one coded code.

    Parameters
    ----------
    long_tsv      : path to {code}.long.tsv (canonical long format, single code)
    config_path   : path to JSON with min_occurrences, missing_as_control, column_prefix
    output_path   : path to write the diag TSV (takes precedence over output_dir)
    all_samples   : optional list of all sample IDs for roster (used in tests)
    output_dir    : if output_path is not given, write {column}.diag.tsv here

    Returns the output column name.
    """
    with open(config_path) as f:
        cfg = json.load(f)

    min_occ = int(cfg.get('min_occurrences', 1))
    missing_as_control = bool(cfg.get('missing_as_control', False))
    column_prefix = cfg.get('column_prefix') or ''

    df = pd.read_csv(long_tsv, sep='\t', dtype=str)

    # Derive code name from the phenotype column
    code = df['phenotype'].iloc[0] if len(df) > 0 else os.path.basename(long_tsv).replace('.long.tsv', '')
    col = f'{column_prefix}{code}'
    safe_col = col.replace('/', '_').replace('\\', '_')
    if output_path is None:
        output_path = os.path.join(output_dir or '.', f'{safe_col}.diag.tsv')

    # Count occurrences per sample
    counts = df.groupby('sample_ID').size()

    # Build roster
    if all_samples is None:
        all_samples = counts.index.tolist()
    result = pd.DataFrame({'sample_ID': list(all_samples)})

    # Assign values (vectorized)
    result[col] = pd.NA
    present_mask = result['sample_ID'].isin(counts[counts >= min_occ].index)
    result.loc[present_mask, col] = 1
    if missing_as_control:
        absent_mask = ~result['sample_ID'].isin(counts.index)
        result.loc[absent_mask, col] = 0

    result[col] = result[col].astype('Int64')
    result.to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path}', file=sys.stderr)

    return col


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True, help='Single-code long-format TSV')
    p.add_argument('--config', required=True, help='Parent coded JSON config')
    out = p.add_mutually_exclusive_group(required=True)
    out.add_argument('--output', help='Output diag TSV path')
    out.add_argument('--output_dir',
                     help='Write {column_prefix}{code}.diag.tsv into this directory')
    p.add_argument('--all_samples', default=None,
                   help='Optional newline-separated file of all sample IDs')
    args = p.parse_args()

    all_samples = None
    if args.all_samples and os.path.isfile(args.all_samples):
        with open(args.all_samples) as f:
            all_samples = [line.strip() for line in f if line.strip()]

    run_basic_coded(args.input, args.config, args.output, all_samples,
                    output_dir=args.output_dir)


if __name__ == '__main__':
    main()
