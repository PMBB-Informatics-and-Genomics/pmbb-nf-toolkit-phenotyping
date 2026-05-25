#!/usr/bin/env python3
"""
diagnostic_basic.py — Binary presence/absence for a single diagnostic code.

Usage:
    diagnostic_basic.py --input N80.0.long.tsv --config diag.json \
                        --output N80.0.diag.tsv [--all_samples samples.txt]
"""

import argparse
import json
import os
import sys
import pandas as pd


def run_diagnostic_basic(long_tsv, config_path, output_path, all_samples=None):
    """
    Compute binary presence for one diagnostic code.

    Parameters
    ----------
    long_tsv     : path to {code}.long.tsv (canonical long format, single code)
    config_path  : path to JSON with min_occurrences, missing_as_control
    output_path  : path to write {code}.diag.tsv
    all_samples  : optional list of all sample IDs for roster (used in tests)
    """
    with open(config_path) as f:
        cfg = json.load(f)

    min_occ = int(cfg.get('min_occurrences', 1))
    missing_as_control = bool(cfg.get('missing_as_control', False))

    df = pd.read_csv(long_tsv, sep='\t', dtype=str)

    # Derive code name from the phenotype column
    code = df['phenotype'].iloc[0] if len(df) > 0 else os.path.basename(long_tsv).replace('.long.tsv', '')

    # Count occurrences per sample
    counts = df.groupby('sample_ID').size()

    # Build roster
    if all_samples is None:
        all_samples = counts.index.tolist()
    result = pd.DataFrame({'sample_ID': list(all_samples)})

    # Assign values (vectorized)
    result[code] = pd.NA
    present_mask = result['sample_ID'].isin(counts[counts >= min_occ].index)
    result.loc[present_mask, code] = 1
    if missing_as_control:
        absent_mask = ~result['sample_ID'].isin(counts.index)
        result.loc[absent_mask, code] = 0

    result[code] = result[code].astype('Int64')
    result.to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True, help='Single-code long-format TSV')
    p.add_argument('--config', required=True, help='Parent diagnostic JSON config')
    p.add_argument('--output', required=True, help='Output diag TSV')
    p.add_argument('--all_samples', default=None,
                   help='Optional newline-separated file of all sample IDs')
    args = p.parse_args()

    all_samples = None
    if args.all_samples and os.path.isfile(args.all_samples):
        with open(args.all_samples) as f:
            all_samples = [line.strip() for line in f if line.strip()]

    run_diagnostic_basic(args.input, args.config, args.output, all_samples)


if __name__ == '__main__':
    main()
