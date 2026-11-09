#!/usr/bin/env python3
"""
apply_mask.py — Filter a phenotype long.tsv using precomputed mask intervals.

For each sample in the long.tsv, removes rows whose occurrence_date falls within
any mask window defined in mask_intervals.tsv for the event types listed in the
phenotype's config.

IMPORTANT: If the phenotype config has date_col=null (no real dates — occurrence_date
was filled from reference_date), masking is SKIPPED and the input is passed through
unchanged. Masking is meaningless when all dates are synthetic. A WARNING is logged.
This behavior is by design — see event_masking docs.

Overlapping windows from multiple event types are merged before filtering.

Usage:
    apply_mask.py \\
        --input HbA1c.long.tsv \\
        --config HbA1c.json \\
        --mask_intervals mask_intervals.tsv \\
        --output HbA1c.masked.long.tsv
"""

import argparse
import json
import os
import shutil
import sys

import pandas as pd


def _merge_intervals(intervals: list) -> list:
    """Merge overlapping (date, date) pairs. Input need not be sorted."""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [list(sorted_ivs[0])]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def _build_sample_windows(mask_df: pd.DataFrame, event_types: list) -> dict:
    """
    Filter mask_df to event_types, then merge all windows per sample.
    Returns {sample_ID: [(start_date, end_date), ...]} with merged intervals.
    Duplicate rows are collapsed via interval merging.
    """
    filtered = mask_df[mask_df['event_type'].isin(event_types)].copy()
    if filtered.empty:
        return {}

    filtered['window_start'] = pd.to_datetime(filtered['window_start'], errors='coerce').dt.date
    filtered['window_end'] = pd.to_datetime(filtered['window_end'], errors='coerce').dt.date
    filtered = filtered.dropna(subset=['window_start', 'window_end'])

    result = {}
    for sid, grp in filtered.groupby('sample_ID'):
        intervals = list(zip(grp['window_start'], grp['window_end']))
        result[sid] = _merge_intervals(intervals)
    return result


def _mask_rows(df: pd.DataFrame, sample_windows: dict) -> pd.DataFrame:
    """
    Return df with rows removed where occurrence_date falls within the sample's
    mask windows (inclusive on both endpoints). Rows with unparseable dates are kept.
    """
    if not sample_windows:
        return df

    # Build a windows DataFrame for vectorized merge
    win_rows = [
        {'sample_ID': sid, '_ws': s, '_we': e}
        for sid, windows in sample_windows.items()
        for s, e in windows
    ]
    win_df = pd.DataFrame(win_rows)

    df2 = df.assign(
        _date=pd.to_datetime(df['occurrence_date'], errors='coerce').dt.date,
        _idx=df.index,
    )
    merged = df2.merge(win_df, on='sample_ID', how='left')
    in_window = (merged['_date'] >= merged['_ws']) & (merged['_date'] <= merged['_we'])
    drop_idx = set(merged.loc[in_window, '_idx'])

    n_unparseable = df2['_date'].isna().sum()
    if n_unparseable > 0:
        print(
            f'WARNING: {n_unparseable} rows have unparseable occurrence_date and will not be masked.',
            file=sys.stderr,
        )

    return df.loc[~df.index.isin(drop_idx)].reset_index(drop=True)


def run_apply_mask(
    long_tsv_path: str,
    config_path: str,
    mask_intervals_path: str,
    output_path: str,
) -> None:
    """
    Apply event mask intervals to a phenotype long.tsv.

    Parameters
    ----------
    long_tsv_path      : path to input long.tsv (canonical long format)
    config_path        : path to per-phenotype JSON (must contain date_col + event_mask)
    mask_intervals_path: path to mask_intervals.tsv from EVENT_MASK
    output_path        : path to write filtered long.tsv
    """
    with open(config_path) as f:
        cfg = json.load(f)

    event_mask_types = cfg.get('event_mask') or []

    # Skip masking entirely if no real dates available
    if not cfg.get('date_col'):
        print(
            f'WARNING: apply_mask skipped for "{cfg.get("phenotype_name", "?")}": '
            f'date_col is null — occurrence_date is synthetic (reference_date). '
            f'Event masking requires real observation dates.',
            file=sys.stderr,
        )
        shutil.copy2(long_tsv_path, output_path)
        return

    try:
        mask_df = pd.read_csv(mask_intervals_path, sep='\t', dtype=str)
    except pd.errors.EmptyDataError:
        mask_df = pd.DataFrame(columns=['sample_ID', 'event_type', 'window_start', 'window_end'])

    required_mask_cols = {'sample_ID', 'event_type', 'window_start', 'window_end'}
    if not required_mask_cols.issubset(mask_df.columns):
        missing = required_mask_cols - set(mask_df.columns)
        raise ValueError(f'mask_intervals.tsv is missing columns: {missing}')

    sample_windows = _build_sample_windows(mask_df, event_mask_types)

    df = pd.read_csv(long_tsv_path, sep='\t', dtype=str)

    required_long_cols = {'sample_ID', 'occurrence_date'}
    if not required_long_cols.issubset(df.columns):
        missing = required_long_cols - set(df.columns)
        raise ValueError(f'Input long.tsv is missing columns: {missing}')

    n_before = len(df)
    df = _mask_rows(df, sample_windows)
    n_removed = n_before - len(df)

    df.to_csv(output_path, sep='\t', index=False)
    print(
        f'apply_mask: {cfg.get("phenotype_name", "?")} — '
        f'{n_removed} rows removed ({n_before - n_removed} remaining) → {output_path}',
        file=sys.stderr,
    )


def main():
    p = argparse.ArgumentParser(description='Filter long.tsv rows within event mask windows')
    p.add_argument('--input', required=True, help='Input long.tsv')
    p.add_argument('--config', required=True, help='Per-phenotype JSON config')
    p.add_argument('--mask_intervals', required=True, help='mask_intervals.tsv from event_mask')
    p.add_argument('--output', required=True, help='Output filtered long.tsv')
    args = p.parse_args()
    run_apply_mask(args.input, args.config, args.mask_intervals, args.output)


if __name__ == '__main__':
    main()
