#!/usr/bin/env python3
"""
eda_profile.py — ydata-profiling report for one profile group.

Usage:
    eda_profile.py --group ICD10 --inputs E11.9.diag.tsv I10.diag.tsv \
                   --selection profile_selection.tsv \
                   --eda_config eda_config.json --output ICD10.html

Loads only the selected columns of the group (from profile_selection.tsv,
written by eda_report.py), outer-joins them on sample_ID, downsamples to
eda.profile_max_samples and writes a minimal-mode ydata report.

Downsampling keeps cases: if any selected binary column is rare, up to half
the sample budget goes to samples that are a case in at least one such
column; the rest is a uniform random draw. Metrics in eda_summary.tsv always
use every sample — only the profile is downsampled.

The profiler is fg-data-profiling (`import data_profiling`, the maintained
successor of ydata-profiling), falling back to ydata-profiling. It is
optional: if neither is installed it exits with EXIT_NO_PROFILER, which
EDA_PROFILE ignores — the rest of the EDA still completes, and since the task
is not cached as a success, -resume retries it once the profiler is installed.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

# Exit status when no profiler is installed; EDA_PROFILE's errorStrategy ignores it
EXIT_NO_PROFILER = 3


def load_group(group, input_paths, selection_path):
    sel = pd.read_csv(selection_path, sep='\t', dtype=str)
    sel = sel[sel['group'] == group]
    by_file = {os.path.basename(p): p for p in input_paths}
    frames = []
    for source_file, sub in sel.groupby('source_file'):
        path = by_file.get(source_file)
        if path is None:
            print(f'WARNING: {source_file} not staged — skipping', file=sys.stderr)
            continue
        cols = ['sample_ID'] + list(dict.fromkeys(sub['column_name']))
        df = pd.read_csv(path, sep='\t', dtype={'sample_ID': str},
                         usecols=lambda c: c in cols)
        frames.append(df.drop_duplicates('sample_ID').set_index('sample_ID'))
    if not frames:
        return pd.DataFrame()
    wide = pd.concat(frames, axis=1, join='outer')
    return wide.loc[:, ~wide.columns.duplicated()]


def downsample(df, max_samples, seed=42):
    """Case-preserving downsample to at most max_samples rows."""
    if max_samples <= 0 or len(df) <= max_samples:
        return df
    rng = np.random.default_rng(seed)
    binary = [c for c in df.columns
              if set(pd.to_numeric(df[c], errors='coerce').dropna().unique()) <= {0, 1}]
    rare = [c for c in binary if (df[c] == 1).mean() < 0.5]
    keep = pd.Index([])
    if rare:
        case_rows = df.index[(df[rare] == 1).any(axis=1)]
        n_case = min(len(case_rows), max_samples // 2)
        keep = pd.Index(rng.choice(case_rows, size=n_case, replace=False))
    rest = df.index.difference(keep)
    fill = rng.choice(rest, size=min(len(rest), max_samples - len(keep)), replace=False)
    return df.loc[keep.append(pd.Index(fill))]


def _profile_report_class():
    """ProfileReport from fg-data-profiling (successor) or ydata-profiling, else None."""
    try:
        from data_profiling import ProfileReport
    except ImportError:
        try:
            from ydata_profiling import ProfileReport
        except ImportError:
            return None
    return ProfileReport


def run_eda_profile(group, input_paths, selection_path, eda_config_path, output_path):
    with open(eda_config_path) as f:
        settings = json.load(f)['settings']

    ProfileReport = _profile_report_class()
    if ProfileReport is None:
        print('WARNING: fg-data-profiling / ydata-profiling not installed — skipping '
              f'profile for group {group!r} (metrics and flags are unaffected)',
              file=sys.stderr)
        sys.exit(EXIT_NO_PROFILER)

    df = load_group(group, input_paths, selection_path)
    if df.empty:
        print(f'WARNING: no data for group {group!r}', file=sys.stderr)
        return False
    n_all = len(df)
    df = downsample(df, int(settings['profile_max_samples']), int(settings.get('seed', 42)))

    # minimal=True already turns correlations off; opt back in for 'auto' only.
    extra = {'correlations': {'auto': {'calculate': True}}} \
        if settings.get('profile_correlations') else {}
    report = ProfileReport(
        df.reset_index(drop=True),
        title=f'Phenotype profile: {group} ({df.shape[1]} columns, '
              f'{len(df):,} of {n_all:,} samples)',
        minimal=True,
        interactions=None,
        samples=None,
        duplicates=None,
        progress_bar=False,
        **extra,
    )
    report.to_file(output_path)
    print(f'Wrote: {output_path}', file=sys.stderr)
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--group', required=True)
    p.add_argument('--inputs', nargs='+', required=True)
    p.add_argument('--selection', required=True)
    p.add_argument('--eda_config', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    run_eda_profile(args.group, args.inputs, args.selection, args.eda_config, args.output)


if __name__ == '__main__':
    main()
