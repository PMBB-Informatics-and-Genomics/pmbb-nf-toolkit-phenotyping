#!/usr/bin/env python3
"""
advanced_coded.py — Named phenotype from multiple coded long.tsvs.

Reads long.tsvs from one or more basic coded phenotypes, pools all codes
per patient, applies case_codes / case_exclude / control_exclude pattern
matching across all sources, and writes a single binary column.

Usage:
    advanced_coded.py \\
        --input E11.9.long.tsv E11.65.long.tsv 250.0.long.tsv ... \\
        --config T2Diab.json \\
        --output T2Diab.diag.tsv \\
        [--all_samples samples.txt]
"""

import argparse
import json
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from icd_utils import sample_has_any_pattern, build_prefix_index, normalize_patterns

HINT_MAX_CODES = 5


def unmatched_pattern_warnings(raw_patterns, all_codes):
    """
    Return a short warning for each config pattern that matches no code in the input data.

    raw_patterns : list of (field, pattern, prefix_matching) as written in the config
    all_codes    : set of every distinct code in the input long.tsvs (before
                   min_occurrences), so a warning means the code is absent, not rare

    Only the input data is consulted. When codes starting with the pattern's stem
    exist (e.g. 'K74.6' vs 'K74.60'), they are listed as a hint:
        "case_codes pattern 'K74.6' matched no codes (data has K74.60, K74.69)"
    """
    warnings = []
    for field, pattern, prefix_matching in raw_patterns:
        pats = normalize_patterns([pattern], prefix_matching=prefix_matching)
        if not pats or sample_has_any_pattern(all_codes, pats):
            continue
        msg = f"{field} pattern '{pattern}' matched no codes"
        stem = pats[0].rstrip('*').rstrip('.')
        similar = sorted(c for c in all_codes if c.startswith(stem))
        if similar:
            shown = ', '.join(similar[:HINT_MAX_CODES])
            more = f', +{len(similar) - HINT_MAX_CODES} more' if len(similar) > HINT_MAX_CODES else ''
            msg += f' (data has {shown}{more})'
        warnings.append(msg)
    return warnings


def run_advanced_coded(long_tsvs, config_path, output_path, all_samples=None,
                       warnings_output=None):
    """
    Compute a single named phenotype column from multiple coded long.tsvs.

    Parameters
    ----------
    long_tsvs    : list of paths to long.tsv files (one or more basic phenotypes)
    config_path  : path to JSON with sources dict containing case_codes/case_exclude/
                   control_exclude per source, plus min_occurrences, missing_as_control
    output_path  : path to write {phenotype_name}.diag.tsv
    all_samples  : optional list of all sample IDs for roster
    warnings_output : optional path; short warnings are also written here, one per
                      line (empty file when there are none), for the manifest's
                      `warnings` column
    """
    with open(config_path) as f:
        cfg = json.load(f)

    # Pool patterns from all sources (cross-source exclude logic)
    sources = cfg.get('sources') or {}
    pheno_prefix_matching = cfg.get('prefix_matching', True)

    case_codes = []
    case_exclude = []
    control_exclude = []
    raw_patterns = []   # (field, pattern, prefix_matching) for unmatched-pattern warnings
    for s in sources.values():
        src_prefix_matching = s.get('prefix_matching', pheno_prefix_matching)
        case_codes      += normalize_patterns(s.get('case_codes')      or [], prefix_matching=src_prefix_matching)
        case_exclude    += normalize_patterns(s.get('case_exclude')    or [], prefix_matching=src_prefix_matching)
        control_exclude += normalize_patterns(s.get('control_exclude') or [], prefix_matching=src_prefix_matching)
        for field in ('case_codes', 'case_exclude', 'control_exclude'):
            raw_patterns += [(field, str(p).strip(), src_prefix_matching)
                             for p in (s.get(field) or []) if str(p).strip()]

    min_occ         = int(cfg.get('min_occurrences', 1))
    missing_as_ctrl = bool(cfg.get('missing_as_control', False))
    output_name     = cfg.get('output_name') or cfg['phenotype_name']

    # Read and concatenate all input long.tsvs
    frames = [pd.read_csv(p, sep='\t', dtype=str) for p in long_tsvs]
    df = pd.concat(frames, ignore_index=True)

    all_codes = set(df['value'].dropna().astype(str).str.strip())
    warnings = unmatched_pattern_warnings(raw_patterns, all_codes)

    # Build per-sample code sets, respecting min_occurrences per (sample, code).
    counts = df.groupby(['sample_ID', 'value']).size().reset_index(name='n')
    qualifying = counts[counts['n'] >= min_occ]
    sample_codes = (
        qualifying.groupby('sample_ID')['value']
        .apply(set)
        .to_dict()
    )

    # Build prefix indexes once for performance
    case_idx      = build_prefix_index(case_codes)
    excl_idx      = build_prefix_index(case_exclude)
    ctrl_excl_idx = build_prefix_index(control_exclude)

    if all_samples is None:
        all_samples = df['sample_ID'].unique().tolist()

    result = pd.DataFrame({'sample_ID': list(all_samples)})
    values = []
    for sid in all_samples:
        codes = sample_codes.get(sid, set())
        if case_exclude and sample_has_any_pattern(codes, case_exclude, excl_idx):
            values.append(pd.NA)       # excluded: overrides case and control
        elif sample_has_any_pattern(codes, case_codes, case_idx):
            values.append(1)           # confirmed case
        else:
            if control_exclude and sample_has_any_pattern(codes, control_exclude, ctrl_excl_idx):
                values.append(pd.NA)   # excluded from controls
            else:
                values.append(0 if missing_as_ctrl else pd.NA)

    result[output_name] = pd.array(values, dtype='Int64')
    if not (result[output_name] == 1).any():
        warnings.append('0 cases')
    for msg in warnings:
        print(f"WARNING: phenotype '{output_name}': {msg}", file=sys.stderr)
    if warnings_output:
        with open(warnings_output, 'w') as f:
            f.writelines(msg + '\n' for msg in warnings)
    result.to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input',       required=True, nargs='+', help='One or more long-format TSV files')
    p.add_argument('--config',      required=True, help='Per-phenotype JSON config')
    p.add_argument('--output',      required=True, help='Output diag TSV')
    p.add_argument('--all_samples', default=None,
                   help='Optional newline-separated file of all sample IDs')
    p.add_argument('--warnings_output', default=None,
                   help='Optional file to also write warnings to, one per line')
    args = p.parse_args()

    all_samples = None
    if args.all_samples and os.path.isfile(args.all_samples):
        with open(args.all_samples) as f:
            all_samples = [line.strip() for line in f if line.strip()]

    run_advanced_coded(args.input, args.config, args.output, all_samples,
                       args.warnings_output)


if __name__ == '__main__':
    main()
