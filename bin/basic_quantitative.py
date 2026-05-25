#!/usr/bin/env python3
"""
Basic Quantitative module for pmbb-nf-toolkit-phenotyping.

Aggregates quantitative long-format data to per-sample summary statistics
and optional quantile binning.

For each quantitative phenotype:
  1. Subset rows where concept == 'quantitative' (or 'measurements') and
     phenotype == <target>
  2. Per sample: compute requested stats (mean, median, std, min, max, count)
  3. Optionally bin the per-sample stat into quantile bins (qcut_bins)
  4. Combine all phenotypes into a wide-format output (one row per sample)

Output columns:
  sample_ID | {output_name}_mean | {output_name}_median | ... | {output_name}_bin

Usage:
  basic_quantitative.py --input long_format.tsv --output summary.tsv [options]

See --help for full option list.
"""

import argparse
import sys
import os
import json

import pandas as pd
import numpy as np

from utils import resolve_sep

# ---------------------------------------------------------------------------
# YAML is optional — fall back gracefully if pyyaml is absent
# ---------------------------------------------------------------------------
try:
    import yaml
    _HAVE_YAML = True
except ImportError:
    _HAVE_YAML = False


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Per-sample quantitative summary statistics from long-format data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--input', required=True,
                   help='Long-format TSV/CSV produced by the pre-processing module.')
    p.add_argument('--params_file', default=None,
                   help='YAML/JSON params file with per-phenotype config (stats, binning, etc.).')
    p.add_argument('--phenotypes', default=None,
                   help='Comma-separated list (or newline-separated file) of phenotypes to process. '
                        'Overrides params_file phenotype list. Null = all quantitative phenotypes.')
    p.add_argument('--stats', default=None,
                   help='Comma-separated stats to compute globally: mean,median,std,min,max,count,squared. '
                        'Overrides global stats in params_file.')
    p.add_argument('--qcut_bins', type=int, default=None,
                   help='Number of quantile bins (0 or omit = no binning). '
                        'Overrides global qcut_bins in params_file.')
    p.add_argument('--output', required=True,
                   help='Output wide-format TSV file (one row per sample).')
    p.add_argument('--report', default='quantitative_report.txt',
                   help='Output plain-text report file.')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Params loading
# ---------------------------------------------------------------------------

GLOBAL_DEFAULTS = {
    'qcut_bins': 0,
    'stats': ['mean', 'median', 'std'],
    'output_name': None,
    'pregnancy_filter': False,
}

SUPPORTED_STATS = {'mean', 'median', 'std', 'min', 'max', 'count', 'squared'}


def load_params(path):
    """Load YAML or JSON params file. Returns dict with 'global' and 'phenotypes' keys."""
    if not path or not os.path.isfile(path):
        return {'global': {}, 'phenotypes': {}}

    with open(path) as f:
        content = f.read().strip()

    if not content:
        return {'global': {}, 'phenotypes': {}}

    # Try YAML first, fall back to JSON
    if _HAVE_YAML:
        try:
            data = yaml.safe_load(content)
            if isinstance(data, dict):
                return data
        except Exception:
            pass

    try:
        data = json.loads(content)
        if isinstance(data, dict):
            return data
    except Exception:
        pass

    print(f"WARNING: Could not parse params file '{path}' as YAML or JSON. Using defaults.",
          file=sys.stderr)
    return {'global': {}, 'phenotypes': {}}


def resolve_pheno_config(pheno_name, params, cli_stats_override, cli_qcut_override):
    """
    Merge global defaults → params_file global → params_file phenotype-specific
    → CLI overrides to produce the effective config for one phenotype.
    """
    global_cfg = {**GLOBAL_DEFAULTS, **params.get('global', {})}
    pheno_cfg  = {**global_cfg, **params.get('phenotypes', {}).get(pheno_name, {})}

    # CLI overrides take final precedence
    if cli_stats_override is not None:
        pheno_cfg['stats'] = cli_stats_override
    if cli_qcut_override is not None:
        pheno_cfg['qcut_bins'] = cli_qcut_override

    # Normalize stats to list of lowercase strings
    stats = pheno_cfg.get('stats', ['mean', 'median', 'std'])
    if isinstance(stats, str):
        stats = [s.strip() for s in stats.split(',') if s.strip()]
    pheno_cfg['stats'] = [s.lower() for s in stats if s.lower() in SUPPORTED_STATS]

    if not pheno_cfg['stats']:
        pheno_cfg['stats'] = ['mean']

    # output_name defaults to the phenotype name
    if not pheno_cfg.get('output_name'):
        pheno_cfg['output_name'] = pheno_name

    # qcut_bins must be a non-negative integer
    try:
        pheno_cfg['qcut_bins'] = int(pheno_cfg.get('qcut_bins', 0))
    except (TypeError, ValueError):
        pheno_cfg['qcut_bins'] = 0

    return pheno_cfg


# ---------------------------------------------------------------------------
# Per-phenotype aggregation
# ---------------------------------------------------------------------------

STAT_FNS = {
    'mean':    lambda s: s.mean(),
    'median':  lambda s: s.median(),
    'std':     lambda s: s.std(ddof=1),
    'min':     lambda s: s.min(),
    'max':     lambda s: s.max(),
    'count':   lambda s: s.count(),
    'squared': lambda s: s.mean() ** 2,
}


def aggregate_phenotype(long_df, pheno_name, cfg):
    """
    Aggregate one quantitative phenotype to per-sample summary statistics.

    Returns a DataFrame with sample_ID as index and stat columns named
    '{output_name}_{stat}', plus optionally '{output_name}_bin'.
    """
    subset = long_df[long_df['phenotype'] == pheno_name].copy()
    subset['_numeric'] = pd.to_numeric(subset['value'], errors='coerce')
    subset = subset.dropna(subset=['_numeric'])

    output_name = cfg['output_name']
    stats       = cfg['stats']
    qcut_bins   = cfg['qcut_bins']

    if subset.empty:
        return pd.DataFrame(columns=['sample_ID'])

    # Per-sample aggregation
    grouped = subset.groupby('sample_ID')['_numeric']
    result = pd.DataFrame(index=grouped.groups.keys())
    result.index.name = 'sample_ID'

    for stat in stats:
        col = f'{output_name}_{stat}'
        result[col] = grouped.agg(STAT_FNS[stat])

    # Quantile binning on the first stat (usually mean)
    if qcut_bins > 1:
        bin_source_col = f'{output_name}_{stats[0]}'
        bin_col        = f'{output_name}_bin'
        try:
            result[bin_col] = pd.qcut(
                result[bin_source_col],
                q=qcut_bins,
                labels=range(1, qcut_bins + 1),
                duplicates='drop',
            ).astype('Int64')
        except ValueError as e:
            print(f"WARNING: qcut failed for {pheno_name} (bins={qcut_bins}): {e}. "
                  "Skipping binning.", file=sys.stderr)

    return result.reset_index()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_list_or_file(value):
    """Parse a comma-separated string or newline-delimited file into a list."""
    if value is None:
        return None
    value = value.strip()
    if os.path.isfile(value):
        with open(value) as f:
            return [line.strip() for line in f if line.strip()]
    return [v.strip() for v in value.split(',') if v.strip()]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

QUANTITATIVE_CONCEPTS = {'quantitative', 'measurements', 'vitals'}


def run_basic_quantitative(long_tsv, config_path, output_path):
    """
    Aggregate a single-phenotype long-format TSV to per-sample stats.
    Reads all params from the per-phenotype JSON config.
    """
    with open(config_path) as f:
        cfg = json.load(f)

    output_name = cfg.get('output_name') or cfg.get('phenotype_name', 'phenotype')
    stats = cfg.get('stats', ['mean', 'median', 'std'])
    if isinstance(stats, str):
        stats = [s.strip() for s in stats.split(',')]
    stats = [s.lower() for s in stats if s.lower() in SUPPORTED_STATS] or ['mean']
    qcut_bins = int(cfg.get('qcut_bins', 0))
    min_occ = int(cfg.get('min_occurrences', 1))

    sep, engine = resolve_sep(filepath=long_tsv)
    df = pd.read_csv(long_tsv, sep=sep, engine=engine, dtype=str)

    df['_numeric'] = pd.to_numeric(df['value'], errors='coerce')
    df = df.dropna(subset=['_numeric'])

    if df.empty:
        pd.DataFrame(columns=['sample_ID']).to_csv(output_path, sep='\t', index=False)
        return

    grouped = df.groupby('sample_ID')
    counts = grouped.size()
    result = pd.DataFrame({'sample_ID': df['sample_ID'].unique()})

    for stat in stats:
        col = f'{output_name}_{stat}'
        per_sample = grouped['_numeric'].agg(STAT_FNS[stat]).rename(col)
        result = result.merge(per_sample.reset_index(), on='sample_ID', how='left')

    # Apply min_occurrences: set all stat cols to NA for under-threshold samples
    stat_cols = [f'{output_name}_{s}' for s in stats]
    low_occ = result['sample_ID'].map(counts).fillna(0) < min_occ
    result.loc[low_occ, stat_cols] = pd.NA

    # Quantile binning on first stat
    if qcut_bins > 1:
        bin_src = f'{output_name}_{stats[0]}'
        try:
            result[f'{output_name}_bin'] = pd.qcut(
                pd.to_numeric(result[bin_src], errors='coerce'),
                q=qcut_bins,
                labels=range(1, qcut_bins + 1),
                duplicates='drop',
            ).astype('Int64')
        except ValueError as e:
            print(f'WARNING: qcut failed ({e})', file=sys.stderr)

    result.to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description='Per-sample quantitative stats from single-phenotype long-format TSV')
    p.add_argument('--input', required=True)
    p.add_argument('--config', required=True, help='Per-phenotype JSON config')
    p.add_argument('--output', required=True)
    p.add_argument('--report', default='quantitative_report.txt')
    args = p.parse_args()
    run_basic_quantitative(args.input, args.config, args.output)
    with open(args.report, 'w') as f:
        f.write(f'basic_quantitative: processed {args.input}\n')


def _legacy_main():
    args = parse_args()

    # ---- Load params ----
    params = load_params(args.params_file)
    cli_stats  = parse_list_or_file(args.stats)
    cli_qcut   = args.qcut_bins

    # ---- Load long-format input ----
    print(f"Reading input: {args.input}", file=sys.stderr)
    sep, engine = resolve_sep(filepath=args.input)
    df = pd.read_csv(args.input, sep=sep, engine=engine, dtype=str)
    print(f"  Shape: {df.shape}", file=sys.stderr)

    # ---- Filter to quantitative rows ----
    if 'concept' not in df.columns:
        print("WARNING: No 'concept' column found. Processing all rows as quantitative.",
              file=sys.stderr)
        quant_df = df.copy()
    else:
        quant_mask = df['concept'].str.lower().isin(QUANTITATIVE_CONCEPTS)
        quant_df   = df[quant_mask].copy()

    if quant_df.empty:
        print("WARNING: No quantitative rows found in input.", file=sys.stderr)

    # ---- Determine phenotypes to process ----
    cli_phenos = parse_list_or_file(args.phenotypes)
    if cli_phenos:
        pheno_list = cli_phenos
    elif params.get('phenotypes'):
        pheno_list = list(params['phenotypes'].keys())
    else:
        pheno_list = sorted(quant_df['phenotype'].dropna().unique().tolist()) if not quant_df.empty else []

    # ---- All samples present in input (for a complete output roster) ----
    all_samples = df['sample_ID'].dropna().unique()

    # ---- Aggregate each phenotype ----
    report_lines = [
        'pmbb-nf-toolkit-phenotyping: basic_quantitative report',
        f'Input: {args.input}',
        f'Input rows (quantitative): {len(quant_df)}',
        f'Phenotypes processed: {len(pheno_list)}',
        '',
    ]

    all_results = [pd.DataFrame({'sample_ID': all_samples})]

    for pheno in pheno_list:
        cfg = resolve_pheno_config(pheno, params, cli_stats, cli_qcut)

        if cfg.get('pregnancy_filter'):
            print(f"NOTE: {pheno} has pregnancy_filter=true. Pregnancy filtering is applied "
                  "by the PREGNANCY_FILTER module upstream; this module processes data as-is.",
                  file=sys.stderr)

        pheno_result = aggregate_phenotype(quant_df, pheno, cfg)
        n_samples = len(pheno_result)

        stat_cols = [c for c in pheno_result.columns if c != 'sample_ID']
        report_lines.append(
            f'  {pheno} → {cfg["output_name"]}: {n_samples} samples with data, '
            f'columns={stat_cols}, bins={cfg["qcut_bins"]}'
        )

        if not pheno_result.empty:
            all_results.append(pheno_result)

    # ---- Merge all phenotypes on sample_ID (outer join to keep all samples) ----
    summary = all_results[0]
    for part in all_results[1:]:
        summary = summary.merge(part, on='sample_ID', how='left')

    # ---- Write output ----
    print(f"Writing output: {args.output} ({len(summary)} rows × {len(summary.columns)} cols)",
          file=sys.stderr)
    summary.to_csv(args.output, sep='\t', index=False)

    report_lines.extend([
        '',
        f'Output: {args.output}',
        f'Output rows (samples): {len(summary)}',
        f'Output columns: {list(summary.columns)}',
    ])
    with open(args.report, 'w') as f:
        f.write('\n'.join(report_lines) + '\n')
    print(f"Report written: {args.report}", file=sys.stderr)


if __name__ == '__main__':
    main()
