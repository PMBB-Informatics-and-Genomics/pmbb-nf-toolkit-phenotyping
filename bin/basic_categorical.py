#!/usr/bin/env python3
"""
Basic Categorical module for pmbb-nf-toolkit-phenotyping.

Encodes categorical/demographic phenotypes from long-format data into a
wide-format matrix (one row per sample).

Three encoding modes per phenotype:
  1. Binary (default): 1 if ≥ min_occurrences rows exist for the sample, else 0 or NaN.
  2. Dictionary: map the phenotype's value strings to numeric codes
     (e.g., sex: male→0, female→1, unknown→NA).
  3. One-hot: expand each category into its own 0/1 column
     (e.g., ancestry → ANCESTRY_AFR, ANCESTRY_EUR, …).

Output columns:
  sample_ID | {output_name} | {output_name}_{cat} | …

Usage:
  basic_categorical.py --input long_format.tsv --output categorical_matrix.tsv [options]

See --help for full option list.
"""

import argparse
import sys
import os
import json

import pandas as pd
import numpy as np

from utils import resolve_sep

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
        description='Categorical encoding of long-format biobank data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--input', required=True,
                   help='Long-format TSV/CSV produced by the pre-processing module.')
    p.add_argument('--params_file', default=None,
                   help='YAML/JSON params file with per-phenotype config.')
    p.add_argument('--phenotypes', default=None,
                   help='Comma-separated list (or newline-separated file) of phenotypes to process. '
                        'Null = all categorical phenotypes in params_file or in data.')
    p.add_argument('--missing_as_control', action='store_true', default=False,
                   help='If set, samples absent from a phenotype are coded 0 rather than NaN.')
    p.add_argument('--min_occurrences', type=int, default=None,
                   help='Minimum row count to code a sample as case in binary mode '
                        '(overrides params_file global; default: 1).')
    p.add_argument('--output', required=True,
                   help='Output wide-format TSV file (one row per sample).')
    p.add_argument('--report', default='categorical_report.txt',
                   help='Output plain-text report file.')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Params loading
# ---------------------------------------------------------------------------

GLOBAL_DEFAULTS = {
    'missing_as_control': False,
    'min_occurrences': 1,
    'one_hot': False,
    'dictionary': None,
    'categories': None,
    'output_name': None,
}


def load_params(path):
    """Load YAML or JSON params file. Returns dict with 'global' and 'phenotypes' keys."""
    if not path or not os.path.isfile(path):
        return {'global': {}, 'phenotypes': {}}

    with open(path) as f:
        content = f.read().strip()

    if not content:
        return {'global': {}, 'phenotypes': {}}

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


def resolve_pheno_config(pheno_name, params, cli_missing_override, cli_min_occ_override):
    """
    Merge global defaults → params_file global → params_file phenotype-specific
    → CLI overrides to produce the effective config for one phenotype.
    """
    global_cfg = {**GLOBAL_DEFAULTS, **params.get('global', {})}
    pheno_cfg  = {**global_cfg, **params.get('phenotypes', {}).get(pheno_name, {})}

    if cli_missing_override:
        pheno_cfg['missing_as_control'] = True
    if cli_min_occ_override is not None:
        pheno_cfg['min_occurrences'] = cli_min_occ_override

    # Coerce min_occurrences to int
    try:
        pheno_cfg['min_occurrences'] = int(pheno_cfg.get('min_occurrences', 1))
    except (TypeError, ValueError):
        pheno_cfg['min_occurrences'] = 1
    if pheno_cfg['min_occurrences'] < 1:
        pheno_cfg['min_occurrences'] = 1

    # output_name defaults to the phenotype name
    if not pheno_cfg.get('output_name'):
        pheno_cfg['output_name'] = pheno_name

    # Normalise dictionary: string "NA" values → actual None (will become pd.NA)
    d = pheno_cfg.get('dictionary')
    if isinstance(d, dict):
        pheno_cfg['dictionary'] = {
            str(k): (None if str(v) == 'NA' else v)
            for k, v in d.items()
        }

    return pheno_cfg


# ---------------------------------------------------------------------------
# Encoding helpers
# ---------------------------------------------------------------------------

def _missing_fill(result_df, col, in_data_mask, missing_as_control):
    """Set out-of-data cells to 0 (if missing_as_control) or pd.NA."""
    not_in_data = ~in_data_mask
    if missing_as_control:
        result_df.loc[not_in_data, col] = 0
    else:
        result_df.loc[not_in_data, col] = pd.NA


def encode_binary(subset, cfg, all_samples):
    """
    Binary encoding: 1 if the sample has ≥ min_occurrences rows, else 0.
    Samples entirely absent from the phenotype → NaN (or 0 if missing_as_control).
    """
    output_name      = cfg['output_name']
    min_occ          = cfg['min_occurrences']
    missing_as_ctrl  = cfg['missing_as_control']

    counts = subset.groupby('sample_ID').size()

    result = pd.DataFrame({'sample_ID': all_samples})
    in_data_mask = result['sample_ID'].isin(counts.index)

    result[output_name] = result['sample_ID'].map(counts)
    result.loc[in_data_mask, output_name] = (
        result.loc[in_data_mask, output_name].ge(min_occ).astype('Int64')
    )
    _missing_fill(result, output_name, in_data_mask, missing_as_ctrl)
    result[output_name] = result[output_name].astype('Int64')

    return result[['sample_ID', output_name]]


def encode_dictionary(subset, cfg, all_samples):
    """
    Dictionary encoding: map the value string through the dictionary to a numeric code.
    If a sample has multiple rows, take the mode (most common mapped value).
    Values not in the dictionary → NaN for that row.
    """
    output_name     = cfg['output_name']
    dictionary      = cfg.get('dictionary') or {}
    missing_as_ctrl = cfg['missing_as_control']

    subset = subset.copy()
    subset['_mapped'] = subset['value'].map(
        {str(k): v for k, v in dictionary.items()}
    )

    def _mode(s):
        valid = s.dropna()
        if valid.empty:
            return pd.NA
        return valid.mode().iloc[0]

    per_sample = subset.groupby('sample_ID')['_mapped'].agg(_mode)

    result = pd.DataFrame({'sample_ID': all_samples})
    in_data_mask = result['sample_ID'].isin(per_sample.index)
    result[output_name] = result['sample_ID'].map(per_sample)
    _missing_fill(result, output_name, in_data_mask, missing_as_ctrl)

    return result[['sample_ID', output_name]]


def encode_one_hot(subset, cfg, all_samples):
    """
    One-hot encoding: create one binary column per category.
    Samples that have data for this phenotype but not a specific category → 0.
    Samples with no data at all → NaN (or 0 if missing_as_control).
    """
    output_name     = cfg['output_name']
    missing_as_ctrl = cfg['missing_as_control']
    categories      = cfg.get('categories') or sorted(
        subset['value'].dropna().unique().tolist()
    )

    result       = pd.DataFrame({'sample_ID': all_samples})
    in_data_set  = set(subset['sample_ID'].dropna().unique())
    in_data_mask = result['sample_ID'].isin(in_data_set)

    for cat in categories:
        col     = f'{output_name}_{cat}'
        has_cat = set(subset[subset['value'] == str(cat)]['sample_ID'].unique())
        result[col] = result['sample_ID'].isin(has_cat).astype('Int64')
        _missing_fill(result, col, in_data_mask, missing_as_ctrl)
        result[col] = result[col].astype('Int64')

    return result


# ---------------------------------------------------------------------------
# Per-phenotype dispatcher
# ---------------------------------------------------------------------------

def aggregate_phenotype(long_df, pheno_name, cfg, all_samples):
    """
    Encode one categorical phenotype. Dispatches to binary / dictionary / one-hot
    based on cfg. Returns a DataFrame with sample_ID plus encoded columns.
    """
    subset = long_df[long_df['phenotype'] == pheno_name].copy()

    if cfg.get('dictionary'):
        return encode_dictionary(subset, cfg, all_samples)
    elif cfg.get('one_hot'):
        return encode_one_hot(subset, cfg, all_samples)
    else:
        return encode_binary(subset, cfg, all_samples)


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

CATEGORICAL_DATA_TYPES = {'categorical', 'survey', 'demographic', 'demographics'}


def run_basic_categorical(long_tsv, config_path, output_path, all_samples=None):
    """Encode a single-phenotype long-format TSV to a per-sample categorical matrix."""
    with open(config_path) as f:
        cfg = json.load(f)

    phenotype_name = cfg.get('phenotype_name', 'phenotype')
    output_name = cfg.get('output_name') or phenotype_name
    min_occ = int(cfg.get('min_occurrences', 1))
    missing_as_control = bool(cfg.get('missing_as_control', False))

    # Normalise dictionary NA values
    d = cfg.get('dictionary')
    if isinstance(d, dict):
        cfg['dictionary'] = {str(k): (None if str(v) == 'NA' else v) for k, v in d.items()}
    cfg['output_name'] = output_name
    cfg['min_occurrences'] = min_occ
    cfg['missing_as_control'] = missing_as_control

    sep, engine = resolve_sep(filepath=long_tsv)
    df = pd.read_csv(long_tsv, sep=sep, engine=engine, dtype=str)

    if all_samples is None:
        all_samples = df['sample_ID'].dropna().unique().tolist()

    # Apply min_occurrences: samples with < N rows in the TSV → treated as absent
    counts = df.groupby('sample_ID').size()
    below_threshold = set(counts[counts < min_occ].index)
    df_filtered = df[~df['sample_ID'].isin(below_threshold)]

    result = aggregate_phenotype(df_filtered, phenotype_name, cfg, all_samples)

    # Samples below threshold get NA regardless of missing_as_control
    if below_threshold:
        below_mask = result['sample_ID'].isin(below_threshold)
        stat_cols = [col for col in result.columns if col != 'sample_ID']
        result.loc[below_mask, stat_cols] = pd.NA

    result.to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--config', required=True, help='Per-phenotype JSON config')
    p.add_argument('--output', required=True)
    p.add_argument('--report', default='categorical_report.txt')
    p.add_argument('--all_samples', default=None,
                   help='Optional newline-separated file of all sample IDs')
    args = p.parse_args()

    all_samples = None
    if args.all_samples and os.path.isfile(args.all_samples):
        with open(args.all_samples) as f:
            all_samples = [line.strip() for line in f if line.strip()]

    run_basic_categorical(args.input, args.config, args.output, all_samples)
    with open(args.report, 'w') as f:
        f.write(f'basic_categorical: processed {args.input}\n')


def _legacy_main():
    args = parse_args()

    # ---- Load params ----
    params      = load_params(args.params_file)
    cli_missing = args.missing_as_control
    cli_min_occ = args.min_occurrences

    # ---- Load long-format input ----
    print(f"Reading input: {args.input}", file=sys.stderr)
    sep, engine = resolve_sep(filepath=args.input)
    df  = pd.read_csv(args.input, sep=sep, engine=engine, dtype=str)
    print(f"  Shape: {df.shape}", file=sys.stderr)

    # ---- Filter to categorical rows ----
    if 'data_type' not in df.columns:
        print("WARNING: No 'data_type' column found. Processing all rows as categorical.",
              file=sys.stderr)
        cat_df = df.copy()
    else:
        cat_mask = df['data_type'].str.lower().isin(CATEGORICAL_DATA_TYPES)
        cat_df   = df[cat_mask].copy()

    if cat_df.empty:
        print("WARNING: No categorical rows found in input.", file=sys.stderr)

    # ---- Determine phenotypes to process ----
    cli_phenos = parse_list_or_file(args.phenotypes)
    if cli_phenos:
        pheno_list = cli_phenos
    elif params.get('phenotypes'):
        pheno_list = list(params['phenotypes'].keys())
    else:
        pheno_list = (
            sorted(cat_df['phenotype'].dropna().unique().tolist())
            if not cat_df.empty else []
        )

    # ---- All samples present in input (for complete output roster) ----
    all_samples = df['sample_ID'].dropna().unique()

    # ---- Encode each phenotype ----
    report_lines = [
        'pmbb-nf-toolkit-phenotyping: basic_categorical report',
        f'Input: {args.input}',
        f'Input rows (categorical): {len(cat_df)}',
        f'Phenotypes processed: {len(pheno_list)}',
        '',
    ]

    all_results = [pd.DataFrame({'sample_ID': all_samples})]

    for pheno in pheno_list:
        cfg = resolve_pheno_config(pheno, params, cli_missing, cli_min_occ)
        mode = ('dictionary' if cfg.get('dictionary')
                else 'one_hot' if cfg.get('one_hot') else 'binary')

        pheno_result = aggregate_phenotype(cat_df, pheno, cfg, all_samples)

        stat_cols = [c for c in pheno_result.columns if c != 'sample_ID']
        report_lines.append(
            f'  {pheno} → {cfg["output_name"]}: mode={mode}, '
            f'min_occurrences={cfg["min_occurrences"]}, columns={stat_cols}'
        )

        if not pheno_result.empty:
            all_results.append(pheno_result)

    # ---- Merge all phenotypes on sample_ID (left join to keep all samples) ----
    matrix = all_results[0]
    for part in all_results[1:]:
        matrix = matrix.merge(part, on='sample_ID', how='left')

    # ---- Write output ----
    print(f"Writing output: {args.output} ({len(matrix)} rows × {len(matrix.columns)} cols)",
          file=sys.stderr)
    matrix.to_csv(args.output, sep='\t', index=False)

    report_lines.extend([
        '',
        f'Output: {args.output}',
        f'Output rows (samples): {len(matrix)}',
        f'Output columns: {list(matrix.columns)}',
    ])
    with open(args.report, 'w') as f:
        f.write('\n'.join(report_lines) + '\n')
    print(f"Report written: {args.report}", file=sys.stderr)


if __name__ == '__main__':
    main()
