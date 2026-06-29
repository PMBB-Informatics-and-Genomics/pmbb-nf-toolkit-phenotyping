#!/usr/bin/env python3
import argparse
import fnmatch
import json
import os
import shutil
import sys
from datetime import date

import pandas as pd
import numpy as np

from utils import resolve_sep


def apply_filter(df, filter_dict):
    """
    Apply row-level filters. filter_dict maps column names to lists of glob patterns.
    AND logic across columns, OR logic within a column's pattern list.
    Missing columns are skipped (not an error).
    Returns filtered DataFrame with reset index.
    """
    if not filter_dict:
        return df
    mask = pd.Series([True] * len(df), index=df.index)
    for col, patterns in filter_dict.items():
        if col not in df.columns:
            continue
        col_mask = df[col].astype(str).apply(
            lambda v: any(fnmatch.fnmatch(v, str(p)) for p in patterns)
        )
        mask &= col_mask
    return df[mask].reset_index(drop=True)


def get_matching_codes(series, subsample):
    """
    Return list of unique values in series matching subsample.
    subsample='all' → all unique non-null values.
    subsample=[list] → values matching any pattern (glob).
    """
    unique_codes = series.dropna().unique().tolist()
    if subsample == 'all':
        return unique_codes
    matched = []
    for code in unique_codes:
        if any(fnmatch.fnmatch(str(code), str(p)) for p in subsample):
            matched.append(code)
    return matched


def _read_table(path, sep_config=None):
    sep, engine = resolve_sep(sep_config, path)
    return pd.read_csv(path, sep=sep, engine=engine, dtype=str)


def _resolve_date(cfg, df):
    """Return occurrence_date series: date_col values if present, else reference_date."""
    date_col = cfg.get('date_col')
    ref = cfg.get('reference_date', 'today')
    if ref == 'today':
        ref = date.today().isoformat()
    if date_col and date_col in df.columns:
        return df[date_col].fillna(ref)
    return pd.Series([ref] * len(df), index=df.index)


def build_long(df, cfg, phenotype_name_override=None):
    """
    Build canonical long-format DataFrame.
    Columns: sample_ID, data_type, phenotype, value, occurrence_date
    """
    sample_id_col = cfg['sample_id_col']
    value_col = cfg['value_col']
    phenotype_name = phenotype_name_override or cfg['phenotype_name']
    data_type = cfg['data_type']

    out = pd.DataFrame()
    out['sample_ID'] = df[sample_id_col].values
    out['data_type'] = data_type
    out['phenotype'] = phenotype_name
    out['value'] = df[value_col].values
    out['occurrence_date'] = _resolve_date(cfg, df).values
    return out.dropna(subset=['value']).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Outlier filtering
# ---------------------------------------------------------------------------

def filter_outliers_iqr(series: pd.Series, multiplier: float = 1.5,
                        mode: str = 'cap') -> pd.Series:
    numeric = pd.to_numeric(series, errors='coerce')
    Q1 = numeric.quantile(0.25)
    Q3 = numeric.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - multiplier * IQR
    upper = Q3 + multiplier * IQR
    if mode == 'cap':
        return numeric.clip(lower=lower, upper=upper).astype(str).where(series.notna(), series)
    else:
        return numeric.where((numeric >= lower) & (numeric <= upper)).astype(str)


def filter_outliers_zscore(series: pd.Series, n_sd: float = 3.0,
                           mode: str = 'cap') -> pd.Series:
    numeric = pd.to_numeric(series, errors='coerce')
    mean = numeric.mean()
    std = numeric.std()
    lower = mean - n_sd * std
    upper = mean + n_sd * std
    if mode == 'cap':
        return numeric.clip(lower=lower, upper=upper).astype(str).where(series.notna(), series)
    else:
        return numeric.where((numeric >= lower) & (numeric <= upper)).astype(str)


def apply_outlier_filter(df: pd.DataFrame, method: str, iqr_mult: float,
                         zscore_sd: float, mode: str,
                         phenotypes: list | None) -> tuple[pd.DataFrame, dict]:
    """
    Apply outlier filtering to quantitative rows in the long-format dataframe.
    Returns modified dataframe and a report dict.
    """
    report = {}
    quant_mask = df['data_type'] == 'quantitative'

    if phenotypes:
        pheno_mask = df['phenotype'].isin(phenotypes)
        target_mask = quant_mask & pheno_mask
    else:
        target_mask = quant_mask

    if not target_mask.any():
        return df, report

    target_df = df[target_mask].copy()
    result_parts = []

    for pheno, grp in target_df.groupby('phenotype'):
        original_values = grp['value'].copy()
        numeric = pd.to_numeric(grp['value'], errors='coerce')

        if numeric.isna().all():
            result_parts.append(grp)
            continue

        if method == 'iqr':
            new_values = filter_outliers_iqr(grp['value'], iqr_mult, mode)
        elif method == 'zscore':
            new_values = filter_outliers_zscore(grp['value'], zscore_sd, mode)
        else:
            result_parts.append(grp)
            continue

        changed = (new_values != original_values) & original_values.notna()
        removed = new_values.isna() & original_values.notna()
        report[pheno] = {
            'n_total': len(grp),
            'n_capped': int(changed.sum()),
            'n_removed': int(removed.sum()),
        }
        grp = grp.copy()
        grp['value'] = new_values.astype(str)
        result_parts.append(grp)

    if result_parts:
        filtered = pd.concat(result_parts, ignore_index=True)
        df = df[~target_mask]
        df = pd.concat([df, filtered], ignore_index=True)

    return df, report


# ---------------------------------------------------------------------------
# Scaling
# ---------------------------------------------------------------------------

def apply_scaling(df: pd.DataFrame, method: str,
                  phenotypes: list | None) -> tuple[pd.DataFrame, dict]:
    """Scale quantitative phenotypes."""
    report = {}
    quant_mask = df['data_type'] == 'quantitative'

    if phenotypes:
        pheno_mask = df['phenotype'].isin(phenotypes)
        target_mask = quant_mask & pheno_mask
    else:
        target_mask = quant_mask

    if not target_mask.any():
        return df, report

    target_df = df[target_mask].copy()
    result_parts = []

    for pheno, grp in target_df.groupby('phenotype'):
        numeric = pd.to_numeric(grp['value'], errors='coerce')
        if numeric.isna().all():
            result_parts.append(grp)
            continue

        if method == 'zscore':
            mean, std = numeric.mean(), numeric.std()
            scaled = (numeric - mean) / std if std > 0 else numeric - mean
            report[pheno] = {'method': 'zscore', 'mean': float(mean), 'std': float(std)}
        elif method == 'minmax':
            lo, hi = numeric.min(), numeric.max()
            scaled = (numeric - lo) / (hi - lo) if hi > lo else numeric - lo
            report[pheno] = {'method': 'minmax', 'min': float(lo), 'max': float(hi)}
        else:
            result_parts.append(grp)
            continue

        grp = grp.copy()
        grp['value'] = scaled.astype(str)
        result_parts.append(grp)

    if result_parts:
        scaled_df = pd.concat(result_parts, ignore_index=True)
        df = df[~target_mask]
        df = pd.concat([df, scaled_df], ignore_index=True)

    return df, report


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def _preprocess_coded(cfg, table_path, output_dir, allowed_samples):
    """Stream coded table in chunks, writing per-code TSVs without loading full table."""
    value_col = cfg['value_col']
    sample_id_col = cfg['sample_id_col']
    subsample = cfg.get('subsample', 'all')
    filter_spec = cfg.get('filter') or {}

    target_codes = None if subsample == 'all' else [str(c) for c in subsample]
    use_glob = target_codes is not None and any(
        '*' in c or '?' in c or '[' in c for c in target_codes
    )
    target_set = None if target_codes is None else set(target_codes)

    sep, engine = resolve_sep(cfg.get('sep'), table_path)
    headers_written = set()
    row_counts = {}

    for chunk in pd.read_csv(table_path, sep=sep, engine=engine, dtype=str, chunksize=100_000):
        chunk = apply_filter(chunk, filter_spec)

        if allowed_samples is not None and sample_id_col in chunk.columns:
            chunk = chunk[chunk[sample_id_col].isin(allowed_samples)]

        if target_codes is not None:
            if use_glob:
                chunk = chunk[chunk[value_col].astype(str).apply(
                    lambda v: any(fnmatch.fnmatch(v, p) for p in target_codes)
                )]
            else:
                chunk = chunk[chunk[value_col].astype(str).isin(target_set)]

        if chunk.empty:
            continue

        for code, code_chunk in chunk.groupby(value_col):
            code_str = str(code)
            code_cfg = dict(cfg)
            code_cfg['phenotype_name'] = code_str
            long_rows = build_long(code_chunk, code_cfg)

            safe_code = code_str.replace('/', '_').replace('\\', '_')
            out_path = os.path.join(output_dir, f'{safe_code}.long.tsv')
            write_header = code_str not in headers_written
            long_rows.to_csv(out_path, sep='\t', index=False,
                             mode='w' if write_header else 'a',
                             header=write_header)
            headers_written.add(code_str)
            row_counts[code_str] = row_counts.get(code_str, 0) + len(long_rows)

    if not headers_written:
        print(f'WARNING: no codes matched subsample for {cfg["phenotype_name"]}', file=sys.stderr)
        return

    for code_str, count in sorted(row_counts.items()):
        safe_code = code_str.replace('/', '_').replace('\\', '_')
        print(f'  Wrote: {os.path.join(output_dir, safe_code + ".long.tsv")} ({count} rows)',
              file=sys.stderr)


def run_preprocess(config_path, output_dir, sample_list_path=None):
    """Core logic: read config JSON, process table, write long-format TSV(s)."""
    with open(config_path) as f:
        cfg = json.load(f)

    os.makedirs(output_dir, exist_ok=True)
    phenotype_name = cfg['phenotype_name']

    # preprocessed_path: pass through
    preprocessed_path = cfg.get('preprocessed_path')
    if preprocessed_path:
        out_path = os.path.join(output_dir, f'{phenotype_name}.long.tsv')
        shutil.copy2(preprocessed_path, out_path)
        print(f'Pass-through: {preprocessed_path} → {out_path}', file=sys.stderr)
        return

    allowed_samples = None
    if sample_list_path:
        with open(sample_list_path) as f:
            allowed_samples = {line.strip() for line in f if line.strip()}

    data_type = cfg['data_type']
    table_path = cfg['table']

    # Diagnostic: stream table in chunks to avoid loading full dataset into memory
    if data_type == 'coded':
        _preprocess_coded(cfg, table_path, output_dir, allowed_samples)
        return

    # Quantitative / Categorical: row-per-patient tables, safe to load fully
    df = _read_table(table_path, sep_config=cfg.get('sep'))
    df = apply_filter(df, cfg.get('filter') or {})

    if allowed_samples is not None:
        sample_id_col = cfg['sample_id_col']
        if sample_id_col in df.columns:
            df = df[df[sample_id_col].isin(allowed_samples)].reset_index(drop=True)

    long_df = build_long(df, cfg)

    # Outlier filtering (quantitative only)
    if data_type == 'quantitative' and cfg.get('outlier_method', 'none') != 'none':
        long_df, _ = apply_outlier_filter(
            long_df,
            cfg['outlier_method'],
            cfg.get('outlier_iqr_multiplier', 1.5),
            cfg.get('outlier_zscore_sd', 3.0),
            cfg.get('outlier_mode', 'cap'),
            None,
        )

    # Scaling (quantitative only)
    if data_type == 'quantitative' and cfg.get('scale'):
        long_df, _ = apply_scaling(long_df, cfg.get('scale_method', 'zscore'), None)

    out_path = os.path.join(output_dir, f'{phenotype_name}.long.tsv')
    long_df.to_csv(out_path, sep='\t', index=False)
    print(f'Wrote: {out_path} ({len(long_df)} rows)', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True, help='Per-phenotype JSON config')
    p.add_argument('--output_dir', required=True, help='Directory for output TSV(s)')
    p.add_argument('--sample_list', default=None,
                   help='File with one sample ID per line; rows not in list are excluded')
    args = p.parse_args()
    run_preprocess(args.config, args.output_dir, sample_list_path=args.sample_list)


if __name__ == '__main__':
    main()
