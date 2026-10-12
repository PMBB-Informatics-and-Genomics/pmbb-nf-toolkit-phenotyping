#!/usr/bin/env python3
"""
risk_classification.py — Apply data-driven risk classification to a gathered wide TSV.

Usage:
    risk_classification.py \\
        --wide_tsv gathered.tsv \\
        --risk_config risk_classification_config.json \\
        --output risk_matrix.tsv \\
        --report risk_report.txt
"""

import argparse
import json
import sys
import numpy as np
import pandas as pd


def _classify_binary(series, high_cut=None, low_cut=None):
    """
    Binary classification on a pandas Series.
    Returns float Series: 1.0 if value > high_cut OR < low_cut, else 0.0. NA → NA.
    """
    not_na = series.notna()
    result = pd.Series(np.nan, index=series.index, dtype=float)
    result[not_na] = 0.0
    if high_cut is not None:
        result[not_na & (series > high_cut)] = 1.0
    if low_cut is not None:
        result[not_na & (series < low_cut)] = 1.0
    return result


def _classify_multilevel(series, cutoff_values, labels=None):
    """
    Multi-level classification.
    N cutoffs → N+1 levels. Level 0: value <= c1.
    Level k: c(k) < value <= c(k+1). Level N: value > cN.
    NA → NA. Labels (list of N+1 strings) replace integer levels if provided.
    """
    not_na = series.notna()
    result = pd.Series(np.nan, index=series.index, dtype=object)
    edges = [-np.inf] + list(cutoff_values) + [np.inf]
    for level, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = not_na & (series > lo) & (series <= hi)
        result[mask] = labels[level] if labels else level
    return result


def _level_counts(series):
    """Dict of value → count for non-NA entries."""
    return {str(k): int(v) for k, v in series.dropna().value_counts().items()}


def apply_rule(series, rule):
    """
    Apply a single risk classification rule to a pandas Series.

    Parameters
    ----------
    series : pd.Series of numeric values (may contain NaN)
    rule   : dict with keys: method, and method-specific keys
             (high, low, cutoffs, n, labels, stat, source_col, output_col, phenotype)

    Returns
    -------
    result : pd.Series aligned with series.index
    report : dict with computed cutoffs and class counts
    """
    method = rule['method']
    vals = series.dropna()
    report = {'method': method, 'stat': rule.get('stat', 'mean'),
              'n_na': int(series.isna().sum())}

    if vals.empty:
        result = pd.Series(np.nan, index=series.index, dtype=float)
        report['counts'] = {}
        return result, report

    if method == 'percentile':
        cutoffs_pct = rule.get('cutoffs')
        high = rule.get('high')
        low = rule.get('low')
        labels = rule.get('labels')

        if cutoffs_pct is not None:
            cutoff_values = [float(np.nanpercentile(vals, c)) for c in cutoffs_pct]
            result = _classify_multilevel(series, cutoff_values, labels)
            report['cutoff_values'] = {str(pct): cv for pct, cv in zip(cutoffs_pct, cutoff_values)}
        else:
            high_cut = float(np.nanpercentile(vals, high)) if high is not None else None
            low_cut = float(np.nanpercentile(vals, low)) if low is not None else None
            result = _classify_binary(series, high_cut, low_cut)
            if high_cut is not None:
                report['high_cutoff'] = high_cut
            if low_cut is not None:
                report['low_cutoff'] = low_cut

    elif method == 'threshold':
        cutoffs_raw = rule.get('cutoffs')
        high = rule.get('high')
        low = rule.get('low')
        labels = rule.get('labels')

        if cutoffs_raw is not None:
            result = _classify_multilevel(series, cutoffs_raw, labels)
            report['cutoff_values'] = cutoffs_raw
        else:
            result = _classify_binary(series, high, low)
            if high is not None:
                report['high_cutoff'] = high
            if low is not None:
                report['low_cutoff'] = low

    elif method == 'quantile_bin':
        n = rule['n']
        high = rule.get('high')
        low = rule.get('low')
        labels = rule.get('labels')

        bins, bin_edges = pd.qcut(vals, q=n, labels=False, retbins=True, duplicates='drop')
        bin_series = pd.Series(np.nan, index=series.index, dtype=float)
        bin_series[vals.index] = bins.values.astype(float)
        actual_n = int(bins.max()) + 1
        report['n_bins'] = actual_n
        report['bin_edges'] = [round(float(e), 4) for e in bin_edges]

        if high is None and low is None:
            result = pd.Series(np.nan, index=series.index, dtype=object)
            not_na = bin_series.notna()
            if labels:
                result[not_na] = bin_series[not_na].apply(lambda x: labels[int(x)])
            else:
                result[not_na] = bin_series[not_na].astype(int)
        else:
            result = pd.Series(np.nan, index=series.index, dtype=float)
            not_na = bin_series.notna()
            result[not_na] = 0.0
            if high is not None:
                result[not_na & (bin_series >= actual_n - high)] = 1.0
            if low is not None:
                result[not_na & (bin_series < low)] = 1.0

    else:
        raise ValueError(f"Unknown risk classification method: {method!r}. "
                         f"Expected one of: percentile, threshold, quantile_bin.")

    report['counts'] = _level_counts(result)
    return result, report


def run_risk_classification(wide_tsv, risk_config_path, output_path, report_path):
    """
    Apply all rules from risk_config to wide_tsv; write risk_matrix.tsv and report.
    """
    wide_df = pd.read_csv(wide_tsv, sep='\t')
    with open(risk_config_path) as f:
        risk_config = json.load(f)

    rules = risk_config.get('rules', [])
    if not rules:
        print('WARNING: no risk classification rules found in config', file=sys.stderr)

    result_df = wide_df[['sample_ID']].copy()
    report_lines = []

    for rule in rules:
        source_col = rule['source_col']
        output_col = rule['output_col']

        if source_col not in wide_df.columns:
            print(
                f"WARNING: source column {source_col!r} not found in wide TSV "
                f"— skipping rule for phenotype '{rule['phenotype']}'",
                file=sys.stderr,
            )
            continue

        series = pd.to_numeric(wide_df[source_col], errors='coerce')
        result, rep = apply_rule(series, rule)
        result_df[output_col] = result.values

        # Build report block
        report_lines.append(
            f"{output_col}  method={rep['method']}  stat={rep['stat']}"
        )
        for key in ('high_cutoff', 'low_cutoff', 'cutoff_values', 'n_bins', 'bin_edges'):
            if key in rep:
                report_lines.append(f'  {key}={rep[key]}')
        report_lines.append(f"  n_na={rep['n_na']}")
        for cls, cnt in sorted(rep.get('counts', {}).items(), key=lambda x: str(x[0])):
            report_lines.append(f'  class={cls}: n={cnt}')
        report_lines.append('')

    result_df.to_csv(output_path, sep='\t', index=False)
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    print(
        f'Wrote {len(result_df.columns) - 1} risk column(s) to {output_path}',
        file=sys.stderr,
    )


def main():
    p = argparse.ArgumentParser(
        description='Apply risk classification rules to a gathered wide TSV'
    )
    p.add_argument('--wide_tsv',     required=True, help='Path to gathered wide TSV')
    p.add_argument('--risk_config',  required=True, help='Path to risk_classification_config.json')
    p.add_argument('--output',       required=True, help='Output path for risk_matrix.tsv')
    p.add_argument('--report',       required=True, help='Output path for risk_report.txt')
    args = p.parse_args()
    run_risk_classification(args.wide_tsv, args.risk_config, args.output, args.report)


if __name__ == '__main__':
    main()
