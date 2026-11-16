import json
import os
import sys
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from risk_classification import apply_rule, run_risk_classification


# ── helpers ───────────────────────────────────────────────────────────────────

def _series(values, name='col'):
    """Build a float Series with sample_ID-like index."""
    return pd.Series(
        [float(v) if v is not None else np.nan for v in values],
        name=name,
    )


def _wide_tsv(tmp_path, data):
    """Write a wide TSV dict→DataFrame and return path."""
    df = pd.DataFrame(data)
    p = tmp_path / 'gathered.tsv'
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _risk_config(tmp_path, rules):
    """Write risk config JSON and return path."""
    p = tmp_path / 'risk_config.json'
    p.write_text(json.dumps({'rules': rules}))
    return str(p)


# ── percentile binary ─────────────────────────────────────────────────────────

def test_percentile_high_only():
    """Values strictly above p90 → 1.0; all others → 0.0."""
    vals = _series(list(range(1, 11)))   # 1..10; p90 ≈ 9.1
    rule = {'method': 'percentile', 'high': 90, 'stat': 'mean',
            'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    # value 10 is the only one > p90
    assert result.iloc[9] == 1.0    # index 9 = value 10
    assert result.iloc[8] == 0.0    # index 8 = value 9
    assert result.iloc[0] == 0.0    # index 0 = value 1


def test_percentile_low_only():
    """Values strictly below p10 → 1.0; all others → 0.0."""
    vals = _series(list(range(1, 11)))   # p10 ≈ 1.9
    rule = {'method': 'percentile', 'low': 10, 'stat': 'mean',
            'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    # value 1 is the only one < p10
    assert result.iloc[0] == 1.0
    assert result.iloc[1] == 0.0
    assert result.iloc[9] == 0.0


def test_percentile_both_tails():
    """Values above p90 OR below p10 → 1.0; middle values → 0.0."""
    vals = _series(list(range(1, 11)))
    rule = {'method': 'percentile', 'high': 90, 'low': 10,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] == 1.0   # value 1 (low tail)
    assert result.iloc[9] == 1.0   # value 10 (high tail)
    assert result.iloc[4] == 0.0   # value 5 (middle)


def test_percentile_multilevel():
    """cutoffs=[33, 67] → 3 levels: 0 (low), 1 (mid), 2 (high)."""
    vals = _series(list(range(1, 11)))
    rule = {'method': 'percentile', 'cutoffs': [33, 67],
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.notna().all()
    # All three levels must be present
    assert 0 in result.values
    assert 1 in result.values
    assert 2 in result.values
    # Highest values must be in level 2
    assert result.iloc[9] == 2    # value 10 → top level
    # Lowest values must be in level 0
    assert result.iloc[0] == 0    # value 1 → bottom level


def test_percentile_multilevel_with_labels():
    """Labels replace integer levels in output."""
    vals = _series(list(range(1, 11)))
    rule = {'method': 'percentile', 'cutoffs': [33, 67],
            'labels': ['low', 'medium', 'high'],
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[9] == 'high'
    assert result.iloc[0] == 'low'
    assert set(result.dropna().unique()).issubset({'low', 'medium', 'high'})


# ── threshold binary ──────────────────────────────────────────────────────────

def test_threshold_high_only():
    """Values strictly above raw threshold → 1.0; rest → 0.0."""
    vals = _series([18.0, 25.0, 30.0, 35.0, 40.0])
    rule = {'method': 'threshold', 'high': 30,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    # 35.0 and 40.0 are > 30
    assert result.iloc[0] == 0.0   # 18.0
    assert result.iloc[2] == 0.0   # exactly 30.0 → not strictly greater
    assert result.iloc[3] == 1.0   # 35.0
    assert result.iloc[4] == 1.0   # 40.0


def test_threshold_low_only():
    """Values strictly below raw threshold → 1.0; rest → 0.0."""
    vals = _series([10.0, 18.0, 18.5, 25.0, 30.0])
    rule = {'method': 'threshold', 'low': 18.5,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] == 1.0   # 10.0 < 18.5
    assert result.iloc[1] == 1.0   # 18.0 < 18.5
    assert result.iloc[2] == 0.0   # exactly 18.5 → not strictly less
    assert result.iloc[3] == 0.0   # 25.0


def test_threshold_both_tails():
    """Values > high OR < low → 1.0; between → 0.0."""
    vals = _series([10.0, 18.5, 25.0, 30.0, 40.0])
    rule = {'method': 'threshold', 'high': 30, 'low': 18.5,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] == 1.0   # 10.0 < 18.5
    assert result.iloc[1] == 0.0   # exactly 18.5 → not < 18.5
    assert result.iloc[2] == 0.0   # 25.0 in range
    assert result.iloc[3] == 0.0   # exactly 30.0 → not > 30
    assert result.iloc[4] == 1.0   # 40.0 > 30


def test_threshold_multilevel():
    """cutoffs=[18.5, 30] (raw values) → 3 levels."""
    vals = _series([10.0, 20.0, 32.0])
    rule = {'method': 'threshold', 'cutoffs': [18.5, 30],
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] == 0    # 10.0 <= 18.5 → level 0
    assert result.iloc[1] == 1    # 20.0 > 18.5 and <= 30 → level 1
    assert result.iloc[2] == 2    # 32.0 > 30 → level 2


# ── quantile_bin binary ───────────────────────────────────────────────────────

def test_quantile_bin_high_tail_binary():
    """Top 1 bin of 4 quartile bins → 1.0; rest → 0.0."""
    # 8 samples with values 1–8; quartile bins: Q1=[1,2], Q2=[3,4], Q3=[5,6], Q4=[7,8]
    vals = _series([1, 2, 3, 4, 5, 6, 7, 8])
    rule = {'method': 'quantile_bin', 'n': 4, 'high': 1,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    # values 7, 8 are in top bin → 1
    assert result.iloc[6] == 1.0   # value 7
    assert result.iloc[7] == 1.0   # value 8
    assert result.iloc[4] == 0.0   # value 5


def test_quantile_bin_low_tail_binary():
    """Bottom 1 bin of 4 quartile bins → 1.0; rest → 0.0."""
    vals = _series([1, 2, 3, 4, 5, 6, 7, 8])
    rule = {'method': 'quantile_bin', 'n': 4, 'low': 1,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] == 1.0   # value 1
    assert result.iloc[1] == 1.0   # value 2
    assert result.iloc[2] == 0.0   # value 3


def test_quantile_bin_both_tails_binary():
    """Bottom 1 AND top 1 bin of 4 → 1.0; middle 2 → 0.0."""
    vals = _series([1, 2, 3, 4, 5, 6, 7, 8])
    rule = {'method': 'quantile_bin', 'n': 4, 'high': 1, 'low': 1,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] == 1.0   # value 1 (bottom bin)
    assert result.iloc[7] == 1.0   # value 8 (top bin)
    assert result.iloc[3] == 0.0   # value 4 (middle)
    assert result.iloc[4] == 0.0   # value 5 (middle)


def test_quantile_bin_multilevel():
    """No high/low → all bin indices (0, 1, 2, 3) are output levels."""
    vals = _series([1, 2, 3, 4, 5, 6, 7, 8])
    rule = {'method': 'quantile_bin', 'n': 4,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, rep = apply_rule(vals, rule)
    assert result.notna().all()
    unique_levels = set(result.values)
    assert unique_levels == {0, 1, 2, 3}


def test_quantile_bin_multilevel_with_labels():
    """Labels replace bin indices in multi-level output."""
    vals = _series([1, 2, 3, 4, 5, 6, 7, 8])
    rule = {'method': 'quantile_bin', 'n': 4,
            'labels': ['Q1', 'Q2', 'Q3', 'Q4'],
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, _ = apply_rule(vals, rule)
    assert result.iloc[0] in ('Q1',)   # lowest bin
    assert result.iloc[7] in ('Q4',)   # highest bin
    assert set(result.dropna().unique()).issubset({'Q1', 'Q2', 'Q3', 'Q4'})


def test_quantile_bin_all_na_returns_na():
    """All-NA series → all-NA result, no crash (Critical #1 fix)."""
    vals = _series([None, None, None])
    rule = {'method': 'quantile_bin', 'n': 4,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, rep = apply_rule(vals, rule)
    assert result.isna().all()
    assert rep['n_na'] == 3
    assert rep['counts'] == {}


# ── NA propagation ────────────────────────────────────────────────────────────

def test_na_excluded_from_cutoff_computation():
    """NA values are excluded from percentile computation but receive NA in output."""
    vals = _series([None, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    rule = {'method': 'percentile', 'high': 90,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, rep = apply_rule(vals, rule)
    assert pd.isna(result.iloc[0])          # NA input → NA output
    assert rep['n_na'] == 1
    assert result.iloc[9] == 1.0            # value 10 → case


def test_na_threshold_preserves_na():
    """NA in source → NA in output, even for threshold method."""
    vals = _series([None, 25.0, 35.0])
    rule = {'method': 'threshold', 'high': 30,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, rep = apply_rule(vals, rule)
    assert pd.isna(result.iloc[0])
    assert result.iloc[1] == 0.0
    assert result.iloc[2] == 1.0
    assert rep['n_na'] == 1


def test_na_quantile_bin_preserves_na():
    """NA in source → NA in quantile_bin output."""
    vals = _series([None, 1, 2, 3, 4, 5, 6, 7, 8])
    rule = {'method': 'quantile_bin', 'n': 4, 'high': 1,
            'stat': 'mean', 'source_col': 'x', 'output_col': 'x_risk', 'phenotype': 'x'}
    result, rep = apply_rule(vals, rule)
    assert pd.isna(result.iloc[0])
    assert rep['n_na'] == 1


# ── run_risk_classification integration ───────────────────────────────────────

def test_run_risk_classification_output_has_risk_column(tmp_path):
    """run_risk_classification writes a TSV with sample_ID and the risk column."""
    data = {'sample_ID': ['P01', 'P02', 'P03', 'P04', 'P05',
                          'P06', 'P07', 'P08', 'P09', 'P10'],
            'HDL_mean': list(range(1, 11))}
    tsv = _wide_tsv(tmp_path, data)
    config = _risk_config(tmp_path, [{
        'phenotype': 'HDL', 'stat': 'mean',
        'source_col': 'HDL_mean', 'output_col': 'HDL_risk',
        'method': 'percentile', 'high': 90,
    }])
    out = str(tmp_path / 'risk_matrix.tsv')
    report = str(tmp_path / 'risk_report.txt')
    run_risk_classification(tsv, config, out, report)
    df = pd.read_csv(out, sep='\t')
    assert 'sample_ID' in df.columns
    assert 'HDL_risk' in df.columns
    assert len(df) == 10


def test_run_risk_classification_missing_source_col_skipped(tmp_path):
    """When source_col is absent from wide TSV, rule is skipped with a warning."""
    data = {'sample_ID': ['P01', 'P02'], 'OTHER_mean': [1.0, 2.0]}
    tsv = _wide_tsv(tmp_path, data)
    config = _risk_config(tmp_path, [{
        'phenotype': 'HDL', 'stat': 'mean',
        'source_col': 'HDL_mean', 'output_col': 'HDL_risk',
        'method': 'percentile', 'high': 90,
    }])
    out = str(tmp_path / 'risk_matrix.tsv')
    report = str(tmp_path / 'risk_report.txt')
    run_risk_classification(tsv, config, out, report)
    df = pd.read_csv(out, sep='\t')
    # risk column must not appear if source was missing
    assert 'HDL_risk' not in df.columns


def test_run_risk_classification_report_written(tmp_path):
    """risk_report.txt is written and non-empty."""
    data = {'sample_ID': ['P01', 'P02', 'P03', 'P04', 'P05',
                          'P06', 'P07', 'P08', 'P09', 'P10'],
            'BMI_mean': list(range(1, 11))}
    tsv = _wide_tsv(tmp_path, data)
    config = _risk_config(tmp_path, [{
        'phenotype': 'BMI', 'stat': 'mean',
        'source_col': 'BMI_mean', 'output_col': 'BMI_risk',
        'method': 'percentile', 'high': 90,
    }])
    out = str(tmp_path / 'risk_matrix.tsv')
    report = str(tmp_path / 'risk_report.txt')
    run_risk_classification(tsv, config, out, report)
    report_text = open(report).read()
    assert 'BMI_risk' in report_text
    assert 'method=percentile' in report_text
