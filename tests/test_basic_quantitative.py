import json
import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from basic_quantitative import run_basic_quantitative


def _long_tsv(tmp_path, rows):
    """rows: list of (sample_ID, value)"""
    data = [{'sample_ID': sid, 'data_type': 'quantitative', 'phenotype': 'BMI',
              'value': str(v), 'occurrence_date': '2020-01-01'} for sid, v in rows]
    df = pd.DataFrame(data)
    p = tmp_path / 'BMI.long.tsv'
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _cfg(tmp_path, overrides=None):
    base = {
        'phenotype_name': 'BMI',
        'data_type': 'quantitative',
        'output_name': 'BMI',
        'stats': ['mean', 'std'],
        'qcut_bins': 0,
        'min_occurrences': 1,
    }
    if overrides:
        base.update(overrides)
    p = tmp_path / 'BMI.json'
    p.write_text(json.dumps(base))
    return str(p)


def test_mean_and_std_computed(tmp_path):
    tsv = _long_tsv(tmp_path, [('P001', 10.0), ('P001', 20.0), ('P002', 15.0)])
    cfg = _cfg(tmp_path)
    out = str(tmp_path / 'BMI.quant.tsv')
    run_basic_quantitative(tsv, cfg, out)
    df = pd.read_csv(out, sep='\t')
    p001 = df[df['sample_ID'] == 'P001'].iloc[0]
    assert abs(float(p001['BMI_mean']) - 15.0) < 0.01
    assert abs(float(p001['BMI_std']) - 7.07) < 0.1


def test_min_occurrences_sets_na(tmp_path):
    """Sample with < min_occurrences rows → NA for all stat columns."""
    tsv = _long_tsv(tmp_path, [('P001', 10.0), ('P002', 20.0), ('P002', 22.0)])
    cfg = _cfg(tmp_path, {'min_occurrences': 2})
    out = str(tmp_path / 'BMI.quant.tsv')
    run_basic_quantitative(tsv, cfg, out)
    df = pd.read_csv(out, sep='\t')
    p001_mean = df[df['sample_ID'] == 'P001']['BMI_mean'].iloc[0]
    p002_mean = df[df['sample_ID'] == 'P002']['BMI_mean'].iloc[0]
    assert pd.isna(p001_mean)
    assert not pd.isna(p002_mean)


def test_qcut_bins(tmp_path):
    rows = [(f'P{i:03d}', float(i * 10)) for i in range(1, 9)]
    tsv = _long_tsv(tmp_path, rows)
    cfg = _cfg(tmp_path, {'stats': ['mean'], 'qcut_bins': 4})
    out = str(tmp_path / 'BMI.quant.tsv')
    run_basic_quantitative(tsv, cfg, out)
    df = pd.read_csv(out, sep='\t')
    assert 'BMI_bin' in df.columns
    assert df['BMI_bin'].notna().any()


def test_output_name_override(tmp_path):
    tsv = _long_tsv(tmp_path, [('P001', 28.4)])
    cfg = _cfg(tmp_path, {'output_name': 'BODY_MASS'})
    out = str(tmp_path / 'BMI.quant.tsv')
    run_basic_quantitative(tsv, cfg, out)
    df = pd.read_csv(out, sep='\t')
    assert 'BODY_MASS_mean' in df.columns


def test_single_stat(tmp_path):
    tsv = _long_tsv(tmp_path, [('P001', 28.4), ('P001', 29.0)])
    cfg = _cfg(tmp_path, {'stats': ['median']})
    out = str(tmp_path / 'BMI.quant.tsv')
    run_basic_quantitative(tsv, cfg, out)
    df = pd.read_csv(out, sep='\t')
    assert 'BMI_median' in df.columns
    assert 'BMI_mean' not in df.columns
