import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from gather import run_gather


def _write_tsv(tmp_path, name, data):
    """data: dict of {col: [values]}"""
    df = pd.DataFrame(data)
    p = tmp_path / name
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def test_wide_format_join(tmp_path):
    bmi = _write_tsv(tmp_path, 'BMI.quant.tsv', {
        'sample_ID': ['P001', 'P002'],
        'BMI_mean': [28.4, 31.7],
    })
    sex = _write_tsv(tmp_path, 'SEX.cat.tsv', {
        'sample_ID': ['P001', 'P003'],
        'SEX': [0, 1],
    })
    out_prefix = str(tmp_path / 'out')
    run_gather([bmi, sex], out_prefix, fmt='wide')
    df = pd.read_csv(out_prefix + '.wide.tsv', sep='\t')
    assert set(df.columns) == {'sample_ID', 'BMI_mean', 'SEX'}
    assert len(df) == 3   # P001, P002, P003 — outer join


def test_wide_missing_filled_with_na(tmp_path):
    bmi = _write_tsv(tmp_path, 'BMI.quant.tsv', {
        'sample_ID': ['P001', 'P002'],
        'BMI_mean': [28.4, 31.7],
    })
    sex = _write_tsv(tmp_path, 'SEX.cat.tsv', {
        'sample_ID': ['P001'],
        'SEX': [0],
    })
    out_prefix = str(tmp_path / 'out')
    run_gather([bmi, sex], out_prefix, fmt='wide')
    df = pd.read_csv(out_prefix + '.wide.tsv', sep='\t')
    p002_sex = df[df['sample_ID'] == 'P002']['SEX'].iloc[0]
    assert pd.isna(p002_sex)


def test_long_format_output(tmp_path):
    bmi = _write_tsv(tmp_path, 'BMI.quant.tsv', {
        'sample_ID': ['P001', 'P002'],
        'BMI_mean': [28.4, 31.7],
    })
    out_prefix = str(tmp_path / 'out')
    run_gather([bmi], out_prefix, fmt='long')
    df = pd.read_csv(out_prefix + '.long.tsv', sep='\t')
    assert set(df.columns) == {'sample_ID', 'phenotype', 'value'}
    assert set(df['phenotype'].unique()) == {'BMI_mean'}
    assert len(df) == 2


def test_both_formats(tmp_path):
    bmi = _write_tsv(tmp_path, 'BMI.quant.tsv', {
        'sample_ID': ['P001'],
        'BMI_mean': [28.4],
    })
    out_prefix = str(tmp_path / 'out')
    run_gather([bmi], out_prefix, fmt='both')
    assert os.path.exists(out_prefix + '.wide.tsv')
    assert os.path.exists(out_prefix + '.long.tsv')


def test_long_format_multi_file_no_cross_nans(tmp_path):
    """Long format should only emit rows where data exists — no NaN cross-product rows."""
    bmi = _write_tsv(tmp_path, 'BMI.quant.tsv', {
        'sample_ID': ['P001', 'P002'],
        'BMI_mean': [28.4, 31.7],
    })
    sex = _write_tsv(tmp_path, 'SEX.cat.tsv', {
        'sample_ID': ['P001', 'P003'],
        'SEX': [0, 1],
    })
    out_prefix = str(tmp_path / 'out')
    run_gather([bmi, sex], out_prefix, fmt='long')
    df = pd.read_csv(out_prefix + '.long.tsv', sep='\t')
    # P001 has both phenotypes, P002 only BMI_mean, P003 only SEX → 4 rows
    assert len(df) == 4
    assert not df['value'].isna().any(), 'long format should not contain NaN values'


def test_single_file(tmp_path):
    bmi = _write_tsv(tmp_path, 'BMI.quant.tsv', {
        'sample_ID': ['P001'],
        'BMI_mean': [28.4],
        'BMI_std': [1.2],
    })
    out_prefix = str(tmp_path / 'out')
    run_gather([bmi], out_prefix, fmt='wide')
    df = pd.read_csv(out_prefix + '.wide.tsv', sep='\t')
    assert 'BMI_mean' in df.columns
    assert 'BMI_std' in df.columns
