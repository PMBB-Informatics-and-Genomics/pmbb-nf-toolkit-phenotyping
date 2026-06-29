import json
import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from basic_categorical import run_basic_categorical


def _long_tsv(tmp_path, phenotype, rows):
    """rows: list of (sample_ID, value)"""
    data = [{'sample_ID': sid, 'data_type': 'categorical', 'phenotype': phenotype,
              'value': str(v), 'occurrence_date': '2020-01-01'} for sid, v in rows]
    df = pd.DataFrame(data)
    p = tmp_path / f'{phenotype}.long.tsv'
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _cfg(tmp_path, phenotype, overrides=None):
    base = {
        'phenotype_name': phenotype,
        'data_type': 'categorical',
        'output_name': phenotype,
        'missing_as_control': False,
        'min_occurrences': 1,
        'one_hot': False,
        'dictionary': None,
        'categories': None,
    }
    if overrides:
        base.update(overrides)
    p = tmp_path / f'{phenotype}.json'
    p.write_text(json.dumps(base))
    return str(p)


def test_binary_encoding(tmp_path):
    tsv = _long_tsv(tmp_path, 'SEX', [('P001', 'female'), ('P002', 'male')])
    cfg = _cfg(tmp_path, 'SEX')
    out = str(tmp_path / 'SEX.cat.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=['P001', 'P002', 'P003'])
    df = pd.read_csv(out, sep='\t')
    assert int(df[df['sample_ID'] == 'P001']['SEX'].iloc[0]) == 1
    assert pd.isna(df[df['sample_ID'] == 'P003']['SEX'].iloc[0])


def test_dictionary_encoding(tmp_path):
    tsv = _long_tsv(tmp_path, 'SEX', [('P001', 'female'), ('P002', 'male'), ('P003', 'unknown')])
    cfg = _cfg(tmp_path, 'SEX', {
        'dictionary': {'female': 0, 'male': 1, 'unknown': 'NA'},
        'output_name': 'SEX',
    })
    out = str(tmp_path / 'SEX.cat.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=['P001', 'P002', 'P003'])
    df = pd.read_csv(out, sep='\t')
    assert int(df[df['sample_ID'] == 'P001']['SEX'].iloc[0]) == 0
    assert int(df[df['sample_ID'] == 'P002']['SEX'].iloc[0]) == 1
    assert pd.isna(df[df['sample_ID'] == 'P003']['SEX'].iloc[0])


def test_one_hot_encoding(tmp_path):
    tsv = _long_tsv(tmp_path, 'ANCESTRY', [('P001', 'AFR'), ('P002', 'EUR'), ('P003', 'AFR')])
    cfg = _cfg(tmp_path, 'ANCESTRY', {
        'one_hot': True,
        'categories': ['AFR', 'EUR', 'AMR'],
        'output_name': 'ANCESTRY',
    })
    out = str(tmp_path / 'ANCESTRY.cat.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=['P001', 'P002', 'P003'])
    df = pd.read_csv(out, sep='\t')
    assert 'ANCESTRY_AFR' in df.columns
    assert 'ANCESTRY_EUR' in df.columns
    assert 'ANCESTRY_AMR' in df.columns
    assert int(df[df['sample_ID'] == 'P001']['ANCESTRY_AFR'].iloc[0]) == 1
    assert int(df[df['sample_ID'] == 'P001']['ANCESTRY_EUR'].iloc[0]) == 0


def test_missing_as_control_false(tmp_path):
    tsv = _long_tsv(tmp_path, 'SEX', [('P001', 'female')])
    cfg = _cfg(tmp_path, 'SEX', {'missing_as_control': False})
    out = str(tmp_path / 'SEX.cat.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    assert pd.isna(df[df['sample_ID'] == 'P002']['SEX'].iloc[0])


def test_missing_as_control_true(tmp_path):
    tsv = _long_tsv(tmp_path, 'SEX', [('P001', 'female')])
    cfg = _cfg(tmp_path, 'SEX', {'missing_as_control': True})
    out = str(tmp_path / 'SEX.cat.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    assert int(df[df['sample_ID'] == 'P002']['SEX'].iloc[0]) == 0


def test_min_occurrences(tmp_path):
    """Sample with only 1 occurrence when min_occurrences=2 → NA."""
    tsv = _long_tsv(tmp_path, 'SEX', [
        ('P001', 'female'),
        ('P002', 'male'), ('P002', 'male'),
    ])
    cfg = _cfg(tmp_path, 'SEX', {'min_occurrences': 2})
    out = str(tmp_path / 'SEX.cat.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    assert pd.isna(df[df['sample_ID'] == 'P001']['SEX'].iloc[0])
    assert not pd.isna(df[df['sample_ID'] == 'P002']['SEX'].iloc[0])
