import json
import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from basic_coded import run_basic_coded


def _make_long_tsv(tmp_path, code, rows):
    """rows: list of (sample_ID, occurrence_date)"""
    data = [{'sample_ID': sid, 'data_type': 'coded', 'phenotype': code,
              'value': code, 'occurrence_date': dt} for sid, dt in rows]
    df = pd.DataFrame(data)
    p = tmp_path / f'{code}.long.tsv'
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _make_cfg(tmp_path, overrides=None):
    base = {'min_occurrences': 1, 'missing_as_control': False}
    if overrides:
        base.update(overrides)
    p = tmp_path / 'diag.json'
    p.write_text(json.dumps(base))
    return str(p)


def test_present_sample_gets_1(tmp_path):
    tsv = _make_long_tsv(tmp_path, 'N80.0', [('P001', '2020-01-15'), ('P003', '2020-05-01')])
    cfg = _make_cfg(tmp_path)
    out = str(tmp_path / 'N80.0.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001', 'P002', 'P003'])
    df = pd.read_csv(out, sep='\t')
    p001 = df[df['sample_ID'] == 'P001']['N80.0'].iloc[0]
    assert int(p001) == 1


def test_absent_sample_na_by_default(tmp_path):
    tsv = _make_long_tsv(tmp_path, 'N80.0', [('P001', '2020-01-15')])
    cfg = _make_cfg(tmp_path, {'missing_as_control': False})
    out = str(tmp_path / 'N80.0.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    p002_val = df[df['sample_ID'] == 'P002']['N80.0'].iloc[0]
    assert pd.isna(p002_val)


def test_absent_sample_zero_when_missing_as_control(tmp_path):
    tsv = _make_long_tsv(tmp_path, 'N80.0', [('P001', '2020-01-15')])
    cfg = _make_cfg(tmp_path, {'missing_as_control': True})
    out = str(tmp_path / 'N80.0.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    p002_val = df[df['sample_ID'] == 'P002']['N80.0'].iloc[0]
    assert int(p002_val) == 0


def test_min_occurrences_enforced(tmp_path):
    """Sample with only 1 occurrence when min_occurrences=2 → NA."""
    tsv = _make_long_tsv(tmp_path, 'E11.9', [
        ('P001', '2020-01-01'),                          # 1 occurrence
        ('P002', '2020-01-01'), ('P002', '2021-03-01'),  # 2 occurrences
    ])
    cfg = _make_cfg(tmp_path, {'min_occurrences': 2})
    out = str(tmp_path / 'E11.9.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    p001_val = df[df['sample_ID'] == 'P001']['E11.9'].iloc[0]
    p002_val = df[df['sample_ID'] == 'P002']['E11.9'].iloc[0]
    assert pd.isna(p001_val)
    assert int(p002_val) == 1


def test_output_column_named_after_code(tmp_path):
    tsv = _make_long_tsv(tmp_path, 'I11.9', [('P001', '2020-01-01')])
    cfg = _make_cfg(tmp_path)
    out = str(tmp_path / 'I11.9.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001'])
    df = pd.read_csv(out, sep='\t')
    assert 'I11.9' in df.columns
