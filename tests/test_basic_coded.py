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


# ── column_prefix ─────────────────────────────────────────────────────────

def test_column_prefix_applied_to_output_column(tmp_path):
    tsv = _make_long_tsv(tmp_path, '401.1', [('P001', '2020-01-15')])
    cfg = _make_cfg(tmp_path, {'column_prefix': 'phe_'})
    out = str(tmp_path / 'out.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001', 'P002'])
    df = pd.read_csv(out, sep='\t')
    assert list(df.columns) == ['sample_ID', 'phe_401.1']
    assert int(df.loc[df['sample_ID'] == 'P001', 'phe_401.1'].iloc[0]) == 1


def test_no_prefix_keeps_raw_code_column(tmp_path):
    tsv = _make_long_tsv(tmp_path, 'N80.0', [('P001', '2020-01-15')])
    cfg = _make_cfg(tmp_path, {'column_prefix': ''})
    out = str(tmp_path / 'out.diag.tsv')
    run_basic_coded(tsv, cfg, out, all_samples=['P001'])
    df = pd.read_csv(out, sep='\t')
    assert list(df.columns) == ['sample_ID', 'N80.0']


def test_output_dir_names_file_after_prefixed_column(tmp_path):
    tsv = _make_long_tsv(tmp_path, '401.1', [('P001', '2020-01-15')])
    cfg = _make_cfg(tmp_path, {'column_prefix': 'icd9_'})
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    run_basic_coded(tsv, cfg, output_dir=str(out_dir), all_samples=['P001'])
    assert (out_dir / 'icd9_401.1.diag.tsv').exists()


def test_prefix_does_not_change_long_tsv_values(tmp_path):
    """advanced_coded reads the raw 'value' column — basic_coded must not rewrite it."""
    tsv = _make_long_tsv(tmp_path, '401.1', [('P001', '2020-01-15')])
    before = open(tsv).read()
    cfg = _make_cfg(tmp_path, {'column_prefix': 'phe_'})
    run_basic_coded(tsv, cfg, str(tmp_path / 'o.diag.tsv'), all_samples=['P001'])
    assert open(tsv).read() == before


def test_cli_all_samples_writes_controls(tmp_path):
    """--all_samples roster (as passed by BASIC_CODED) gives 0s for non-cases."""
    import subprocess
    tsv = _make_long_tsv(tmp_path, 'BI_160.21', [('P001', '2020-01-01')])
    cfg = _make_cfg(tmp_path, {'missing_as_control': True})
    roster = tmp_path / 'samples.txt'
    roster.write_text('P001\nP002\nP003\n')
    out = tmp_path / 'o.diag.tsv'
    script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'basic_coded.py')
    subprocess.run([sys.executable, script, '--input', tsv, '--config', cfg,
                    '--output', str(out), '--all_samples', str(roster)], check=True)
    s = pd.read_csv(out, sep='\t').set_index('sample_ID')['BI_160.21']
    assert s.to_dict() == {'P001': 1, 'P002': 0, 'P003': 0}


def test_without_roster_only_cases_written(tmp_path):
    """No roster → roster is the cases only, so no 0s even with missing_as_control."""
    tsv = _make_long_tsv(tmp_path, 'BI_160.21', [('P001', '2020-01-01')])
    cfg = _make_cfg(tmp_path, {'missing_as_control': True})
    out = str(tmp_path / 'o.diag.tsv')
    run_basic_coded(tsv, cfg, out)
    s = pd.read_csv(out, sep='\t').set_index('sample_ID')['BI_160.21']
    assert s.to_dict() == {'P001': 1}
