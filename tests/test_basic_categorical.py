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
    cfg = _cfg(tmp_path, 'SEX', {'binary': True})
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


# ── auto encoding (default when no dictionary / one_hot / binary) ─────────────

def _run(tmp_path, phenotype, rows, overrides=None, all_samples=None):
    tsv = _long_tsv(tmp_path, phenotype, rows)
    cfg = _cfg(tmp_path, phenotype, overrides)
    out = str(tmp_path / f'{phenotype}.cat.tsv')
    codes = str(tmp_path / f'{phenotype}.codes.tsv')
    run_basic_categorical(tsv, cfg, out, all_samples=all_samples, codes_path=codes)
    df = pd.read_csv(out, sep='\t', dtype=str).set_index('sample_ID')[phenotype]
    return df, codes


def test_auto_is_default_and_passes_numeric_values_through(tmp_path):
    df, codes = _run(tmp_path, 'BATCH', [('P001', '1'), ('P002', '4'), ('P003', '2.5')],
                     all_samples=['P001', 'P002', 'P003', 'P004'])
    assert float(df['P001']) == 1
    assert float(df['P002']) == 4
    assert float(df['P003']) == 2.5
    assert pd.isna(df['P004'])
    assert not os.path.exists(codes)


def test_auto_integer_values_stay_integers(tmp_path):
    df, _ = _run(tmp_path, 'BATCH', [('P001', '1'), ('P002', '4')],
                 all_samples=['P001', 'P002', 'P003'])
    assert df['P001'] == '1'
    assert df['P002'] == '4'


def test_auto_codes_strings_in_natural_order(tmp_path):
    df, codes = _run(tmp_path, 'BATCH', [
        ('P001', 'Freeze 10'), ('P002', 'Freeze 2'), ('P003', 'Freeze 1'),
    ])
    assert df['P003'] == '1'
    assert df['P002'] == '2'
    assert df['P001'] == '3'
    mapping = pd.read_csv(codes, sep='\t', dtype=str)
    assert list(mapping.columns) == ['value', 'code']
    assert mapping.values.tolist() == [['Freeze 1', '1'], ['Freeze 2', '2'], ['Freeze 10', '3']]


def test_auto_mixed_numeric_and_strings_codes_everything(tmp_path):
    df, _ = _run(tmp_path, 'BATCH', [('P001', '5'), ('P002', 'unknown')])
    assert df['P001'] == '1'
    assert df['P002'] == '2'


def test_auto_categories_set_code_order(tmp_path):
    df, codes = _run(tmp_path, 'SEX', [('P001', 'male'), ('P002', 'female'), ('P003', 'other')],
                     {'categories': ['male', 'female']})
    assert df['P001'] == '1'
    assert df['P002'] == '2'
    assert pd.isna(df['P003'])
    mapping = pd.read_csv(codes, sep='\t', dtype=str)
    assert mapping.values.tolist() == [['male', '1'], ['female', '2']]


def test_auto_multiple_rows_takes_mode_then_lowest_code(tmp_path):
    df, _ = _run(tmp_path, 'BATCH', [
        ('P001', 'b'), ('P001', 'c'), ('P001', 'c'),
        ('P002', 'c'), ('P002', 'a'),
    ])
    assert df['P001'] == '3'
    assert df['P002'] == '1'


def test_auto_missing_as_control(tmp_path):
    df, _ = _run(tmp_path, 'BATCH', [('P001', 'x')], {'missing_as_control': True},
                 all_samples=['P001', 'P002'])
    assert df['P001'] == '1'
    assert df['P002'] == '0'


def test_precedence_dictionary_beats_one_hot_and_binary(tmp_path):
    df, codes = _run(tmp_path, 'SEX', [('P001', 'female'), ('P002', 'male')],
                     {'dictionary': {'female': 7, 'male': 8}, 'one_hot': True, 'binary': True})
    assert float(df['P001']) == 7
    assert float(df['P002']) == 8
    assert not os.path.exists(codes)


def test_precedence_one_hot_beats_binary(tmp_path):
    tsv = _long_tsv(tmp_path, 'ANC', [('P001', 'AFR')])
    cfg = _cfg(tmp_path, 'ANC', {'one_hot': True, 'binary': True})
    out = str(tmp_path / 'ANC.cat.tsv')
    run_basic_categorical(tsv, cfg, out)
    assert 'ANC_AFR' in pd.read_csv(out, sep='\t').columns


# ── categorical_encoding (shared precedence rule) ─────────────────────────────

from utils import categorical_encoding


@pytest.mark.parametrize('cfg,expected', [
    ({}, 'auto'),
    ({'dictionary': None, 'one_hot': False, 'binary': False}, 'auto'),
    ({'binary': True}, 'binary'),
    ({'one_hot': True, 'binary': True}, 'one_hot'),
    ({'dictionary': {'a': 1}, 'one_hot': True, 'binary': True}, 'dictionary'),
])
def test_categorical_encoding_precedence(cfg, expected):
    assert categorical_encoding(cfg) == expected
