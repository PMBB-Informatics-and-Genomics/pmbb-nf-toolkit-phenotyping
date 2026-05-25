import json
import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from diagnostic_advanced import run_diagnostic_advanced


def _make_long_tsvs(tmp_path, phenotype_name, rows_by_source):
    """
    rows_by_source: dict of {source_label: [(sample_ID, icd_code, date), ...]}
    Returns list of tsv paths.
    """
    paths = []
    for label, rows in rows_by_source.items():
        data = [
            {'sample_ID': sid, 'concept': 'diagnostic', 'phenotype': phenotype_name,
             'value': code, 'occurrence_date': dt}
            for sid, code, dt in rows
        ]
        p = tmp_path / f'{phenotype_name}_{label}.long.tsv'
        pd.DataFrame(data).to_csv(p, sep='\t', index=False)
        paths.append(str(p))
    return paths


def _make_cfg(tmp_path, sources_override=None, **kwargs):
    """Build a minimal advanced_diagnostic config JSON."""
    sources = sources_override or {
        'ICD10': {
            'from_phenotype': 'ICD10',
            'case_codes': ['E11'],
            'case_exclude': [],
            'control_exclude': [],
        }
    }
    base = {
        'phenotype_name': 'T2Diab',
        'output_name': 'T2Diab',
        'concept': 'advanced_diagnostic',
        'min_occurrences': 1,
        'missing_as_control': False,
        'sources': sources,
    }
    base.update(kwargs)
    p = tmp_path / 'T2Diab.json'
    p.write_text(json.dumps(base))
    return str(p)


def _read_result(path, col='T2Diab'):
    return pd.read_csv(path, sep='\t').set_index('sample_ID')[col]


# ---------------------------------------------------------------------------
# Basic case assignment
# ---------------------------------------------------------------------------

def test_case_gets_1(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001', 'P002'])
    s = _read_result(out)
    assert int(s['P001']) == 1


def test_absent_sample_na_by_default(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, missing_as_control=False)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001', 'P002'])
    s = _read_result(out)
    assert pd.isna(s['P002'])


def test_absent_sample_zero_when_missing_as_control(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, missing_as_control=True)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001', 'P002'])
    s = _read_result(out)
    assert int(s['P002']) == 0


# ---------------------------------------------------------------------------
# case_exclude
# ---------------------------------------------------------------------------

def test_case_exclude_overrides_case_to_0(tmp_path):
    """Patient has E11.9 (case) AND E10.9 (case_exclude) → 0."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [
        ('P001', 'E11.9', '2020-01-01'),
        ('P001', 'E10.9', '2020-01-01'),
    ]})
    cfg = _make_cfg(tmp_path, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11'], 'case_exclude': ['E10']}
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 0


def test_case_without_exclude_still_gets_1(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11'], 'case_exclude': ['E10']}
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 1


# ---------------------------------------------------------------------------
# control_exclude
# ---------------------------------------------------------------------------

def test_control_exclude_sets_na(tmp_path):
    """Patient has no case code but has control_exclude code → NA even with missing_as_control."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E10.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, missing_as_control=True, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11'], 'control_exclude': ['E10']}
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert pd.isna(_read_result(out)['P001'])


def test_control_exclude_does_not_affect_cases(tmp_path):
    """Patient has case code AND control_exclude code — case wins (→ 1)."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [
        ('P001', 'E11.9', '2020-01-01'),
        ('P001', 'E10.9', '2020-01-01'),
    ]})
    cfg = _make_cfg(tmp_path, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11'], 'control_exclude': ['E10']}
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 1


# ---------------------------------------------------------------------------
# Prefix matching
# ---------------------------------------------------------------------------

def test_bare_code_prefix_matches_subcode(tmp_path):
    """'E11' in case_codes matches 'E11.9' in patient data."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 1


def test_explicit_wildcard_still_works(tmp_path):
    """'E11.**' wildcard syntax also matches 'E11.9'."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11.**']}
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 1


def test_prefix_does_not_cross_chapter(tmp_path):
    """'E11' should NOT match 'E110' (no dot separator)."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E110', '2020-01-01')]})
    cfg = _make_cfg(tmp_path)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert pd.isna(_read_result(out)['P001'])


# ---------------------------------------------------------------------------
# min_occurrences
# ---------------------------------------------------------------------------

def test_min_occurrences_not_met_is_not_a_case(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, min_occurrences=2)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert pd.isna(_read_result(out)['P001'])


def test_min_occurrences_met_is_a_case(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [
        ('P001', 'E11.9', '2020-01-01'),
        ('P001', 'E11.9', '2021-03-01'),
    ]})
    cfg = _make_cfg(tmp_path, min_occurrences=2)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 1


def test_min_occurrences_per_code_not_pooled(tmp_path):
    """Two different qualifying codes each once; min_occurrences=2 → not a case
    (each code counted individually, not pooled across codes)."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [
        ('P001', 'E11.9', '2020-01-01'),
        ('P001', 'E11.65', '2021-03-01'),
    ]})
    cfg = _make_cfg(tmp_path, min_occurrences=2)
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert pd.isna(_read_result(out)['P001'])


# ---------------------------------------------------------------------------
# Multi-source
# ---------------------------------------------------------------------------

def test_multi_source_case_from_second_source(tmp_path):
    """Patient has no ICD10 case code but has an ICD9 case code → still a case."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {
        'ICD10': [('P001', 'E10.9', '2020-01-01')],   # ICD10: type 1 (not a case)
        'ICD9':  [('P001', '250.0', '2020-01-01')],   # ICD9: T2D → case
    })
    cfg = _make_cfg(tmp_path, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11']},
        'ICD9':  {'from_phenotype': 'ICD9',  'case_codes': ['250.0', '250.2']},
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 1


def test_multi_source_cross_source_exclude(tmp_path):
    """ICD9 case code overridden by ICD10 case_exclude (global exclude logic)."""
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {
        'ICD9':  [('P001', '250.0', '2020-01-01')],   # ICD9 case
        'ICD10': [('P001', 'E10.9', '2020-01-01')],   # ICD10 exclude code
    })
    cfg = _make_cfg(tmp_path, sources_override={
        'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11'], 'case_exclude': ['E10']},
        'ICD9':  {'from_phenotype': 'ICD9',  'case_codes': ['250.0']},
    })
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    assert int(_read_result(out)['P001']) == 0


# ---------------------------------------------------------------------------
# output_name
# ---------------------------------------------------------------------------

def test_output_column_uses_output_name(tmp_path):
    tsvs = _make_long_tsvs(tmp_path, 'T2Diab', {'ICD10': [('P001', 'E11.9', '2020-01-01')]})
    cfg = _make_cfg(tmp_path, output_name='T2D_clean')
    out = str(tmp_path / 'T2Diab.diag.tsv')
    run_diagnostic_advanced(tsvs, cfg, out, all_samples=['P001'])
    df = pd.read_csv(out, sep='\t')
    assert 'T2D_clean' in df.columns
    assert 'T2Diab' not in df.columns
