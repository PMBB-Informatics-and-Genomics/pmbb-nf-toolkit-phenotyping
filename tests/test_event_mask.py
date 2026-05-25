import json
import os
import sys
import pytest
import pandas as pd
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from event_mask import run_event_mask


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_long_tsv(tmp_path, filename, rows):
    """rows: list of (sample_ID, code, date_str)"""
    df = pd.DataFrame([
        {'sample_ID': sid, 'concept': 'diagnostic', 'phenotype': code,
         'value': code, 'occurrence_date': dt}
        for sid, code, dt in rows
    ])
    p = tmp_path / filename
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _write_event_config(tmp_path, events: dict) -> str:
    cfg = {'events': events}
    p = tmp_path / 'event_mask.yaml'
    p.write_text(yaml.dump(cfg))
    return str(p)


def _read_intervals(path):
    return pd.read_csv(path, sep='\t', dtype=str)


# ---------------------------------------------------------------------------
# Window computation
# ---------------------------------------------------------------------------

def test_window_start_and_end_computed_correctly(tmp_path):
    """Event on 2022-01-01 with -270d before and +90d after → window 2021-04-06 to 2022-04-01."""
    tsv = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [('P001', 'O83.3', '2022-01-01')])
    cfg = _write_event_config(tmp_path, {
        'third_trimester': {
            'source': 'ICD10',
            'codes': ['O83.3'],
            'offset_before': -270,
            'offset_after': 90,
        }
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    assert len(df) == 1
    assert df.iloc[0]['sample_ID'] == 'P001'
    assert df.iloc[0]['event_type'] == 'third_trimester'
    assert df.iloc[0]['window_start'] == '2021-04-06'
    assert df.iloc[0]['window_end'] == '2022-04-01'


def test_prefix_match_bare_code(tmp_path):
    """Event code 'O83' (bare) prefix-matches long.tsv named 'O83.3.long.tsv'."""
    tsv = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [('P001', 'O83.3', '2022-06-01')])
    cfg = _write_event_config(tmp_path, {
        'delivery': {
            'source': 'ICD10',
            'codes': ['O83'],
            'offset_before': -5,
            'offset_after': 5,
        }
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    assert len(df) == 1
    assert df.iloc[0]['sample_ID'] == 'P001'


def test_exact_code_does_not_match_sibling(tmp_path):
    """Event code 'O83.3' (exact) does NOT match 'O83.4.long.tsv'."""
    tsv = _write_long_tsv(tmp_path, 'O83.4.long.tsv', [('P001', 'O83.4', '2022-06-01')])
    cfg = _write_event_config(tmp_path, {
        'delivery': {
            'source': 'ICD10',
            'codes': ['O83.3'],
            'offset_before': -5,
            'offset_after': 5,
        }
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    assert len(df) == 0


def test_overlapping_windows_merged(tmp_path):
    """Two events for same sample close together → one merged window."""
    tsv = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [
        ('P001', 'O83.3', '2021-01-01'),  # window: 2020-04-06 to 2021-04-01
        ('P001', 'O83.3', '2021-06-01'),  # window: 2020-09-03 to 2021-08-30 — overlaps
    ])
    cfg = _write_event_config(tmp_path, {
        'preg': {'source': 'ICD10', 'codes': ['O83.3'], 'offset_before': -270, 'offset_after': 90}
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    # Both occurrences for P001 should be merged into 1 row
    assert len(df[df['sample_ID'] == 'P001']) == 1


def test_non_overlapping_windows_kept_separate(tmp_path):
    """Two events far apart → two separate windows in output."""
    tsv = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [
        ('P001', 'O83.3', '2019-01-01'),  # window ends ~2019-04-01
        ('P001', 'O83.3', '2022-01-01'),  # window starts ~2021-04-06 — no overlap
    ])
    cfg = _write_event_config(tmp_path, {
        'preg': {'source': 'ICD10', 'codes': ['O83.3'], 'offset_before': -270, 'offset_after': 90}
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    assert len(df[df['sample_ID'] == 'P001']) == 2


def test_wrong_source_skipped(tmp_path):
    """TSVs are from source 'ICD10' but event requests source 'ICD9' → empty output."""
    tsv = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [('P001', 'O83.3', '2022-01-01')])
    cfg = _write_event_config(tmp_path, {
        'delivery': {
            'source': 'ICD9',
            'codes': ['O83'],
            'offset_before': -10,
            'offset_after': 10,
        }
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)  # source_name='ICD10', event wants 'ICD9'
    df = _read_intervals(out)
    assert len(df) == 0


def test_no_matching_codes_outputs_empty_with_headers(tmp_path):
    """No TSVs match any event code → output file with headers and zero data rows."""
    tsv = _write_long_tsv(tmp_path, 'E11.9.long.tsv', [('P001', 'E11.9', '2022-01-01')])
    cfg = _write_event_config(tmp_path, {
        'preg': {'source': 'ICD10', 'codes': ['O83'], 'offset_before': -10, 'offset_after': 10}
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    assert list(df.columns) == ['sample_ID', 'event_type', 'window_start', 'window_end']
    assert len(df) == 0


def test_multiple_samples_independent(tmp_path):
    """Windows for P001 and P002 are computed independently."""
    tsv = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [
        ('P001', 'O83.3', '2021-01-01'),
        ('P002', 'O83.3', '2022-06-01'),
    ])
    cfg = _write_event_config(tmp_path, {
        'preg': {'source': 'ICD10', 'codes': ['O83.3'], 'offset_before': -270, 'offset_after': 90}
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv], out)
    df = _read_intervals(out)
    assert set(df['sample_ID'].tolist()) == {'P001', 'P002'}
    assert len(df) == 2


def test_multiple_events_same_source(tmp_path):
    """Two events from the same source both get windows computed."""
    tsv1 = _write_long_tsv(tmp_path, 'O83.3.long.tsv', [('P001', 'O83.3', '2021-01-01')])
    tsv2 = _write_long_tsv(tmp_path, 'O03.long.tsv',   [('P001', 'O03',   '2020-01-01')])
    cfg = _write_event_config(tmp_path, {
        'third_tri': {
            'source': 'ICD10', 'codes': ['O83.3'], 'offset_before': -270, 'offset_after': 90
        },
        'first_tri_loss': {
            'source': 'ICD10', 'codes': ['O03'], 'offset_before': -90, 'offset_after': 30
        },
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv1, tsv2], out)
    df = _read_intervals(out)
    assert set(df['event_type'].tolist()) == {'third_tri', 'first_tri_loss'}


def test_bare_code_matches_bare_filename_not_sibling(tmp_path):
    """Event code 'O03' (bare) matches 'O03.long.tsv' but NOT 'O04.long.tsv'."""
    tsv_match    = _write_long_tsv(tmp_path, 'O03.long.tsv', [('P001', 'O03',  '2020-01-01')])
    tsv_no_match = _write_long_tsv(tmp_path, 'O04.long.tsv', [('P002', 'O04',  '2020-01-01')])
    cfg = _write_event_config(tmp_path, {
        'loss': {'source': 'ICD10', 'codes': ['O03'], 'offset_before': -90, 'offset_after': 30}
    })
    out = str(tmp_path / 'out.tsv')
    run_event_mask(cfg, 'ICD10', [tsv_match, tsv_no_match], out)
    df = _read_intervals(out)
    assert set(df['sample_ID'].tolist()) == {'P001'}
