import json
import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from apply_mask import run_apply_mask


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_long_tsv(tmp_path, rows, filename='pheno.long.tsv'):
    """rows: list of dicts with keys sample_ID, concept, phenotype, value, occurrence_date"""
    df = pd.DataFrame(rows)
    p = tmp_path / filename
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _write_mask_intervals(tmp_path, rows, filename='mask_intervals.tsv'):
    """rows: list of dicts with keys sample_ID, event_type, window_start, window_end"""
    df = pd.DataFrame(rows)
    p = tmp_path / filename
    df.to_csv(p, sep='\t', index=False)
    return str(p)


def _write_config(tmp_path, event_mask=None, date_col='ENCOUNTER_DATE'):
    cfg = {
        'phenotype_name': 'HbA1c',
        'concept': 'quantitative',
        'date_col': date_col,
        'event_mask': event_mask or ['preg'],
    }
    p = tmp_path / 'HbA1c.json'
    p.write_text(json.dumps(cfg))
    return str(p)


def _make_row(sid, date_str, value='7.5'):
    return {
        'sample_ID': sid, 'concept': 'quantitative', 'phenotype': 'HbA1c',
        'value': value, 'occurrence_date': date_str,
    }


def _read_result(path):
    return pd.read_csv(path, sep='\t', dtype=str)


# ---------------------------------------------------------------------------
# Core filtering
# ---------------------------------------------------------------------------

def test_row_within_window_filtered(tmp_path):
    """Row with date inside the mask window is removed."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P001', '2021-06-15')])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 0


def test_row_outside_window_kept(tmp_path):
    """Row with date outside all mask windows is kept."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P001', '2023-03-01')])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 1


def test_window_boundaries_are_inclusive(tmp_path):
    """Rows on window_start and window_end dates are filtered (inclusive bounds)."""
    long_tsv = _write_long_tsv(tmp_path, [
        _make_row('P001', '2021-01-01'),   # == window_start → filtered
        _make_row('P001', '2021-12-31'),   # == window_end   → filtered
        _make_row('P001', '2020-12-31'),   # one day before  → kept
        _make_row('P001', '2022-01-01'),   # one day after   → kept
    ])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 2
    assert set(df['occurrence_date'].tolist()) == {'2020-12-31', '2022-01-01'}


def test_sample_not_in_mask_unaffected(tmp_path):
    """Samples with no mask windows pass through untouched."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P002', '2021-06-15')])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 1


def test_no_real_dates_copies_input_unchanged(tmp_path):
    """date_col is null → masking skipped, input copied to output as-is with a warning."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P001', '2024-01-01')])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2023-01-01', 'window_end': '2025-01-01',
    }])
    cfg = _write_config(tmp_path, date_col=None)  # no real dates
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    original = _read_result(long_tsv)
    result = _read_result(out)
    pd.testing.assert_frame_equal(result, original)


def test_event_type_not_in_phenotype_mask_ignored(tmp_path):
    """Intervals for event types not listed in phenotype's event_mask are ignored."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P001', '2021-06-15')])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'OTHER_EVENT',  # phenotype only masks 'preg'
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])
    cfg = _write_config(tmp_path, event_mask=['preg'])
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 1


def test_multiple_event_types_pooled_and_merged(tmp_path):
    """Windows from two event types both apply; adjacent windows merge."""
    long_tsv = _write_long_tsv(tmp_path, [
        _make_row('P001', '2021-03-01'),  # in preg window → filtered
        _make_row('P001', '2021-09-01'),  # in loss window → filtered
        _make_row('P001', '2023-01-01'),  # outside all → kept
    ])
    intervals = _write_mask_intervals(tmp_path, [
        {'sample_ID': 'P001', 'event_type': 'preg',
         'window_start': '2021-01-01', 'window_end': '2021-06-30'},
        {'sample_ID': 'P001', 'event_type': 'loss',
         'window_start': '2021-07-01', 'window_end': '2021-12-31'},
    ])
    cfg = _write_config(tmp_path, event_mask=['preg', 'loss'])
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 1
    assert df.iloc[0]['occurrence_date'] == '2023-01-01'


def test_output_schema_matches_input(tmp_path):
    """Output TSV has the same columns as the input long.tsv."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P001', '2023-01-01')])
    intervals = _write_mask_intervals(tmp_path, [])  # empty mask
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    input_cols = list(pd.read_csv(long_tsv, sep='\t', nrows=0).columns)
    output_cols = list(pd.read_csv(out, sep='\t', nrows=0).columns)
    assert input_cols == output_cols


def test_duplicate_windows_in_mask_intervals_handled(tmp_path):
    """Duplicate (sample_ID, event_type, window) rows in mask_intervals are merged correctly."""
    long_tsv = _write_long_tsv(tmp_path, [_make_row('P001', '2021-06-15')])
    # Same window appears twice — should still filter the row exactly once
    intervals = _write_mask_intervals(tmp_path, [
        {'sample_ID': 'P001', 'event_type': 'preg',
         'window_start': '2021-01-01', 'window_end': '2021-12-31'},
        {'sample_ID': 'P001', 'event_type': 'preg',
         'window_start': '2021-01-01', 'window_end': '2021-12-31'},
    ])
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 0


def test_self_masking_still_filters_correctly(tmp_path):
    """
    Solution 1 design: even if a phenotype is both the event source and is masked,
    apply_mask runs normally and filters correctly. No crash, correct output.
    """
    long_tsv = _write_long_tsv(tmp_path, [
        _make_row('P001', '2021-06-15'),   # inside window → filtered
        _make_row('P001', '2023-01-01'),   # outside → kept
    ])
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])
    cfg = _write_config(tmp_path)
    out = str(tmp_path / 'out.long.tsv')
    run_apply_mask(long_tsv, cfg, intervals, out)
    df = _read_result(out)
    assert len(df) == 1
    assert df.iloc[0]['occurrence_date'] == '2023-01-01'


def test_non_default_index_filtered_correctly(tmp_path):
    """_mask_rows must work even if df has a non-default index."""
    from apply_mask import _mask_rows, _build_sample_windows

    rows = [_make_row('P001', '2021-06-15'), _make_row('P001', '2023-01-01')]
    long_tsv = _write_long_tsv(tmp_path, rows)
    intervals = _write_mask_intervals(tmp_path, [{
        'sample_ID': 'P001', 'event_type': 'preg',
        'window_start': '2021-01-01', 'window_end': '2021-12-31',
    }])

    # Read the dataframe and simulate non-default index
    df = pd.read_csv(long_tsv, sep='\t', dtype=str)
    df.index = [10, 20]  # non-default index

    # Read mask intervals
    mask_df = pd.read_csv(intervals, sep='\t', dtype=str)
    sample_windows = _build_sample_windows(mask_df, ['preg'])

    # Call _mask_rows directly
    result = _mask_rows(df, sample_windows)

    # Should have filtered the 2021-06-15 row (inside window) and kept 2023-01-01
    assert len(result) == 1
    assert result.iloc[0]['occurrence_date'] == '2023-01-01'
