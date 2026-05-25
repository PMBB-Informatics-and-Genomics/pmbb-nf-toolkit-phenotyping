"""
Tests for bin/merge_wide.py
"""

import os
import subprocess
import sys
import tempfile

import pandas as pd
import pytest

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))

from merge_wide import merge_wide  # noqa: E402

PYTHON = sys.executable
SCRIPT = os.path.join(BIN_DIR, 'merge_wide.py')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False)


def _chunk(n_samples, start_id, cols):
    """Build a minimal wide-format chunk DataFrame."""
    rows = []
    for i in range(start_id, start_id + n_samples):
        row = {'sample_ID': f'P{i:04d}'}
        row.update(cols)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# TestMergeWide (library function)
# ---------------------------------------------------------------------------

class TestMergeWide:
    def test_merges_two_chunks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f1 = os.path.join(tmpdir, 'chunk_001.tsv')
            f2 = os.path.join(tmpdir, 'chunk_002.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            _write_tsv(_chunk(5,  1, {'BMI_mean': '25.0'}), f1)
            _write_tsv(_chunk(5,  6, {'BMI_mean': '27.0'}), f2)
            merged = merge_wide([f1, f2], out)
            assert len(merged) == 10
            assert set(merged['sample_ID']) == {f'P{i:04d}' for i in range(1, 11)}

    def test_all_samples_present(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            files = []
            for batch_start in range(1, 51, 10):
                f = os.path.join(tmpdir, f'chunk_{batch_start:03d}.tsv')
                _write_tsv(_chunk(10, batch_start, {'X': '1'}), f)
                files.append(f)
            out = os.path.join(tmpdir, 'merged.tsv')
            merged = merge_wide(files, out)
            assert len(merged) == 50

    def test_sorted_by_sample_id(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f1 = os.path.join(tmpdir, 'chunk_001.tsv')
            f2 = os.path.join(tmpdir, 'chunk_002.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            # Write in reversed order so merge must sort
            _write_tsv(_chunk(3, 4, {'X': '1'}), f1)
            _write_tsv(_chunk(3, 1, {'X': '1'}), f2)
            merged = merge_wide([f1, f2], out)
            ids = merged['sample_ID'].tolist()
            assert ids == sorted(ids)

    def test_no_sort_preserves_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f1 = os.path.join(tmpdir, 'chunk_001.tsv')
            f2 = os.path.join(tmpdir, 'chunk_002.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            _write_tsv(_chunk(3, 4, {'X': '1'}), f1)
            _write_tsv(_chunk(3, 1, {'X': '1'}), f2)
            merged = merge_wide([f1, f2], out, no_sort=True)
            # First rows should be from f1 (P0004, P0005, P0006)
            assert merged['sample_ID'].iloc[0] == 'P0004'

    def test_missing_column_filled_with_nan(self):
        """A phenotype absent from one chunk must appear as NaN in merged output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            f1 = os.path.join(tmpdir, 'chunk_001.tsv')
            f2 = os.path.join(tmpdir, 'chunk_002.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            # chunk 1 has ALT_mean; chunk 2 does not
            _write_tsv(
                pd.DataFrame([{'sample_ID': 'P0001', 'BMI_mean': '25', 'ALT_mean': '32'}]),
                f1,
            )
            _write_tsv(
                pd.DataFrame([{'sample_ID': 'P0002', 'BMI_mean': '27'}]),
                f2,
            )
            merged = merge_wide([f1, f2], out)
            assert 'ALT_mean' in merged.columns
            p2_alt = merged.loc[merged['sample_ID'] == 'P0002', 'ALT_mean'].iloc[0]
            assert pd.isna(p2_alt)

    def test_single_file_is_passthrough(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, 'chunk_001.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            df = _chunk(5, 1, {'BMI_mean': '25.0'})
            _write_tsv(df, f)
            merged = merge_wide([f], out)
            assert len(merged) == 5
            assert list(merged.columns) == list(df.columns)

    def test_output_written_to_disk(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, 'chunk_001.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            _write_tsv(_chunk(3, 1, {'X': '1'}), f)
            merge_wide([f], out)
            assert os.path.isfile(out)
            disk_df = pd.read_csv(out, sep='\t', dtype=str)
            assert len(disk_df) == 3

    def test_nonexistent_file_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, 'merged.tsv')
            with pytest.raises(FileNotFoundError):
                merge_wide(['/nonexistent/chunk.tsv'], out)

    def test_empty_input_list_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, 'merged.tsv')
            with pytest.raises(FileNotFoundError):
                merge_wide([], out)

    def test_many_chunks_correct_total(self):
        n_chunks, per_chunk = 20, 50
        with tempfile.TemporaryDirectory() as tmpdir:
            files = []
            for i in range(n_chunks):
                f = os.path.join(tmpdir, f'chunk_{i+1:03d}.tsv')
                _write_tsv(_chunk(per_chunk, i * per_chunk + 1, {'X': '1'}), f)
                files.append(f)
            out = os.path.join(tmpdir, 'merged.tsv')
            merged = merge_wide(files, out)
            assert len(merged) == n_chunks * per_chunk

    def test_roundtrip_with_split(self):
        """Split then merge should reproduce the original DataFrame."""
        from split_samples import split_samples

        n = 30
        rows = [{'sample_ID': f'P{i:04d}', 'concept': 'quantitative',
                  'phenotype': 'BMI', 'value': str(20 + i),
                  'occurrence_date': '2021-01-01'}
                for i in range(1, n + 1)]
        original = pd.DataFrame(rows)

        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'original.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            _write_tsv(original, inp)

            chunk_files = split_samples(inp, chunk_size=10, output_dir=tmpdir)
            merge_wide(chunk_files, out)

            merged = pd.read_csv(out, sep='\t', dtype=str)
            assert set(merged['sample_ID']) == set(original['sample_ID'])
            assert len(merged) == len(original)


# ---------------------------------------------------------------------------
# TestMergeWideCLI
# ---------------------------------------------------------------------------

class TestMergeWideCLI:
    def test_basic_cli_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f1 = os.path.join(tmpdir, 'chunk_001.tsv')
            f2 = os.path.join(tmpdir, 'chunk_002.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')
            _write_tsv(_chunk(5, 1, {'BMI_mean': '25'}), f1)
            _write_tsv(_chunk(5, 6, {'BMI_mean': '27'}), f2)
            proc = subprocess.run(
                [PYTHON, SCRIPT, '--inputs', f1, f2, '--output', out],
                capture_output=True, text=True,
            )
            assert proc.returncode == 0
            assert os.path.isfile(out)
            merged = pd.read_csv(out, sep='\t', dtype=str)
            assert len(merged) == 10

    def test_missing_output_arg_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, 'chunk.tsv')
            _write_tsv(_chunk(3, 1, {'X': '1'}), f)
            proc = subprocess.run(
                [PYTHON, SCRIPT, '--inputs', f],
                capture_output=True, text=True,
            )
            assert proc.returncode != 0

    def test_missing_inputs_arg_fails(self):
        proc = subprocess.run(
            [PYTHON, SCRIPT, '--output', '/tmp/out.tsv'],
            capture_output=True, text=True,
        )
        assert proc.returncode != 0
