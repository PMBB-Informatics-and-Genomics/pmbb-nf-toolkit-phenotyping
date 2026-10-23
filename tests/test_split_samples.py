"""
Tests for bin/split_samples.py
"""

import os
import subprocess
import sys
import tempfile

import pandas as pd
import pytest

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))

from split_samples import split_samples  # noqa: E402

PYTHON = sys.executable
SCRIPT = os.path.join(BIN_DIR, 'split_samples.py')

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_data')
LONG_FILE = os.path.join(TEST_DATA_DIR, 'input_long_format.tsv')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_long_df(n_samples, rows_per_sample=3):
    """Build a minimal long-format DataFrame with deterministic sample IDs."""
    rows = []
    for i in range(1, n_samples + 1):
        sid = f'P{i:04d}'
        for j in range(rows_per_sample):
            rows.append({'sample_ID': sid, 'data_type': 'quantitative',
                         'phenotype': 'BMI', 'value': str(20 + i),
                         'occurrence_date': '2021-01-01'})
    return pd.DataFrame(rows)


def _write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False)


# ---------------------------------------------------------------------------
# TestSplitSamples (library function)
# ---------------------------------------------------------------------------

class TestSplitSamples:
    def test_correct_number_of_chunks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(25), inp)
            files = split_samples(inp, chunk_size=10, output_dir=tmpdir)
            # 25 samples / 10 per chunk → 3 chunks
            assert len(files) == 3

    def test_single_chunk_when_chunk_size_exceeds_samples(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(5), inp)
            files = split_samples(inp, chunk_size=100, output_dir=tmpdir)
            assert len(files) == 1

    def test_all_samples_present_across_chunks(self):
        n = 30
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            df = _make_long_df(n)
            _write_tsv(df, inp)
            files = split_samples(inp, chunk_size=10, output_dir=tmpdir)
            collected = pd.concat(
                [pd.read_csv(f, sep='\t', dtype=str) for f in files],
                ignore_index=True,
            )
            assert set(collected['sample_ID'].unique()) == set(df['sample_ID'].unique())

    def test_no_sample_split_across_chunks(self):
        """Every row for a given sample_ID must be in the same chunk."""
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(20, rows_per_sample=5), inp)
            files = split_samples(inp, chunk_size=7, output_dir=tmpdir)
            for f in files:
                chunk = pd.read_csv(f, sep='\t', dtype=str)
                for sid, grp in chunk.groupby('sample_ID'):
                    # All rows for this sample_ID must be in the same file
                    assert len(grp) == 5, f"Sample {sid} has split rows in {f}"

    def test_total_row_count_preserved(self):
        n, r = 20, 4
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            df = _make_long_df(n, rows_per_sample=r)
            _write_tsv(df, inp)
            files = split_samples(inp, chunk_size=8, output_dir=tmpdir)
            total_rows = sum(
                len(pd.read_csv(f, sep='\t', dtype=str)) for f in files
            )
            assert total_rows == n * r

    def test_output_files_named_correctly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(15), inp)
            files = split_samples(inp, chunk_size=5, output_dir=tmpdir)
            basenames = [os.path.basename(f) for f in files]
            assert basenames == ['chunk_1.tsv', 'chunk_2.tsv', 'chunk_3.tsv']

    def test_zero_padded_names_for_many_chunks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(100), inp)
            files = split_samples(inp, chunk_size=1, output_dir=tmpdir)
            basenames = [os.path.basename(f) for f in files]
            # 100 chunks → 3 digits → chunk_001.tsv … chunk_100.tsv
            assert basenames[0]  == 'chunk_001.tsv'
            assert basenames[-1] == 'chunk_100.tsv'

    def test_chunk_size_of_1(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(3), inp)
            files = split_samples(inp, chunk_size=1, output_dir=tmpdir)
            assert len(files) == 3
            for f in files:
                chunk = pd.read_csv(f, sep='\t', dtype=str)
                assert chunk['sample_ID'].nunique() == 1

    def test_invalid_chunk_size_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(5), inp)
            with pytest.raises(ValueError):
                split_samples(inp, chunk_size=0, output_dir=tmpdir)

    def test_missing_sample_col_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            pd.DataFrame({'id': ['A', 'B'], 'val': [1, 2]}).to_csv(
                inp, sep='\t', index=False
            )
            with pytest.raises(ValueError, match="sample_ID"):
                split_samples(inp, chunk_size=1, output_dir=tmpdir)

    def test_creates_output_dir_if_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(5), inp)
            out_dir = os.path.join(tmpdir, 'new_subdir')
            assert not os.path.exists(out_dir)
            split_samples(inp, chunk_size=3, output_dir=out_dir)
            assert os.path.isdir(out_dir)

    def test_real_test_data(self):
        """Smoke test against the actual test data file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            files = split_samples(LONG_FILE, chunk_size=3, output_dir=tmpdir)
            assert len(files) >= 1
            for f in files:
                df = pd.read_csv(f, sep='\t', dtype=str)
                assert 'sample_ID' in df.columns
                assert df['sample_ID'].notna().all()


# ---------------------------------------------------------------------------
# TestSplitSamplesCLI
# ---------------------------------------------------------------------------

class TestSplitSamplesCLI:
    def test_basic_cli_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, 'input.tsv')
            _write_tsv(_make_long_df(10), inp)
            out_dir = os.path.join(tmpdir, 'chunks')
            proc = subprocess.run(
                [PYTHON, SCRIPT, '--input', inp,
                 '--chunk_size', '4', '--output_dir', out_dir],
                capture_output=True, text=True,
            )
            assert proc.returncode == 0
            chunk_files = [f for f in os.listdir(out_dir) if f.endswith('.tsv')]
            assert len(chunk_files) == 3   # 10 samples / 4 = 3 chunks

    def test_missing_input_fails(self):
        proc = subprocess.run(
            [PYTHON, SCRIPT, '--chunk_size', '10', '--output_dir', '/tmp'],
            capture_output=True, text=True,
        )
        assert proc.returncode != 0

    def test_missing_chunk_size_fails(self):
        proc = subprocess.run(
            [PYTHON, SCRIPT, '--input', LONG_FILE, '--output_dir', '/tmp'],
            capture_output=True, text=True,
        )
        assert proc.returncode != 0
