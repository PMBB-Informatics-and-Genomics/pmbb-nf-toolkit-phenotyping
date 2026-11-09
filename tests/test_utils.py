import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from utils import resolve_sep


# ── Named aliases ─────────────────────────────────────────────────────────────

def test_tsv_alias():
    assert resolve_sep('tsv') == ('\t', 'c')

def test_tab_alias():
    assert resolve_sep('tab') == ('\t', 'c')

def test_csv_alias():
    assert resolve_sep('csv') == (',', 'c')

def test_comma_alias():
    assert resolve_sep('comma') == (',', 'c')

def test_pipe_alias():
    assert resolve_sep('pipe') == ('|', 'c')

def test_space_alias():
    assert resolve_sep('space') == (' ', 'c')

def test_whitespace_alias():
    assert resolve_sep('whitespace') == (r'\s+', 'python')


# ── Raw / escape edge cases ───────────────────────────────────────────────────

def test_raw_backslash_t_string():
    """Unquoted \\t in YAML arrives as the two-char string r'\\t'."""
    assert resolve_sep(r'\t') == ('\t', 'c')

def test_raw_backslash_s_plus_string():
    assert resolve_sep(r'\s+') == (r'\s+', 'python')

def test_raw_single_char_pipe():
    assert resolve_sep('|') == ('|', 'c')

def test_raw_single_char_space():
    assert resolve_sep(' ') == (' ', 'c')

def test_raw_multichar_uses_python_engine():
    """Any unrecognized multi-char value uses the python engine."""
    assert resolve_sep('::') == ('::', 'python')


# ── Auto-detect from extension ────────────────────────────────────────────────

def test_no_config_tsv_extension():
    assert resolve_sep(filepath='/data/vitals.tsv') == ('\t', 'c')

def test_no_config_csv_extension():
    assert resolve_sep(filepath='/data/labs.csv') == (',', 'c')

def test_no_config_txt_extension():
    """.txt defaults to tab, not whitespace — PMBB .txt files are tab-delimited."""
    assert resolve_sep(filepath='/data/codes.txt') == ('\t', 'c')

def test_no_config_unknown_extension():
    assert resolve_sep(filepath='/data/file.dat') == (',', 'c')

def test_no_config_no_filepath():
    assert resolve_sep() == (',', 'c')

def test_config_overrides_extension():
    """Explicit sep_config wins over file extension."""
    assert resolve_sep('tsv', filepath='/data/file.csv') == ('\t', 'c')
