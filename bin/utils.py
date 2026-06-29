#!/usr/bin/env python3
"""
Shared utilities for pmbb-nf-toolkit-phenotyping.
"""
import os

_SEP_ALIASES = {
    'tsv':        ('\t',    'c'),
    'tab':        ('\t',    'c'),
    'csv':        (',',     'c'),
    'comma':      (',',     'c'),
    'pipe':       ('|',     'c'),
    'space':      (' ',     'c'),
    'whitespace': (r'\s+',  'python'),
}

_EXT_DEFAULTS = {
    '.tsv': ('\t', 'c'),
    '.csv': (',',  'c'),
    '.txt': ('\t', 'c'),
}


def resolve_sep(sep_config=None, filepath=None):
    """
    Resolve the pandas separator and engine for reading a delimited file.

    sep_config : str or None
        Value from config. Recognized named aliases: tsv, tab, csv, comma,
        pipe, space, whitespace. Any other value is used as the raw separator
        character. Unrecognized separator characters must be double-quoted in
        YAML configs (e.g. ``sep: "|"`` not ``sep: |``).
    filepath : str or None
        Used for extension-based auto-detection when sep_config is None.
        Extensions: .tsv -> tab, .csv -> comma, .txt -> tab,
        anything else -> comma.

    Returns
    -------
    (sep, engine) : tuple[str, str]
        Pass both to pd.read_csv:
        ``pd.read_csv(path, sep=sep, engine=engine, ...)``

    Notes
    -----
    ``whitespace`` (and ``\\s+``) require ``engine='python'`` because the
    C engine does not support arbitrary multi-character regex separators.
    Do NOT use ``whitespace`` for TSV files whose field values may contain spaces.
    """
    if sep_config is not None:
        raw = str(sep_config)
        s = raw.strip()
        if s in _SEP_ALIASES:
            return _SEP_ALIASES[s]
        if s == r'\t':
            return '\t', 'c'
        if s == r'\s+':
            return r'\s+', 'python'
        engine = 'python' if len(raw) > 1 else 'c'
        return raw, engine

    if filepath:
        ext = os.path.splitext(filepath)[1].lower()
        if ext in _EXT_DEFAULTS:
            return _EXT_DEFAULTS[ext]

    return ',', 'c'
