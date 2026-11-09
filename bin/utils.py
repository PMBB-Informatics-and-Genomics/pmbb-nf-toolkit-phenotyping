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


def categorical_encoding(cfg):
    """
    Return the encoding a categorical phenotype config resolves to.

    The first match wins: ``dictionary`` → ``one_hot: true`` → ``binary: true``
    → ``auto`` (the default). Shared by resolve_config (which records it in
    each categorical JSON) and basic_categorical (which dispatches on it).
    """
    if cfg.get('dictionary'):
        return 'dictionary'
    if cfg.get('one_hot'):
        return 'one_hot'
    if cfg.get('binary'):
        return 'binary'
    return 'auto'


# Quantitative per-sample stats. `first` / `last` pick the value at the
# earliest / latest occurrence_date, so they need real dates (see needs_dates).
QUANT_STATS = ('mean', 'median', 'std', 'min', 'max', 'count', 'squared', 'first', 'last')
DATED_STATS = ('first', 'last')


def parse_stats(stats):
    """Return stats as a list of lowercase names; a string is comma-split."""
    if isinstance(stats, str):
        stats = stats.split(',')
    return [str(s).strip().lower() for s in (stats or []) if str(s).strip()]


def unknown_stats_error(stats):
    """Return an error string naming any unsupported stats, or None."""
    bad = [s for s in parse_stats(stats) if s not in QUANT_STATS]
    if not bad:
        return None
    return (f"stats: unknown {', '.join(repr(s) for s in bad)} — "
            f"valid: {', '.join(sorted(QUANT_STATS))}")


def needs_dates(cfg):
    """True if a quantitative config requests first/last, which need real dates."""
    return (cfg.get('data_type') == 'quantitative'
            and any(s in DATED_STATS for s in parse_stats(cfg.get('stats'))))
