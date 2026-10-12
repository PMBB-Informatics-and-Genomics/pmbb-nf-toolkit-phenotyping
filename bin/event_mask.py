#!/usr/bin/env python3
"""
event_mask.py — Compute per-sample date-range mask windows from event ICD codes.

Reads long.tsvs produced by preprocess (one file per ICD code), finds rows
matching event code patterns, and computes [occurrence_date + offset_before,
occurrence_date + offset_after] windows. Overlapping windows per sample per
event type are merged.

Output: mask_intervals.tsv with columns:
    sample_ID, event_type, window_start, window_end

NOTE: This module always reads raw (unmasked) long.tsvs. Event window computation
is intentionally based on unmasked data — see Solution 1 in design docs.

Usage:
    event_mask.py \\
        --config event_mask.yaml \\
        --source ICD10 \\
        --input O83.3.long.tsv Z37.0.long.tsv O82.0.long.tsv ... \\
        --output ICD10.mask_intervals.tsv
"""

import argparse
import os
import sys
from collections import defaultdict
from datetime import timedelta

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(__file__))
from icd_utils import normalize_patterns, code_matches_pattern, build_prefix_index


_SCHEMA = ['sample_ID', 'event_type', 'window_start', 'window_end']


def _code_from_path(tsv_path: str) -> str:
    """
    Extract ICD code from a long.tsv filename.
    'path/to/O83.3.long.tsv' → 'O83.3'
    'path/to/E11.long.tsv'   → 'E11'
    """
    basename = os.path.basename(tsv_path)
    if basename.lower().endswith('.long.tsv'):
        return basename[: -len('.long.tsv')]
    return basename.split('.')[0]


def _merge_intervals(intervals: list) -> list:
    """
    Merge overlapping or adjacent (start, end) date pairs.
    Input list of (date, date) tuples; returns sorted, merged list.
    """
    if not intervals:
        return []
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [list(sorted_ivs[0])]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def run_event_mask(
    config_path: str,
    source_name: str,
    input_tsvs: list,
    output_path: str,
) -> None:
    """
    Compute mask intervals for all events in config_path that reference source_name.

    Parameters
    ----------
    config_path  : path to event_mask.yaml
    source_name  : phenotype name of the source (e.g., 'ICD10')
    input_tsvs   : list of long.tsv paths from that source phenotype
    output_path  : path to write mask_intervals.tsv
    """
    with open(config_path) as f:
        config = yaml.safe_load(f) or {}

    events = config.get('events') or {}

    # Only process events that reference this source
    relevant_events = {
        name: cfg for name, cfg in events.items()
        if cfg.get('source') == source_name
    }

    if not relevant_events:
        _write_empty(output_path)
        print(f'No events reference source "{source_name}"; wrote empty intervals.', file=sys.stderr)
        return

    # Pre-index input TSVs by their code (extracted from filename).
    # Use a list per code so duplicate basenames (e.g., two runs producing
    # the same O83.3.long.tsv) are all retained rather than silently dropped.
    tsv_by_code: dict = defaultdict(list)
    for p in input_tsvs:
        tsv_by_code[_code_from_path(p)].append(p)

    all_rows = []

    for event_name, event_cfg in relevant_events.items():
        raw_codes = event_cfg.get('codes') or []
        offset_before = int(event_cfg.get('offset_before', 0))  # days, typically negative
        offset_after = int(event_cfg.get('offset_after', 0))    # days, typically positive

        patterns = normalize_patterns(raw_codes)
        prefix_idx = build_prefix_index(patterns)

        # Also build a raw-code prefix index for bare-code filenames.
        # normalize_patterns('O03') → 'O03.*' which prefix-matches 'O03.4' but not 'O03'
        # (a file named 'O03.long.tsv' has code 'O03', no dot).  We handle this by also
        # checking the raw event codes directly: filename code must start with a raw code.
        # NOTE: using an explicit wildcard like 'O03.*' in the config also fails to match
        # bare-code filenames — the prefix becomes 'O03.' which 'O03' doesn't start with.
        # Users must specify bare 'O03' (not 'O03.*') to match bare-code filenames.
        raw_prefix_idx = build_prefix_index(raw_codes)

        # Find matching TSVs by code prefix matching against filename-extracted codes
        matching_tsvs = []
        for code, paths in tsv_by_code.items():
            key = code[:2] if len(code) >= 2 else code
            # Try normalized patterns first (covers dotted codes like 'O83.3')
            candidate_pats = prefix_idx.get(key, [])
            matched = candidate_pats and any(
                code_matches_pattern(code, pat) for pat in candidate_pats
            )
            if not matched:
                # Fallback: match raw codes as prefixes of the filename code.
                # This handles bare event codes ('O03') matching 'O03.long.tsv' (code='O03').
                raw_candidates = raw_prefix_idx.get(key, [])
                matched = any(code.startswith(raw_code) for raw_code in raw_candidates)
            if matched:
                matching_tsvs.extend(paths)

        if not matching_tsvs:
            print(
                f'WARNING: event "{event_name}": no TSVs matched codes {raw_codes} '
                f'in source "{source_name}".',
                file=sys.stderr,
            )
            continue

        # Read matching TSVs, compute windows per sample
        sample_intervals: dict = {}  # {sample_ID: [(start, end), ...]}
        for tsv_path in matching_tsvs:
            df = pd.read_csv(tsv_path, sep='\t', dtype=str)
            if 'sample_ID' not in df.columns or 'occurrence_date' not in df.columns:
                print(f'WARNING: {tsv_path} missing required columns; skipping.', file=sys.stderr)
                continue

            df['_date'] = pd.to_datetime(df['occurrence_date'], errors='coerce')
            df = df.dropna(subset=['_date'])
            if df.empty:
                continue
            df['_start'] = df['_date'].dt.date.apply(lambda d: d + timedelta(days=offset_before))
            df['_end']   = df['_date'].dt.date.apply(lambda d: d + timedelta(days=offset_after))
            for sid, grp in df.groupby('sample_ID'):
                sample_intervals.setdefault(sid, []).extend(zip(grp['_start'], grp['_end']))

        # Merge overlapping windows per sample for this event type
        for sid, intervals in sample_intervals.items():
            for start, end in _merge_intervals(intervals):
                all_rows.append({
                    'sample_ID': sid,
                    'event_type': event_name,
                    'window_start': start.isoformat(),
                    'window_end': end.isoformat(),
                })

    result = pd.DataFrame(all_rows, columns=_SCHEMA) if all_rows else pd.DataFrame(columns=_SCHEMA)
    result.to_csv(output_path, sep='\t', index=False)
    print(f'Wrote {len(result)} mask intervals to: {output_path}', file=sys.stderr)


def _write_empty(output_path: str) -> None:
    pd.DataFrame(columns=_SCHEMA).to_csv(output_path, sep='\t', index=False)


def main():
    p = argparse.ArgumentParser(description='Compute event mask intervals from ICD long.tsvs')
    p.add_argument('--config', required=True, help='event_mask.yaml config file')
    p.add_argument('--source', required=True, help='Source phenotype name (e.g., ICD10)')
    p.add_argument('--input', required=True, nargs='+', help='Long.tsv files from the source phenotype')
    p.add_argument('--output', required=True, help='Output mask_intervals.tsv path')
    args = p.parse_args()
    run_event_mask(args.config, args.source, args.input, args.output)


if __name__ == '__main__':
    main()
