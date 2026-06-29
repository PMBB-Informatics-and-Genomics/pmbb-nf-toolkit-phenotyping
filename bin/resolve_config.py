#!/usr/bin/env python3
"""
resolve_config.py — Resolve biobank YAML into per-phenotype JSON files.

Usage:
    resolve_config.py --config pmbb_config.yaml --output_dir ./configs/
"""

import argparse
import fnmatch
import json
import os
import sys
from copy import deepcopy

import pandas as pd
import yaml

from utils import resolve_sep

DATA_TYPES = {'quantitative', 'categorical', 'coded', 'advanced_coded'}

QUANT_DEFAULTS = {
    'stats': ['mean', 'median', 'std'],
    'qcut_bins': 0,
    'scale': False,
    'scale_method': 'zscore',
    'outlier_method': 'iqr',
    'outlier_mode': 'cap',
    'outlier_iqr_multiplier': 1.5,
    'outlier_zscore_sd': 3.0,
    'min_occurrences': 1,
}

CAT_DEFAULTS = {
    'missing_as_control': False,
    'one_hot': False,
    'min_occurrences': 1,
    'dictionary': None,
    'categories': None,
}

CODED_DEFAULTS = {
    'missing_as_control': False,
    'min_occurrences': 1,
    'subsample': 'all',
    'prefix_matching': True,
}

ADVANCED_CODED_DEFAULTS = {
    'missing_as_control': False,
    'min_occurrences': 1,
    'prefix_matching': True,
}

GLOBAL_DEFAULTS = {
    'sample_id_col': 'sample_ID',
    'date_col': None,
    'reference_date': 'today',
    'filter': {},
    'preprocessed_path': None,
    'output_name': None,
    'sep': None,
    'event_mask': [],
}


def merge_params(global_params, data_type_params, pheno_params, data_type):
    """Merge: built-in data_type defaults → global → data_type → phenotype."""
    if data_type == 'quantitative':
        base = deepcopy(QUANT_DEFAULTS)
    elif data_type == 'categorical':
        base = deepcopy(CAT_DEFAULTS)
    elif data_type == 'advanced_coded':
        base = deepcopy(ADVANCED_CODED_DEFAULTS)
    else:
        base = deepcopy(CODED_DEFAULTS)

    base.update(deepcopy(GLOBAL_DEFAULTS))
    base.update(deepcopy(global_params))
    base.update(deepcopy(data_type_params))
    base.update(deepcopy(pheno_params))
    return base


def validate(name, params):
    """Return list of error strings for a resolved phenotype config."""
    errors = []
    data_type = params.get('data_type', '')
    if data_type not in DATA_TYPES:
        errors.append(f"'data_type' must be one of {sorted(DATA_TYPES)}, got '{data_type}'")
        return errors

    pm = params.get('prefix_matching')
    if pm is not None and not isinstance(pm, bool):
        errors.append(f"'prefix_matching' must be a bool, got {type(pm).__name__!r}")

    if data_type == 'advanced_coded':
        sources = params.get('sources') or {}
        if not sources:
            errors.append("'sources' is required for advanced_coded phenotypes")
        for src_name, src in sources.items():
            if not src.get('from_phenotype'):
                errors.append(f"source '{src_name}': 'from_phenotype' is required")
            if not src.get('case_codes'):
                errors.append(f"source '{src_name}': 'case_codes' is required")
            src_pm = src.get('prefix_matching')
            if src_pm is not None and not isinstance(src_pm, bool):
                errors.append(
                    f"source '{src_name}': 'prefix_matching' must be a bool, got {type(src_pm).__name__!r}"
                )
        if '.' in name:
            print(
                f"WARNING: advanced_coded phenotype '{name}' contains a dot. "
                f"The name becomes an output column name — consider using underscores instead.",
                file=sys.stderr,
            )
    else:
        if not params.get('preprocessed_path'):
            if not params.get('table'):
                errors.append("'table' is required when 'preprocessed_path' is not set")
            if not params.get('value_col'):
                errors.append("'value_col' is required when 'preprocessed_path' is not set")
    return errors


def _get_coded_codes(resolved):
    """Read the coded table and return unique codes after applying filter and subsample.

    Reads in chunks to avoid loading the full table into memory.
    """
    value_col = resolved['value_col']
    filter_spec = resolved.get('filter') or {}
    filter_cols = list(filter_spec.keys())
    needed_cols = list(dict.fromkeys([value_col] + filter_cols))

    unique_codes = set()
    sep, engine = resolve_sep(resolved.get('sep'), resolved['table'])
    for chunk in pd.read_csv(resolved['table'], sep=sep, engine=engine, dtype=str,
                             usecols=needed_cols, chunksize=100_000):
        for col, patterns in filter_spec.items():
            if col not in chunk.columns:
                continue
            chunk = chunk[chunk[col].astype(str).apply(
                lambda v: any(fnmatch.fnmatch(v, str(p)) for p in patterns)
            )]
        unique_codes.update(chunk[value_col].dropna().tolist())

    unique_codes = list(unique_codes)

    subsample = resolved.get('subsample', 'all')
    if subsample == 'all':
        return unique_codes

    return [c for c in unique_codes
            if any(fnmatch.fnmatch(str(c), str(p)) for p in subsample)]


def resolve(config_path, output_dir, coded_chunk_size=0):
    """Read biobank YAML and write one JSON per phenotype."""
    config_path = os.path.abspath(config_path)
    config_dir = os.path.dirname(config_path)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    global_params = config.get('global', {}) or {}
    data_types_block = config.get('data_types', {}) or {}
    phenotypes = config.get('phenotypes', {}) or {}

    if not phenotypes:
        print('WARNING: no phenotypes defined in config', file=sys.stderr)
        return

    risk_rules = []

    for name, pheno_params in phenotypes.items():
        pheno_params = pheno_params or {}
        data_type = pheno_params.get('data_type') or global_params.get('data_type')
        if not data_type:
            print(f"ERROR: phenotype '{name}' has no data_type", file=sys.stderr)
            sys.exit(1)

        data_type_params = data_types_block.get(data_type, {}) or {}
        resolved = merge_params(global_params, data_type_params, pheno_params, data_type)
        resolved['phenotype_name'] = name
        resolved['data_type'] = data_type
        if not resolved.get('output_name'):
            resolved['output_name'] = name

        # Resolve relative table path to absolute using YAML file's directory
        if resolved.get('table') and not os.path.isabs(resolved['table']):
            resolved['table'] = os.path.abspath(
                os.path.join(config_dir, resolved['table'])
            )

        # Same for preprocessed_path
        if resolved.get('preprocessed_path') and not os.path.isabs(resolved['preprocessed_path']):
            resolved['preprocessed_path'] = os.path.abspath(
                os.path.join(config_dir, resolved['preprocessed_path'])
            )

        errors = validate(name, resolved)
        if errors:
            for e in errors:
                print(f"ERROR: phenotype '{name}': {e}", file=sys.stderr)
            sys.exit(1)

        out_subdir = os.path.join(output_dir, data_type)
        os.makedirs(out_subdir, exist_ok=True)

        # advanced_coded: no table paths to resolve, no chunking — write JSON as-is
        if data_type == 'advanced_coded':
            out_path = os.path.join(out_subdir, f'{name}.json')
            with open(out_path, 'w') as f:
                json.dump(resolved, f, indent=2, default=str)
            print(f'Wrote: {out_path}', file=sys.stderr)
            continue

        # Diagnostic chunking: split unique codes across N JSON files
        if data_type == 'coded' and coded_chunk_size > 0 and not resolved.get('preprocessed_path'):
            codes = _get_coded_codes(resolved)
            chunks = [codes[i:i + coded_chunk_size]
                      for i in range(0, len(codes), coded_chunk_size)]
            for idx, chunk in enumerate(chunks, start=1):
                chunk_resolved = deepcopy(resolved)
                chunk_resolved['subsample'] = chunk
                out_path = os.path.join(out_subdir, f'{name}_chunk_{idx:03d}.json')
                with open(out_path, 'w') as f:
                    json.dump(chunk_resolved, f, indent=2, default=str)
                print(f'Wrote: {out_path} ({len(chunk)} codes)', file=sys.stderr)
            continue

        out_path = os.path.join(out_subdir, f'{name}.json')
        with open(out_path, 'w') as f:
            json.dump(resolved, f, indent=2, default=str)
        print(f'Wrote: {out_path}', file=sys.stderr)

        # Collect risk classification rule for quantitative phenotypes
        if data_type == 'quantitative':
            rc = resolved.get('risk_classification')
            if rc and isinstance(rc, dict):
                method = rc.get('method')
                if method:
                    stat = rc.get('stat', 'mean')
                    output_name = resolved['output_name']
                    rule = {
                        'phenotype': name,
                        'stat': stat,
                        'source_col': f'{output_name}_{stat}',
                        'output_col': f'{output_name}_risk',
                        'method': method,
                    }
                    for field in ('high', 'low', 'cutoffs', 'n', 'labels'):
                        if field in rc:
                            rule[field] = rc[field]
                    risk_rules.append(rule)

    if risk_rules:
        risk_config_path = os.path.join(output_dir, 'risk_classification_config.json')
        with open(risk_config_path, 'w') as f:
            json.dump({'rules': risk_rules}, f, indent=2)
        print(f'Wrote: {risk_config_path}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description='Resolve biobank YAML config into per-phenotype JSONs')
    p.add_argument('--config', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--coded_chunk_size', type=int, default=0,
                   help='Split coded phenotypes into chunks of this size (0 = disabled)')
    args = p.parse_args()
    resolve(args.config, args.output_dir, coded_chunk_size=args.coded_chunk_size)


if __name__ == '__main__':
    main()
