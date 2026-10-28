#!/usr/bin/env python3
"""
resolve_config.py — Resolve biobank YAML into per-phenotype JSON files.

Usage:
    resolve_config.py --config pmbb_config.yaml --output_dir ./configs/
"""

import argparse
import fnmatch
import itertools
import json
import os
import sys
from copy import deepcopy

import pandas as pd
import yaml

from utils import categorical_encoding, needs_dates, resolve_sep, unknown_stats_error

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
    'eda': True,
}

CAT_DEFAULTS = {
    'missing_as_control': False,
    'one_hot': False,
    'binary': False,
    'min_occurrences': 1,
    'dictionary': None,
    'categories': None,
    'eda': True,
}

CODED_DEFAULTS = {
    'missing_as_control': False,
    'min_occurrences': 1,
    'subsample': 'all',
    'prefix_matching': True,
    # Prepended to each exploded basic_coded output column (and its .diag.tsv
    # filename) only. Long-format files keep raw codes, so advanced_coded and
    # event_mask are unaffected.
    'column_prefix': '',
    # Coded phenotypes explode into thousands of columns, so ydata profiles are
    # opt-in: eda: true, or a list of code globs (eda: ["E11*", I10]).
    'eda': False,
}

ADVANCED_CODED_DEFAULTS = {
    'missing_as_control': False,
    'min_occurrences': 1,
    'prefix_matching': True,
    'eda': False,
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
    'restrict_to': None,
}

# Top-level `eda:` block — report settings (used only with --run_eda). The
# per-phenotype `eda` key above only decides which columns get a ydata
# profile; metrics and flags are always computed for every column.
EDA_DEFAULTS = {
    'covariates': [],
    'covariates_file': None,
    'covariates_file_cols': None,
    'age_bins': 5,
    'max_covariate_levels': 10,
    'max_missing': 0.5,
    'min_cases': 100,
    'separation_threshold': 0.95,
    'separation_min_cases': 10,
    'max_subset_levels': 3,
    'separation_interactions': [],
    'min_level_n': 50,
    'min_expected': 5,
    'missing_gap': 0.8,
    'profile_flagged_top_k': 50,
    'max_profile_columns': 200,
    'profile_group_size': 50,
    'profile_max_samples': 20000,
    'profile_correlations': False,
    'flagged_table_max_rows': 5000,
    'pairwise': False,
    'pairwise_jaccard': 0.95,
    'pairwise_containment': 0.99,
    'pairwise_top_k': 1000,
    'pairwise_scope': 'all',
    'seed': 42,
}

# Removed / renamed eda: keys → hint shown with the unknown-key error.
EDA_RENAMED = {
    'quasi_threshold': "renamed to 'separation_threshold'",
    'sex_covariate': 'removed: restrictions are data-driven (any covariate) or restrict_to',
    'restriction_rules': 'removed: see community/sex_specific/ for code lists to use in a filter config',
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

    if data_type == 'quantitative':
        stats_error = unknown_stats_error(params.get('stats'))
        if stats_error:
            errors.append(stats_error)
        # first/last order rows by date. Without a real date_col, PREPROCESS
        # fills every row with reference_date and the order is meaningless.
        if needs_dates(params) and not params.get('date_col') and not params.get('preprocessed_path'):
            errors.append("stats 'first'/'last' require 'date_col' (the measurement date column)")

    eda = params.get('eda')
    if eda is not None and not isinstance(eda, bool) and not (
            isinstance(eda, list) and all(isinstance(x, (str, int, float)) for x in eda)):
        errors.append(f"'eda' must be true, false, or a list of code/column globs, got {eda!r}")

    rt = params.get('restrict_to')
    if rt is not None and not (isinstance(rt, dict) and len(rt) == 1):
        errors.append(f"'restrict_to' must be a one-entry mapping {{covariate: level}}, got {rt!r}")

    pm = params.get('prefix_matching')
    if pm is not None and not isinstance(pm, bool):
        errors.append(f"'prefix_matching' must be a bool, got {type(pm).__name__!r}")

    cp = params.get('column_prefix')
    if data_type == 'coded' and cp is not None and not isinstance(cp, str):
        errors.append(f"'column_prefix' must be a string, got {type(cp).__name__!r}")

    if data_type == 'advanced_coded':
        if params.get('column_prefix'):
            errors.append(
                "'column_prefix' is not supported on advanced_coded phenotypes — "
                "set it on the source 'coded' phenotype instead; advanced_coded "
                "case_codes always use raw codes"
            )
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


def resolve_eda(eda_block, eda_phenotypes, global_params, config_dir):
    """
    Resolve the top-level `eda:` block into the eda_config.json payload.

    eda_phenotypes : {name: resolved phenotype params} for every phenotype.
    Exits with an error message on invalid settings.
    """
    eda_block = eda_block or {}
    if not isinstance(eda_block, dict):
        print("ERROR: 'eda' block must be a mapping", file=sys.stderr)
        sys.exit(1)
    unknown = sorted(set(eda_block) - set(EDA_DEFAULTS))
    if unknown:
        print(f"ERROR: eda: unknown setting(s) {unknown} — valid: {sorted(EDA_DEFAULTS)}",
              file=sys.stderr)
        for key in unknown:
            if key in EDA_RENAMED:
                print(f"  {key}: {EDA_RENAMED[key]}", file=sys.stderr)
        sys.exit(1)

    settings = deepcopy(EDA_DEFAULTS)
    settings.update(deepcopy(eda_block))

    errors = []
    covariates = settings['covariates'] or []
    if isinstance(covariates, str):
        covariates = [covariates]
    settings['covariates'] = covariates
    for cov in covariates:
        p = eda_phenotypes.get(cov)
        if p is None:
            errors.append(f"covariate '{cov}' is not a phenotype in this config")
        elif p['data_type'] not in ('quantitative', 'categorical'):
            errors.append(f"covariate '{cov}' must be quantitative or categorical, "
                          f"got {p['data_type']}")

    cov_file = settings.get('covariates_file')
    if cov_file and not os.path.isabs(cov_file):
        settings['covariates_file'] = os.path.abspath(os.path.join(config_dir, cov_file))
    cov_cols = settings.get('covariates_file_cols')
    if cov_cols is not None and not isinstance(cov_cols, list):
        errors.append("'covariates_file_cols' must be a list of column names")
    elif cov_file and cov_cols and os.path.isfile(settings['covariates_file']):
        sep, engine = resolve_sep(filepath=settings['covariates_file'])
        header = pd.read_csv(settings['covariates_file'], sep=sep, engine=engine,
                             nrows=0).columns
        missing = [c for c in cov_cols if c not in header]
        if missing:
            errors.append(f"covariates_file_cols {missing} not found in {cov_file}")

    if settings['pairwise_scope'] not in ('all', 'within_phenotype'):
        errors.append("pairwise_scope must be 'all' or 'within_phenotype'")
    for key in ('max_missing', 'separation_threshold', 'missing_gap',
                'pairwise_jaccard', 'pairwise_containment'):
        v = settings[key]
        if not isinstance(v, (int, float)) or not 0 <= v <= 1:
            errors.append(f"'{key}' must be a number in [0, 1], got {v!r}")
    for key, low in (('max_subset_levels', 1), ('separation_min_cases', 1)):
        k = settings[key]
        if isinstance(k, bool) or not isinstance(k, int) or k < low:
            errors.append(f"'{key}' must be an integer >= {low}, got {k!r}")

    # Covariate names / levels are joined with these in sep_level / sep_detail
    all_covs = list(covariates) + list(cov_cols or [])
    for cov in all_covs:
        if any(ch in str(cov) for ch in '&|=;:'):
            errors.append(f"covariate name '{cov}' must not contain any of & | = ; :")
    inter = settings['separation_interactions']
    if inter == 'auto':
        inter = [list(pair) for pair in itertools.combinations(all_covs, 2)]
    if not isinstance(inter, list) or not all(
            isinstance(c, list) and len(c) >= 2 and len(set(c)) == len(c) for c in inter):
        errors.append("'separation_interactions' must be 'auto' or a list of covariate "
                      "lists, e.g. [[COHORT, SEX]]")
    else:
        for combo in inter:
            unknown_covs = [c for c in combo if c not in all_covs]
            if unknown_covs:
                errors.append(f"separation_interactions {combo}: {unknown_covs} not in "
                              "eda.covariates / covariates_file_cols")
        settings['separation_interactions'] = inter

    known_covariates = set(covariates) | set(cov_cols or [])
    for name, p in eda_phenotypes.items():
        rt = p.get('restrict_to')
        if rt:
            cov = next(iter(rt))
            if cov_file and not cov_cols:
                continue  # all file columns are covariates; checked at runtime
            if cov not in known_covariates:
                errors.append(f"phenotype '{name}': restrict_to covariate '{cov}' "
                              f"is not listed in eda.covariates / covariates_file_cols")

    if errors:
        for e in errors:
            print(f"ERROR: eda: {e}", file=sys.stderr)
        sys.exit(1)

    dictionaries = {
        cov: eda_phenotypes[cov].get('dictionary')
        for cov in covariates if eda_phenotypes[cov].get('dictionary')
    }
    return {
        'settings': settings,
        'sample_id_col': global_params.get('sample_id_col', GLOBAL_DEFAULTS['sample_id_col']),
        'covariate_dictionaries': dictionaries,
        'phenotypes': {
            name: {
                'data_type': p['data_type'],
                'output_name': p.get('output_name') or name,
                'column_prefix': p.get('column_prefix') or '',
                'eda': p.get('eda', False),
                'restrict_to': p.get('restrict_to'),
            }
            for name, p in eda_phenotypes.items()
        },
    }


def resolve(config_path, output_dir, coded_chunk_size=0):
    """Read biobank YAML and write one JSON per phenotype."""
    # realpath (not abspath): Nextflow stages the config as a symlink in the
    # work dir, and relative table paths must resolve against the real file.
    config_path = os.path.realpath(config_path)
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
    eda_phenotypes = {}

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
        if data_type == 'categorical':
            # Record which encoding basic_categorical will apply, so the
            # published configs/ JSON shows the effective choice.
            resolved['encoding'] = categorical_encoding(resolved)

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
        eda_phenotypes[name] = resolved

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

    # Always written (cheap); only read by the EDA_* processes under --run_eda.
    eda_config = resolve_eda(config.get('eda'), eda_phenotypes, global_params, config_dir)
    eda_config_path = os.path.join(output_dir, 'eda_config.json')
    with open(eda_config_path, 'w') as f:
        json.dump(eda_config, f, indent=2, default=str)
    print(f'Wrote: {eda_config_path}', file=sys.stderr)


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
