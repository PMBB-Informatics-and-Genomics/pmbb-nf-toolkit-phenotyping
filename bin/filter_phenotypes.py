#!/usr/bin/env python3
"""
filter_phenotypes.py — Opt-in post-processing filter for the gathered table.

Usage:
    filter_phenotypes.py --config filter.yaml \
                         --input phenotype_table.wide.tsv \
                         --output_dir filtered \
                         [--covariates eda/eda_covariates.tsv] \
                         [--summary eda/eda_summary.tsv] \
                         [--manifest phenotype_manifest.tsv]

The pipeline's own outputs stay raw; this step writes NEW files only:

    phenotype_table.filtered.wide.tsv   the input minus dropped columns, with
                                        out-of-scope samples set to NA
    filter_log.tsv                      one row per affected column: action,
                                        reason, values / cases set to NA
    phenotype_manifest.filtered.tsv     manifest rows for kept columns (--manifest)

Config keys (start from eda/filter_suggested.yaml, written by --run_eda):

    drop_columns     list of column names / globs to remove
    keep_columns     list of names / globs that are always kept (wins over drops)
    drop_flags       drop every column carrying any of these EDA flags (--summary)
    restrictions     {column or glob: {covariate: C, level: L}} — samples whose
                     covariate C is not L are set to NA (cases included, counted
                     in the log as likely miscodes). L may be a list of levels.
    restrict_separated {min_share: S, covariates: [C, ...]} — restrict every column
                     whose cases the EDA found confined to a level (or subset of
                     levels) of C with share >= S (--summary sep_detail). Omit
                     `covariates` for all. Skipped for a covariate a column
                     already has an explicit `restrictions` entry for. The EDA
                     only records separation at or above its separation_threshold.
    restrict_unknown na (default) | keep — samples with no value for C

Covariate levels come from eda/eda_covariates.tsv, the discretized table the
EDA used, so `level` is the same code shown in restriction_suggestions.tsv
(e.g. "0" for SEX female when its dictionary maps female: 0).

Values are read and written as text, so kept cells are byte-identical to the
input; NA is written as an empty cell, as GATHER does.
"""

import argparse
import fnmatch
import os
import sys

import pandas as pd
import yaml

from eda_stats import LEVEL_SEP, parse_sep_detail

CONFIG_KEYS = {'drop_columns', 'keep_columns', 'drop_flags', 'restrictions',
               'restrict_separated', 'restrict_unknown'}
CASE_VALUES = {'1', '1.0'}


def load_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    unknown = sorted(set(cfg) - CONFIG_KEYS)
    if unknown:
        raise ValueError(f'filter config: unknown key(s) {unknown} — valid: {sorted(CONFIG_KEYS)}')
    # Keys whose entries are all commented out parse as None
    for key in ('drop_columns', 'keep_columns', 'drop_flags'):
        cfg[key] = [str(x) for x in (cfg.get(key) or [])]
    cfg['restrictions'] = cfg.get('restrictions') or {}
    cfg['restrict_unknown'] = cfg.get('restrict_unknown') or 'na'
    if cfg['restrict_unknown'] not in ('na', 'keep'):
        raise ValueError("filter config: restrict_unknown must be 'na' or 'keep'")
    for pattern, rule in cfg['restrictions'].items():
        if not isinstance(rule, dict) or set(rule) != {'covariate', 'level'}:
            raise ValueError(f'filter config: restrictions[{pattern!r}] must be '
                             '{covariate: ..., level: ...}')
        levels = rule['level'] if isinstance(rule['level'], list) else [rule['level']]
        if not levels:
            raise ValueError(f'filter config: restrictions[{pattern!r}].level is empty')
        rule['level'] = [str(x) for x in levels]
    rs = cfg.get('restrict_separated')
    if rs is not None:
        if not isinstance(rs, dict) or 'min_share' not in rs or set(rs) - {'min_share', 'covariates'}:
            raise ValueError('filter config: restrict_separated must be '
                             '{min_share: S, covariates: [...]} (covariates optional)')
        s = rs['min_share']
        if isinstance(s, bool) or not isinstance(s, (int, float)) or not 0 < s <= 1:
            raise ValueError(f'filter config: restrict_separated.min_share must be in (0, 1], got {s!r}')
        covs = rs.get('covariates') or []
        if not isinstance(covs, list):
            raise ValueError('filter config: restrict_separated.covariates must be a list')
        cfg['restrict_separated'] = {'min_share': float(s), 'covariates': [str(c) for c in covs]}
    return cfg


def _match(columns, patterns):
    """{pattern: [matching columns]} using exact names or fnmatch globs."""
    return {p: [c for c in columns if c == p or fnmatch.fnmatchcase(c, p)] for p in patterns}


def run_filter(config_path, input_path, output_dir, covariates_path=None,
               summary_path=None, manifest_path=None):
    cfg = load_config(config_path)
    os.makedirs(output_dir, exist_ok=True)

    wide = pd.read_csv(input_path, sep='\t', dtype=str, keep_default_na=False)
    if 'sample_ID' not in wide.columns:
        raise ValueError(f'{input_path}: no sample_ID column')
    pheno_cols = [c for c in wide.columns if c != 'sample_ID']
    log = {}

    def warn_unmatched(kind, matches):
        for p, cols in matches.items():
            if not cols:
                print(f'WARNING: {kind} {p!r} matches no column', file=sys.stderr)

    # ── Drops ─────────────────────────────────────────────────────────────────
    drop = {}
    m = _match(pheno_cols, cfg['drop_columns'])
    warn_unmatched('drop_columns', m)
    for p, cols in m.items():
        for c in cols:
            drop.setdefault(c, 'drop_columns' if c == p else f'drop_columns ~ {p}')

    if cfg['drop_flags']:
        if not summary_path:
            raise ValueError('drop_flags needs --summary (eda/eda_summary.tsv)')
        summary = pd.read_csv(summary_path, sep='\t', dtype=str, keep_default_na=False)
        wanted = set(cfg['drop_flags'])
        for col, flags in zip(summary['column_name'], summary['flags']):
            hit = sorted(wanted & set(filter(None, flags.split(','))))
            if hit and col in wide.columns:
                drop.setdefault(col, 'drop_flags: ' + ','.join(hit))

    keep_m = _match(pheno_cols, cfg['keep_columns'])
    warn_unmatched('keep_columns', keep_m)
    for c in {c for cols in keep_m.values() for c in cols}:
        if c in drop:
            log[c] = {'action': 'kept', 'reason': f'keep_columns overrides {drop.pop(c)}'}

    for c, reason in drop.items():
        log[c] = {'action': 'dropped', 'reason': reason}
    out = wide.drop(columns=list(drop))

    # ── Restrictions ──────────────────────────────────────────────────────────
    remaining = [c for c in out.columns if c != 'sample_ID']
    rules = {}   # column → [(covariate, {levels}, reason)]
    m = _match(remaining, list(cfg['restrictions']))
    warn_unmatched('restrictions', m)
    for p, cols in m.items():
        rule = cfg['restrictions'][p]
        for c in cols:
            rules.setdefault(c, []).append(
                (str(rule['covariate']), set(rule['level']),
                 f"{rule['covariate']} == {LEVEL_SEP.join(rule['level'])}"))

    rs = cfg.get('restrict_separated')
    if rs:
        if not summary_path:
            raise ValueError('restrict_separated needs --summary (eda/eda_summary.tsv)')
        summary = pd.read_csv(summary_path, sep='\t', dtype=str, keep_default_na=False)
        if 'sep_detail' not in summary.columns:
            raise ValueError(f'{summary_path}: no sep_detail column (not an eda_summary.tsv?)')
        wanted = set(rs['covariates'])
        for col, detail in zip(summary['column_name'], summary['sep_detail']):
            if col not in remaining:
                continue
            explicit = {cov for cov, _, _ in rules.get(col, [])}
            for s in parse_sep_detail(detail or None):
                if (s['outcome'] != 'cases' or s['share'] < rs['min_share']
                        or (wanted and s['covariate'] not in wanted)
                        or s['covariate'] in explicit):
                    continue
                rules.setdefault(col, []).append(
                    (s['covariate'], set(s['levels']),
                     f"{s['covariate']} == {LEVEL_SEP.join(s['levels'])} "
                     f"(separated, share {s['share']:.3f})"))
                explicit.add(s['covariate'])

    if rules:
        if not covariates_path:
            raise ValueError('restrictions need --covariates (eda/eda_covariates.tsv)')
        cov = pd.read_csv(covariates_path, sep='\t', dtype=str, keep_default_na=False)
        cov = cov.drop_duplicates('sample_ID').set_index('sample_ID')
        for c, col_rules in rules.items():
            out_of_scope = pd.Series(False, index=out.index)
            for cov_name, levels, _ in col_rules:
                if cov_name not in cov.columns:
                    raise ValueError(f'restriction on {c!r}: covariate {cov_name!r} not in '
                                     f'{covariates_path} (has: {list(cov.columns)})')
                lv = out['sample_ID'].map(cov[cov_name]).fillna('')
                unknown = lv == ''
                out_of_scope |= ~lv.isin(levels) & ~unknown
                if cfg['restrict_unknown'] == 'na':
                    out_of_scope |= unknown
            hit = out_of_scope & (out[c] != '')
            n_cases = int((hit & out[c].isin(CASE_VALUES)).sum())
            out.loc[out_of_scope, c] = ''
            log[c] = {'action': 'restricted',
                      'reason': '; '.join(r for _, _, r in col_rules)
                                + f' (restrict_unknown: {cfg["restrict_unknown"]})',
                      'n_set_na': int(hit.sum()), 'n_cases_set_na': n_cases}
            if n_cases:
                print(f'NOTE: {c}: {n_cases} cases outside '
                      f'{"; ".join(r for _, _, r in col_rules)} set to NA (likely miscodes)',
                      file=sys.stderr)

    # ── Outputs ───────────────────────────────────────────────────────────────
    table_path = os.path.join(output_dir, 'phenotype_table.filtered.wide.tsv')
    out.to_csv(table_path, sep='\t', index=False)

    log_df = pd.DataFrame(
        [{'column_name': c, **v} for c, v in log.items()],
        columns=['column_name', 'action', 'reason', 'n_set_na', 'n_cases_set_na'])
    log_df = log_df.sort_values(['action', 'column_name'])
    for col in ('n_set_na', 'n_cases_set_na'):
        log_df[col] = log_df[col].astype('Int64')
    log_df.to_csv(os.path.join(output_dir, 'filter_log.tsv'), sep='\t', index=False)

    if manifest_path:
        manifest = pd.read_csv(manifest_path, sep='\t', dtype=str, keep_default_na=False)
        manifest[manifest['column_name'].isin(out.columns)].to_csv(
            os.path.join(output_dir, 'phenotype_manifest.filtered.tsv'), sep='\t', index=False)

    n_restricted = int((log_df['action'] == 'restricted').sum())
    print(f'Filter: {len(pheno_cols)} columns in, {len(drop)} dropped, '
          f'{n_restricted} restricted → {table_path}', file=sys.stderr)
    return out, log_df


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True, help='Filter YAML (see eda/filter_suggested.yaml)')
    p.add_argument('--input', required=True, help='phenotype_table.wide.tsv')
    p.add_argument('--output_dir', default='.')
    p.add_argument('--covariates', default=None, help='eda/eda_covariates.tsv (for restrictions)')
    p.add_argument('--summary', default=None, help='eda/eda_summary.tsv (for drop_flags)')
    p.add_argument('--manifest', default=None, help='phenotype_manifest.tsv')
    args = p.parse_args()
    try:
        run_filter(args.config, args.input, args.output_dir, args.covariates,
                   args.summary, args.manifest)
    except ValueError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
