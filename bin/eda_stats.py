#!/usr/bin/env python3
"""
eda_stats.py — Per-column EDA metrics for a batch of per-phenotype result TSVs.

Usage:
    eda_stats.py --inputs BMI.quant.tsv E11.9.diag.tsv ... \
                 --eda_config eda_config.json \
                 --output batch.eda.tsv \
                 [--covariates SEX.cat.tsv age.quant.tsv] \
                 [--universe sample_list.txt] \
                 [--cases_output batch.cases.tsv] \
                 [--covariates_output eda_covariates.tsv]

Writes one row per non-sample_ID column of every input (see STATS_COLUMNS).
Metrics are cheap and computed for every column; eda_report.py turns them into
flags and decides which columns get a (costly) ydata profile.

Covariates come from the per-phenotype results of phenotypes named in
eda.covariates (file basename = phenotype name) and/or eda.covariates_file.
Numeric covariates with more than max_covariate_levels values are binned into
age_bins quantiles. Each binary phenotype is cross-tabulated against each
covariate to detect (quasi-)complete separation — cases confined to one level,
or to a subset of up to max_subset_levels levels (sep_level "A|B") — and every
column's missingness is compared across covariate levels.

The sample universe (denominator for missingness) is --universe if given,
else the union of covariate samples, else each file's own rows — quantitative
results only list measured samples, so without a universe their missingness
reads as 0.

--covariates_output writes the discretized covariate table (sample_ID + one
level column per covariate) — the same level codes the report uses, so the
filter step (filter_phenotypes.py) applies restrictions against identical levels.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

from utils import resolve_sep

STATS_COLUMNS = [
    'column_name', 'source_file', 'phenotype', 'data_type', 'kind', 'is_covariate',
    'n_total', 'n', 'n_missing', 'missing_frac',
    'n_cases', 'n_controls', 'case_frac',
    'mean', 'sd', 'min', 'p01', 'p50', 'p99', 'max', 'n_unique',
    'n_levels', 'min_level_count',
    'sep_type', 'sep_covariate', 'sep_level', 'sep_outcome', 'sep_share', 'sep_detail',
    'miss_gap', 'miss_gap_covariate', 'miss_gap_detail',
]

_SUFFIX_TYPES = (
    ('.quant.tsv', 'quantitative'),
    ('.cat.tsv', 'categorical'),
    ('.diag.tsv', 'coded'),
)


def file_data_type(path):
    """Return (stem, data_type) from a result filename; data_type None if unknown."""
    base = os.path.basename(path)
    for suffix, data_type in _SUFFIX_TYPES:
        if base.endswith(suffix):
            return base[:-len(suffix)], data_type
    return base, None


def read_result(path):
    """Read a result TSV indexed by sample_ID (all value columns numeric)."""
    df = pd.read_csv(path, sep='\t', dtype={'sample_ID': str})
    if 'sample_ID' not in df.columns:
        return None
    df = df.drop_duplicates('sample_ID').set_index('sample_ID')
    return df.apply(pd.to_numeric, errors='coerce')


# Joins the levels of a multi-level separation in sep_level / sep_detail (e.g. "A|B").
LEVEL_SEP = '|'
# Joins the covariates of an interaction and their levels (COHORT&SEX = "1&0").
INTERACTION_SEP = '&'


def bin_covariate(series, settings):
    """Discretize a covariate: keep few-level columns, quantile-bin the rest."""
    s = series.dropna()
    if s.nunique() <= settings['max_covariate_levels']:
        out = series.map(lambda v: v if pd.isna(v) else _fmt_level(v))
        return out.astype(object)
    numeric = pd.to_numeric(series, errors='coerce')
    binned = pd.qcut(numeric, q=settings['age_bins'], duplicates='drop')
    return binned.astype(str).where(binned.notna(), None).astype(object)


def _fmt_level(v):
    """Render a level label: 1.0 → '1', 'female' stays."""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def load_covariates(cov_paths, eda_cfg):
    """
    Build the covariate table (sample_ID index, one discretized column each).

    Returns (covariates DataFrame, {covariate column: source result basename}).
    """
    settings = eda_cfg['settings']
    names = set(settings.get('covariates') or [])
    frames, sources = [], {}

    for path in cov_paths or []:
        stem, data_type = file_data_type(path)
        if stem not in names:
            continue
        df = read_result(path)
        if df is None or df.shape[1] == 0:
            continue
        if data_type == 'quantitative':
            df = df.iloc[:, :1]  # first stat column (e.g. age_mean)
        if df.shape[1] == 1:
            df.columns = [stem]
        for col in df.columns:
            sources[col] = os.path.basename(path)
        frames.append(df)

    cov_file = settings.get('covariates_file')
    if cov_file:
        sep, engine = resolve_sep(filepath=cov_file)
        id_col = eda_cfg.get('sample_id_col', 'sample_ID')
        ext = pd.read_csv(cov_file, sep=sep, engine=engine, dtype=str)
        if id_col not in ext.columns:
            raise ValueError(f"covariates_file {cov_file}: no '{id_col}' column")
        ext = ext.drop_duplicates(id_col).set_index(id_col)
        ext.index.name = 'sample_ID'
        cols = settings.get('covariates_file_cols') or list(ext.columns)
        ext = ext[cols]
        for col in cols:
            numeric = pd.to_numeric(ext[col], errors='coerce')
            if numeric.notna().sum() == ext[col].notna().sum():
                ext[col] = numeric
        frames.append(ext)

    if not frames:
        return pd.DataFrame(), sources
    cov = pd.concat(frames, axis=1, join='outer')
    cov = cov.loc[:, ~cov.columns.duplicated()]
    for col in cov.columns:
        cov[col] = bin_covariate(cov[col], settings)
    return add_interactions(cov, settings), sources


def add_interactions(cov, settings):
    """
    Append one column per eda.separation_interactions entry (resolve_config
    expands 'auto' to every pair): name and levels joined with INTERACTION_SEP,
    e.g. COHORT&SEX = "1&0"; missing if any part is missing. Interactions are
    used for separation only, not missing-by-covariate.
    """
    for combo in settings.get('separation_interactions') or []:
        missing = [c for c in combo if c not in cov.columns]
        if missing:
            print(f'WARNING: separation_interactions {combo}: no covariate {missing} — skipped',
                  file=sys.stderr)
            continue
        parts = cov[list(combo)]
        joined = parts.astype(str).agg(INTERACTION_SEP.join, axis=1)
        cov[INTERACTION_SEP.join(combo)] = joined.where(parts.notna().all(axis=1), None)
    return cov


def separation(n_by, cases_by, settings):
    """
    Detect case/control confinement in one 1×k covariate table.

    n_by, cases_by : arrays of non-missing count and case count per level.
    Returns (sep_type, level indices, outcome, share) or None.

    Looks for the smallest subset of at most max_subset_levels levels (the
    levels with the most outcomes first) holding >= separation_threshold of the
    cases (or controls):

    complete : all of them sit in the subset
    quasi    : >= separation_threshold do, while the subset holds less than
               separation_threshold of samples (so it is not just a dominant level)

    Both require at least separation_min_cases cases (controls) in total and
    min_expected of them expected outside the subset, so rare phenotypes are
    not flagged just for being rare, and at least one populated level must
    remain outside it.
    """
    N = n_by.sum()
    C = cases_by.sum()
    if N == 0 or C == 0 or C == N or (n_by >= settings['min_level_n']).sum() < 2:
        return None
    threshold = settings['separation_threshold']
    max_k = min(int(settings['max_subset_levels']), int((n_by > 0).sum()) - 1)
    best = None
    for outcome, counts, total in (('cases', cases_by, C), ('controls', n_by - cases_by, N - C)):
        if total < settings['separation_min_cases']:
            continue
        order = np.argsort(-counts, kind='stable')
        for k in range(1, max_k + 1):
            idx = order[:k]
            share = counts[idx].sum() / total
            if share < threshold:
                continue
            n_in = n_by[idx].sum()
            if (N - n_in) * total / N < settings['min_expected']:
                break
            if share >= 1.0:
                cand = ('complete', sorted(int(i) for i in idx), outcome, float(share))
            elif n_in / N < threshold:
                cand = ('quasi', sorted(int(i) for i in idx), outcome, float(share))
            else:
                break
            if best is None or _rank(cand) > _rank(best):
                best = cand
            break
    return best


def _explained_by_parts(sep, seps):
    """
    True for an interaction separation that adds nothing over a separation of
    one of its parts: its levels are exactly every present combination whose
    part falls in that part's separated levels (e.g. COHORT&SEX = 1&0|1&1 when
    COHORT = 1 alone already holds the cases).
    """
    cov = sep[4]
    if INTERACTION_SEP not in cov:
        return False
    parts = cov.split(INTERACTION_SEP)
    levels, present = set(sep[6]), sep[7]
    for other in seps:
        if other[4] in parts and other[2] == sep[2]:
            i = parts.index(other[4])
            covered = {p for p in present if p.split(INTERACTION_SEP)[i] in set(other[6])}
            if covered == levels:
                return True
    return False


def parse_sep_detail(detail):
    """
    Split a sep_detail cell ("SEX=0:complete:cases(1.000);COHORT=A|B:quasi:cases(0.970)")
    into dicts with covariate, levels (list), sep_type, outcome and share.
    """
    out = []
    if not isinstance(detail, str):
        return out
    for entry in filter(None, detail.split(';')):
        head, sep_type, tail = entry.rsplit(':', 2)
        covariate, level = head.split('=', 1)
        outcome, share = tail.rstrip(')').split('(')
        out.append({'covariate': covariate, 'levels': level.split(LEVEL_SEP),
                    'sep_type': sep_type, 'outcome': outcome, 'share': float(share)})
    return out


def _rank(sep):
    """Sort key: complete beats quasi, cases beat controls, fewer levels, then share."""
    sep_type, idx, outcome, share = sep[:4]
    return ({'complete': 2, 'quasi': 1}[sep_type], outcome == 'cases', -len(idx), share)


def compute_file_stats(df, source_file, data_type, cov, cov_sources, settings,
                       universe=None):
    """Return a list of stats dicts (one per column) for one result DataFrame."""
    if universe is not None:
        df = df.reindex(df.index.union(pd.Index(universe)))
    stem, _ = file_data_type(source_file)
    rows = []
    cov_aligned = cov.reindex(df.index) if len(cov.columns) else cov
    own_covs = {c for c, src in cov_sources.items() if src == source_file}

    for col in df.columns:
        y = df[col]
        nonnull = y.dropna()
        n_total = len(y)
        n = int(len(nonnull))
        row = dict.fromkeys(STATS_COLUMNS)
        row.update({
            'column_name': col, 'source_file': source_file, 'phenotype': stem,
            'data_type': data_type, 'is_covariate': bool(own_covs),
            'n_total': n_total, 'n': n, 'n_missing': n_total - n,
            'missing_frac': (n_total - n) / n_total if n_total else np.nan,
            'n_unique': int(nonnull.nunique()),
        })

        values = set(nonnull.unique())
        if data_type == 'coded' or (data_type == 'categorical' and values <= {0, 1}):
            kind = 'binary'
        elif data_type == 'categorical':
            kind = 'categorical'
        else:
            kind = 'quantitative'
        row['kind'] = kind

        if kind == 'binary':
            cases = int((nonnull == 1).sum())
            row.update({'n_cases': cases, 'n_controls': n - cases,
                        'case_frac': cases / n if n else np.nan})
        elif kind == 'categorical':
            counts = nonnull.value_counts()
            row.update({'n_levels': int(len(counts)),
                        'min_level_count': int(counts.min()) if len(counts) else 0})
        if kind == 'quantitative' and n:
            q = nonnull.quantile([0.01, 0.5, 0.99])
            row.update({'mean': nonnull.mean(), 'sd': nonnull.std(),
                        'min': nonnull.min(), 'p01': q[0.01], 'p50': q[0.5],
                        'p99': q[0.99], 'max': nonnull.max()})

        if len(cov_aligned.columns) and not own_covs:
            _covariate_metrics(row, y, kind, cov_aligned, settings)
        rows.append(row)
    return rows


def _covariate_metrics(row, y, kind, cov, settings):
    """Fill separation and missing-by-covariate fields of row in place."""
    seps, gaps = [], []
    for c in cov.columns:
        levels = cov[c]
        mask = levels.notna()
        if not mask.any():
            continue
        lv = levels[mask]
        yv = y[mask]

        # Missingness gap across levels large enough to trust
        size_by = lv.groupby(lv).size()
        miss_by = yv.isna().groupby(lv).mean()
        big = size_by[size_by >= settings['min_level_n']].index
        if len(big) >= 2 and INTERACTION_SEP not in str(c):
            m = miss_by.loc[big]
            gaps.append((float(m.max() - m.min()), c,
                         ';'.join(f'{k}={v:.3f}' for k, v in m.items())))

        if kind == 'binary' and row['n_controls'] and row['n_cases']:
            obs = yv.notna()
            n_by = obs.groupby(lv).sum()
            cases_by = (yv == 1).groupby(lv).sum().reindex(n_by.index)
            res = separation(n_by.to_numpy(float), cases_by.to_numpy(float), settings)
            if res:
                sep_type, idx, outcome, share = res
                levels = [str(n_by.index[i]) for i in idx]
                present = [str(l) for l, n in n_by.items() if n > 0]
                seps.append((sep_type, idx, outcome, share, c, LEVEL_SEP.join(levels),
                             levels, present))

    seps = [s[:6] for s in seps if not _explained_by_parts(s, seps)]
    if seps:
        seps.sort(key=_rank, reverse=True)
        _, _, outcome, share, c, level = seps[0]
        row.update({'sep_type': seps[0][0], 'sep_covariate': c, 'sep_level': level,
                    'sep_outcome': outcome, 'sep_share': share,
                    'sep_detail': ';'.join(f'{c}={lvl}:{t}:{o}({s:.3f})'
                                           for t, _, o, s, c, lvl in seps)})
    if gaps:
        gaps.sort(reverse=True)
        row.update({'miss_gap': gaps[0][0], 'miss_gap_covariate': gaps[0][1],
                    'miss_gap_detail': gaps[0][2]})


def run_eda_stats(input_paths, eda_config_path, output_path, covariate_paths=None,
                  universe_path=None, cases_output=None, covariates_output=None):
    with open(eda_config_path) as f:
        eda_cfg = json.load(f)
    settings = eda_cfg['settings']

    cov, cov_sources = load_covariates(covariate_paths, eda_cfg)
    if not len(cov.columns):
        print('WARNING: no covariates — separation and missing-by-covariate '
              'metrics skipped', file=sys.stderr)
    if covariates_output:
        cov.rename_axis('sample_ID').reset_index().to_csv(
            covariates_output, sep='\t', index=False)

    universe = None
    if universe_path:
        with open(universe_path) as f:
            universe = [line.strip() for line in f if line.strip()]
    elif len(cov.columns):
        universe = list(cov.index)

    rows = []
    cases_lines = []
    for path in input_paths:
        source_file = os.path.basename(path)
        _, data_type = file_data_type(path)
        if data_type is None:
            print(f'WARNING: {path}: unrecognized result suffix — skipping', file=sys.stderr)
            continue
        df = read_result(path)
        if df is None:
            print(f'WARNING: {path} has no sample_ID column — skipping', file=sys.stderr)
            continue
        file_rows = compute_file_stats(df, source_file, data_type, cov, cov_sources,
                                       settings, universe)
        rows.extend(file_rows)

        if cases_output and settings.get('pairwise'):
            for r in file_rows:
                if r['kind'] == 'binary' and (r['n_cases'] or 0) >= settings['min_cases'] \
                        and r['n_controls']:
                    ids = df.index[df[r['column_name']] == 1]
                    cases_lines.append(f"{r['column_name']}\t{source_file}\t{','.join(ids)}")

    pd.DataFrame(rows, columns=STATS_COLUMNS).to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path} ({len(rows)} columns)', file=sys.stderr)

    if cases_output:
        with open(cases_output, 'w') as f:
            f.write('column_name\tsource_file\tcase_ids\n')
            for line in cases_lines:
                f.write(line + '\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inputs', nargs='+', required=True, help='Per-phenotype result TSVs')
    p.add_argument('--eda_config', required=True, help='eda_config.json from resolve_config')
    p.add_argument('--output', required=True, help='Output stats TSV')
    p.add_argument('--covariates', nargs='*', default=[],
                   help='Result TSVs of covariate phenotypes')
    p.add_argument('--universe', default=None, help='Sample list (one ID per line)')
    p.add_argument('--cases_output', default=None,
                   help='Case-ID lists for the pairwise step (written if eda.pairwise)')
    p.add_argument('--covariates_output', default=None,
                   help='Discretized covariate table used for the metrics')
    args = p.parse_args()
    run_eda_stats(args.inputs, args.eda_config, args.output, args.covariates,
                  args.universe, args.cases_output, args.covariates_output)


if __name__ == '__main__':
    main()
