#!/usr/bin/env python3
"""
eda_report.py — Flags, profile selection and the HTML triage page.

Usage:
    eda_report.py --stats eda_stats.tsv --manifest phenotype_manifest.tsv \
                  --eda_config eda_config.json --output_dir eda \
                  [--cases a.cases.tsv b.cases.tsv] \
                  [--manifest_output phenotype_manifest.tsv]

Reads the per-column metrics from eda_stats.py (all columns, all phenotypes)
and writes to --output_dir:

    eda_summary.tsv              every column, every metric, plus `flags`
    eda_flags.tsv                long format: column_name, flag, value
    eda_report.html              self-contained triage page (aggregates only)
    profile_selection.tsv        columns chosen for ydata profiles + group
    restriction_suggestions.tsv  restrict_to entries + data-driven separation (suggestions only)
    pairwise_overlap.tsv         near-duplicate / nested binary pairs (eda.pairwise)
    filter_suggested.yaml        editable filter config for filter_phenotypes.py

and, with --manifest_output, a copy of the phenotype manifest with a `flags`
column added. Nothing here changes phenotype data: filtering is a separate,
opt-in step (filter_phenotypes.py) that writes new files.

Profile selection is bounded: opted-in columns (per-phenotype `eda` cascade)
plus the profile_flagged_top_k most severe flagged columns, capped at
max_profile_columns, split into groups of profile_group_size per phenotype.
"""

import argparse
import fnmatch
import html
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

from eda_stats import INTERACTION_SEP, LEVEL_SEP, parse_sep_detail

# Higher = worse; drives sorting, top-k and the profile cap.
FLAG_SEVERITY = {
    'separation_complete': 6,
    'separation_quasi': 5,
    'empty': 4,
    'no_controls': 4,
    'low_cases': 4,
    'high_missing': 3,
    'missing_by_covariate': 2,
    'pairwise_duplicate': 2,
    'pairwise_nested': 1,
    'zero_variance': 1,
    'single_level': 1,
}

FLAG_HELP = {
    'separation_complete': 'All cases (or controls) fall in one covariate level (or a small subset of levels)',
    'separation_quasi': 'At least separation_threshold of cases (or controls) fall in one covariate level (or a small subset)',
    'empty': 'No non-missing values',
    'no_controls': 'Binary column with cases but no controls (missing_as_control off / no roster?)',
    'low_cases': 'Binary column with fewer than min_cases cases',
    'high_missing': 'Missing fraction above max_missing',
    'missing_by_covariate': 'Missingness differs strongly between covariate levels',
    'pairwise_duplicate': 'Case set nearly identical to another column',
    'pairwise_nested': 'Case set (nearly) contained in another column',
    'zero_variance': 'Quantitative column with a single value',
    'single_level': 'Categorical column with a single level',
}


# ── Flags ──────────────────────────────────────────────────────────────────────

def compute_flags(stats, settings):
    """Return long-format flags DataFrame (column_name, source_file, flag, value)."""
    out = []

    def add(mask, flag, values):
        sub = stats.loc[mask]
        vals = values.loc[mask] if isinstance(values, pd.Series) else [values] * len(sub)
        for (_, r), v in zip(sub.iterrows(), vals):
            out.append((r['column_name'], r['source_file'], flag, v))

    binary = stats['kind'] == 'binary'
    n = stats['n'].fillna(0)
    cases = stats['n_cases'].fillna(0)
    controls = stats['n_controls'].fillna(0)

    # Cases but no controls: non-cases are NA (missing_as_control off or no
    # roster), so "missingness" is just the case pattern — don't flag it twice.
    no_controls = binary & (n > 0) & (controls == 0) & (cases > 0)

    add(n == 0, 'empty', 0)
    add(no_controls, 'no_controls', stats['n_cases'])
    add(binary & (n > 0) & (controls > 0) & (cases < settings['min_cases']),
        'low_cases', stats['n_cases'])
    add((stats['missing_frac'] > settings['max_missing']) & ~no_controls, 'high_missing',
        stats['missing_frac'].round(4))
    add(stats['sep_type'] == 'complete', 'separation_complete', stats['sep_detail'])
    add(stats['sep_type'] == 'quasi', 'separation_quasi', stats['sep_detail'])
    add((stats['miss_gap'] > settings['missing_gap']) & ~no_controls, 'missing_by_covariate',
        stats['miss_gap_covariate'].astype(str) + ': ' + stats['miss_gap_detail'].astype(str))
    add((stats['kind'] == 'quantitative') & (n > 0) & (stats['n_unique'] <= 1),
        'zero_variance', stats['n_unique'])
    add((stats['kind'] == 'categorical') & (n > 0) & (stats['n_unique'] <= 1),
        'single_level', stats['n_unique'])

    return pd.DataFrame(out, columns=['column_name', 'source_file', 'flag', 'value'])


def attach_flags(stats, flags):
    """Add `flags` (comma list) and `severity` (max) columns to stats."""
    key = ['column_name', 'source_file']
    if flags.empty:
        stats['flags'] = ''
        stats['severity'] = 0
        return stats
    f = flags.assign(sev=flags['flag'].map(FLAG_SEVERITY).fillna(0))
    agg = f.groupby(key).agg(
        flags=('flag', lambda s: ','.join(sorted(set(s), key=lambda x: -FLAG_SEVERITY[x]))),
        severity=('sev', 'max'),
    ).reset_index()
    stats = stats.drop(columns=[c for c in ('flags', 'severity') if c in stats], errors='ignore')
    stats = stats.merge(agg, on=key, how='left')
    stats['flags'] = stats['flags'].fillna('')
    stats['severity'] = stats['severity'].fillna(0).astype(int)
    return stats


# ── Profile selection ─────────────────────────────────────────────────────────

def opted_in(row, eda_cfg):
    """True if the per-phenotype `eda` cascade opts this column into profiling."""
    pheno = eda_cfg.get('phenotypes', {}).get(row['phenotype'])
    if not pheno:
        return False
    eda = pheno.get('eda')
    if isinstance(eda, bool):
        return eda
    if isinstance(eda, list):
        targets = [str(row.get('code') or ''), str(row['column_name'])]
        return any(fnmatch.fnmatch(t, str(p)) for p in eda for t in targets if t)
    return False


def select_profiles(stats, eda_cfg):
    """Return profile_selection DataFrame (column_name, source_file, phenotype, group, reason)."""
    settings = eda_cfg['settings']
    stats = stats.copy()
    stats['opt_in'] = stats.apply(lambda r: opted_in(r, eda_cfg), axis=1)

    flagged = stats[(stats['severity'] > 0) & ~stats['opt_in']]
    flagged = flagged.sort_values(['severity', 'column_name'], ascending=[False, True])
    top = flagged.head(int(settings['profile_flagged_top_k'])).index

    chosen = stats.loc[stats['opt_in'] | stats.index.isin(top)].copy()
    chosen['reason'] = np.where(chosen['opt_in'], 'opt_in', 'flagged')
    chosen = chosen.sort_values(['severity', 'opt_in', 'column_name'],
                                ascending=[False, False, True])

    cap = int(settings['max_profile_columns'])
    if len(chosen) > cap:
        print(f'WARNING: {len(chosen)} columns selected for profiling; keeping the '
              f'{cap} most severe (eda.max_profile_columns)', file=sys.stderr)
        chosen = chosen.head(cap)

    size = max(1, int(settings['profile_group_size']))
    groups = []
    for pheno, sub in chosen.sort_values(['phenotype', 'column_name']).groupby('phenotype', sort=True):
        n_chunks = math.ceil(len(sub) / size)
        safe = re.sub(r'[^A-Za-z0-9_.-]', '_', str(pheno))
        for i in range(n_chunks):
            part = sub.iloc[i * size:(i + 1) * size].copy()
            part['group'] = safe if n_chunks == 1 else f'{safe}_{i + 1:03d}'
            groups.append(part)
    cols = ['column_name', 'source_file', 'phenotype', 'group', 'reason']
    if not groups:
        return pd.DataFrame(columns=cols)
    return pd.concat(groups)[cols].reset_index(drop=True)


# ── Restriction suggestions ───────────────────────────────────────────────────

def _level_code(dictionaries, covariate, label):
    """Map a level label to its encoded value via the covariate's dictionary."""
    d = dictionaries.get(covariate) or {}
    for k, v in d.items():
        if str(k).lower() == str(label).lower():
            return _fmt(v)
    return str(label)


def _level_label(dictionaries, covariate, code):
    d = dictionaries.get(covariate) or {}
    for k, v in d.items():
        if _fmt(v) == str(code):
            return str(k)
    return str(code)


def _labels(dictionaries, covariate, level):
    """Readable label for a sep_level: subsets (|) and interactions (&) mapped part by part."""
    parts = str(covariate).split(INTERACTION_SEP)
    return LEVEL_SEP.join(
        INTERACTION_SEP.join(_level_label(dictionaries, c, l)
                             for c, l in zip(parts, combo.split(INTERACTION_SEP)))
        for combo in str(level).split(LEVEL_SEP))


def _fmt(v):
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def restriction_suggestions(stats, eda_cfg):
    """
    One row per (column, covariate, level) from two evidence sources:
    yaml (restrict_to) and data (cases confined to a level, or to a subset of
    levels, of any covariate — every entry of sep_detail, not just the top one).
    """
    dictionaries = eda_cfg.get('covariate_dictionaries') or {}
    phenos = eda_cfg.get('phenotypes') or {}
    rows = []
    candidates = stats[~stats['is_covariate'].astype(bool)]

    for _, r in candidates.iterrows():
        p = phenos.get(r['phenotype']) or {}
        rt = p.get('restrict_to')
        if rt:
            cov, label = next(iter(rt.items()))
            rows.append((r, cov, _level_code(dictionaries, cov, label), 'yaml',
                         f'restrict_to: {{{cov}: {label}}}', None))
        for s in parse_sep_detail(r.get('sep_detail')):
            if s['outcome'] == 'cases':
                rows.append((r, s['covariate'], LEVEL_SEP.join(s['levels']), 'data',
                             f"{s['sep_type']}: {s['share']:.3f} of cases", s['share']))

    cols = ['column_name', 'source_file', 'phenotype', 'covariate', 'level', 'level_label',
            'sources', 'n_sources', 'share', 'evidence', 'n_cases']
    if not rows:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame([{
        'column_name': r['column_name'], 'source_file': r['source_file'],
        'phenotype': r['phenotype'], 'covariate': cov, 'level': str(level),
        'source': src, 'evidence': ev, 'share': share, 'n_cases': r.get('n_cases'),
    } for r, cov, level, src, ev, share in rows])
    out = df.groupby(['column_name', 'source_file', 'phenotype', 'covariate', 'level'],
                     sort=True).agg(
        sources=('source', lambda s: ','.join(sorted(set(s)))),
        n_sources=('source', 'nunique'),
        share=('share', 'max'),
        evidence=('evidence', lambda s: ' | '.join(dict.fromkeys(s))),
        n_cases=('n_cases', 'first'),
    ).reset_index()
    out['level_label'] = [_labels(dictionaries, c, level)
                          for c, level in zip(out['covariate'], out['level'])]
    return out.sort_values(['n_sources', 'column_name'], ascending=[False, True])[cols]


# ── Pairwise overlap ──────────────────────────────────────────────────────────

def pairwise_overlap(case_paths, stats, settings):
    """Near-duplicate / nested binary columns from sparse case co-occurrence."""
    cols = ['column_a', 'column_b', 'n_cases_a', 'n_cases_b', 'n_both',
            'jaccard', 'containment', 'relation']
    names, sources, case_lists = [], [], []
    for path in case_paths or []:
        df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        for _, r in df.iterrows():
            ids = [x for x in r['case_ids'].split(',') if x]
            if ids:
                names.append(r['column_name'])
                sources.append(r['source_file'])
                case_lists.append(ids)
    if len(names) < 2:
        return pd.DataFrame(columns=cols)

    from scipy import sparse

    index = {}
    rows, colidx = [], []
    for j, ids in enumerate(case_lists):
        for sid in ids:
            rows.append(index.setdefault(sid, len(index)))
            colidx.append(j)
    X = sparse.csc_matrix((np.ones(len(rows), dtype=np.int32), (rows, colidx)),
                          shape=(len(index), len(names)))
    co = (X.T @ X).tocoo()
    size = np.asarray(X.sum(axis=0)).ravel()

    pheno_of = dict(zip(zip(stats['column_name'], stats['source_file']), stats['phenotype']))
    within = settings.get('pairwise_scope') == 'within_phenotype'
    out = []
    for i, j, both in zip(co.row, co.col, co.data):
        if i >= j:
            continue
        if within and pheno_of.get((names[i], sources[i])) != pheno_of.get((names[j], sources[j])):
            continue
        a, b = size[i], size[j]
        jac = both / (a + b - both)
        cont = both / min(a, b)
        if jac >= settings['pairwise_jaccard']:
            rel = 'duplicate'
        elif cont >= settings['pairwise_containment']:
            rel = 'nested'
        else:
            continue
        # nested: column_a is the smaller set, contained in column_b
        if rel == 'nested' and a > b:
            i, j, a, b = j, i, b, a
        out.append((names[i], names[j], int(a), int(b), int(both), jac, cont, rel))
    df = pd.DataFrame(out, columns=cols)
    df = df.sort_values(['jaccard', 'containment'], ascending=False)
    return df.head(int(settings['pairwise_top_k'])).reset_index(drop=True)


def pairwise_flags(pairs, stats):
    """Flags for columns in pairwise relations (value = partner column)."""
    src = dict(zip(stats['column_name'], stats['source_file']))
    out = []
    for _, p in pairs.iterrows():
        flag = 'pairwise_duplicate' if p['relation'] == 'duplicate' else 'pairwise_nested'
        out.append((p['column_a'], src.get(p['column_a']), flag, p['column_b']))
        if flag == 'pairwise_duplicate':
            out.append((p['column_b'], src.get(p['column_b']), flag, p['column_a']))
    return pd.DataFrame(out, columns=['column_name', 'source_file', 'flag', 'value'])


# ── HTML ──────────────────────────────────────────────────────────────────────

def _histogram_svg(values, title, xlabel, bins, log_note=''):
    """Inline-SVG histogram with per-bar hover titles (single series)."""
    values = pd.Series(values).dropna()
    if values.empty:
        return f'<figure class="chart"><figcaption>{html.escape(title)}</figcaption><p class="muted">No data.</p></figure>'
    counts, edges = np.histogram(values, bins=bins)
    W, H, pl, pb, pt = 520, 180, 44, 34, 10
    pw, ph = W - pl - 10, H - pb - pt
    ymax = max(int(counts.max()), 1)
    bw = pw / len(counts)
    bars = []
    for k, c in enumerate(counts):
        h = ph * c / ymax
        x = pl + k * bw + 1
        y = pt + ph - h
        w = max(bw - 2, 1)
        tip = f'{edges[k]:.3g} – {edges[k + 1]:.3g}: {c:,} columns'
        r = min(4, w / 2, h)
        path = (f'M{x:.1f},{pt + ph:.1f} V{y + r:.1f} Q{x:.1f},{y:.1f} {x + r:.1f},{y:.1f} '
                f'H{x + w - r:.1f} Q{x + w:.1f},{y:.1f} {x + w:.1f},{y + r:.1f} V{pt + ph:.1f} Z')
        bars.append(f'<g class="bar"><rect x="{pl + k * bw:.1f}" y="{pt}" width="{bw:.1f}" '
                    f'height="{ph}" class="hit"/>'
                    + (f'<path d="{path}" class="mark"/>' if c else '')
                    + f'<title>{html.escape(tip)}</title></g>')
    ticks = []
    for frac in (0, 0.5, 1):
        v = ymax * frac
        y = pt + ph - ph * frac
        ticks.append(f'<line x1="{pl}" x2="{W - 10}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/>'
                     f'<text x="{pl - 6}" y="{y + 4:.1f}" class="tick" text-anchor="end">{v:,.0f}</text>')
    xt = []
    for k in (0, len(counts) // 2, len(counts)):
        x = pl + k * bw
        xt.append(f'<text x="{x:.1f}" y="{H - 16}" class="tick" text-anchor="middle">{edges[k]:.3g}</text>')
    return (f'<figure class="chart"><figcaption>{html.escape(title)}</figcaption>'
            f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}">'
            + ''.join(ticks) + ''.join(bars) + ''.join(xt)
            + f'<text x="{pl + pw / 2:.1f}" y="{H - 2}" class="tick" text-anchor="middle">'
            f'{html.escape(xlabel)}</text></svg>'
            + (f'<p class="muted small">{html.escape(log_note)}</p>' if log_note else '')
            + '</figure>')


def _table(df, table_id, max_rows=None):
    if df.empty:
        return '<p class="muted">None.</p>'
    note = ''
    if max_rows and len(df) > max_rows:
        note = (f'<p class="muted small">Showing {max_rows:,} of {len(df):,} rows — '
                f'see the TSV for all.</p>')
        df = df.head(max_rows)
    head = ''.join(f'<th data-k="{i}">{html.escape(str(c))}</th>' for i, c in enumerate(df.columns))
    body = []
    for row in df.itertuples(index=False):
        cells = []
        for v in row:
            if isinstance(v, float):
                txt = '' if math.isnan(v) else f'{v:.4g}'
                cells.append(f'<td class="num">{txt}</td>')
            elif isinstance(v, (int, np.integer)):
                cells.append(f'<td class="num">{v:,}</td>')
            else:
                cells.append(f'<td>{html.escape("" if v is None else str(v))}</td>')
        body.append('<tr>' + ''.join(cells) + '</tr>')
    return (f'<div class="tablewrap"><table id="{table_id}" class="sortable"><thead><tr>{head}'
            f'</tr></thead><tbody>{"".join(body)}</tbody></table></div>{note}')


_CSS = """
:root{--surface:#fcfcfb;--surface-2:#f3f2ef;--text:#0b0b0b;--text-2:#52514e;--muted:#7a7974;
--border:#e2e1dc;--series-1:#2a78d6;--accent:#2a78d6;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--surface:#1a1a19;--surface-2:#242422;
--text:#f0efec;--text-2:#c3c2b7;--muted:#8f8e86;--border:#3a3936;--series-1:#3987e5;--accent:#6aa7ee;color-scheme:dark}}
:root[data-theme="dark"]{--surface:#1a1a19;--surface-2:#242422;--text:#f0efec;--text-2:#c3c2b7;
--muted:#8f8e86;--border:#3a3936;--series-1:#3987e5;--accent:#6aa7ee;color-scheme:dark}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--text);font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1200px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 8px;border-bottom:1px solid var(--border);padding-bottom:4px}
.muted{color:var(--muted)}.small{font-size:12px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:16px 0}
.tile{background:var(--surface-2);border-radius:8px;padding:12px}.tile b{display:block;font-size:22px}
.tile span{color:var(--text-2);font-size:12px}
.charts{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:16px}
.chart{margin:0}.chart figcaption{font-weight:600;margin-bottom:4px}.chart svg{width:100%;height:auto}
.chart .mark{fill:var(--series-1)}.chart .hit{fill:transparent}.chart .bar:hover .mark{opacity:.8}
.chart .grid{stroke:var(--border);stroke-width:1}.chart .tick{fill:var(--muted);font-size:11px}
.tablewrap{overflow-x:auto;border:1px solid var(--border);border-radius:8px;max-height:640px;overflow-y:auto}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:5px 8px;border-bottom:1px solid var(--border);text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:var(--surface-2);cursor:pointer;user-select:none;color:var(--text-2)}
td.num{text-align:right;font-variant-numeric:tabular-nums}
input.filter{width:100%;max-width:420px;padding:6px 10px;margin:4px 0 8px;border:1px solid var(--border);
border-radius:6px;background:var(--surface);color:var(--text);font:inherit}
a{color:var(--accent)}ul.links{columns:3 220px}
"""

_JS = """
document.querySelectorAll('table.sortable').forEach(function(t){
  t.querySelectorAll('th').forEach(function(th,i){th.addEventListener('click',function(){
    var tb=t.tBodies[0],rows=Array.from(tb.rows),asc=th.dataset.asc!=='1';
    t.querySelectorAll('th').forEach(function(h){delete h.dataset.asc});th.dataset.asc=asc?'1':'0';
    rows.sort(function(a,b){var x=a.cells[i].textContent.replace(/,/g,''),y=b.cells[i].textContent.replace(/,/g,'');
      var nx=parseFloat(x),ny=parseFloat(y);var r=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);return asc?r:-r});
    rows.forEach(function(r){tb.appendChild(r)});});});});
document.querySelectorAll('input.filter').forEach(function(inp){
  var t=document.getElementById(inp.dataset.table);if(!t)return;
  inp.addEventListener('input',function(){var q=inp.value.toLowerCase();
    Array.from(t.tBodies[0].rows).forEach(function(r){r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?'':'none'});});});
"""


def render_html(stats, flags, suggestions, pairs, selection, settings):
    n_cols = len(stats)
    flagged = stats[stats['severity'] > 0].sort_values(
        ['severity', 'column_name'], ascending=[False, True])
    flag_counts = flags.groupby('flag')['column_name'].nunique() if not flags.empty else pd.Series(dtype=int)

    tiles = [('Columns', f'{n_cols:,}'), ('Flagged', f'{len(flagged):,}')]
    for kind in ('binary', 'quantitative', 'categorical'):
        tiles.append((kind.capitalize(), f'{int((stats["kind"] == kind).sum()):,}'))
    tiles.append(('Profiled', f'{len(selection):,}'))
    tiles_html = ''.join(f'<div class="tile"><b>{v}</b><span>{html.escape(k)}</span></div>'
                         for k, v in tiles)

    flag_rows = pd.DataFrame({
        'flag': list(FLAG_SEVERITY),
        'columns': [int(flag_counts.get(f, 0)) for f in FLAG_SEVERITY],
        'meaning': [FLAG_HELP[f] for f in FLAG_SEVERITY],
    })
    flag_rows = flag_rows[flag_rows['columns'] > 0]

    binary = stats[stats['kind'] == 'binary']
    charts = (
        _histogram_svg(stats['missing_frac'], 'Missing fraction, all columns',
                       'missing fraction', bins=np.linspace(0, 1, 21))
        + _histogram_svg(np.log10(binary['n_cases'].clip(lower=1)),
                         'Cases per binary column (log10)', 'log10(n_cases)', bins=20,
                         log_note=f'low_cases threshold: min_cases = {settings["min_cases"]} '
                         f'(log10 = {math.log10(max(settings["min_cases"], 1)):.2f})'
                         if len(binary) else '')
    )

    show = ['column_name', 'phenotype', 'kind', 'flags', 'n', 'missing_frac',
            'n_cases', 'case_frac', 'sep_detail', 'miss_gap_covariate', 'miss_gap']
    flagged_tbl = _table(flagged[show], 'flagged', settings['flagged_table_max_rows'])

    sugg_tbl = _table(suggestions.drop(columns=['source_file'], errors='ignore'), 'sugg',
                      settings['flagged_table_max_rows'])
    pair_tbl = _table(pairs, 'pairs', 500) if pairs is not None else \
        '<p class="muted">Pairwise checks are off (eda.pairwise: false).</p>'

    groups = sorted(selection['group'].unique()) if len(selection) else []
    if groups:
        links = '<ul class="links">' + ''.join(
            f'<li><a href="profiles/{html.escape(g)}.html">{html.escape(g)}</a> '
            f'<span class="muted small">({int((selection["group"] == g).sum())})</span></li>'
            for g in groups) + '</ul>'
    else:
        links = '<p class="muted">No columns selected for profiling.</p>'

    def filt(tid):
        return f'<input class="filter" data-table="{tid}" placeholder="Filter rows…" aria-label="Filter">'

    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Phenotype EDA Report</title><style>{_CSS}</style></head><body><main>
<h1>Phenotype EDA report</h1>
<p class="muted">Metrics for every output column; ydata profiles for a bounded subset.
Thresholds: max_missing {settings['max_missing']}, min_cases {settings['min_cases']},
separation_threshold {settings['separation_threshold']}, separation_min_cases {settings['separation_min_cases']},
max_subset_levels {settings['max_subset_levels']},
min_expected {settings['min_expected']},
missing_gap {settings['missing_gap']}. Covariates: {html.escape(', '.join(map(str, (settings.get('covariates') or []) + (settings.get('covariates_file_cols') or []) + [INTERACTION_SEP.join(c) for c in settings.get('separation_interactions') or []])) or ('file' if settings.get('covariates_file') else 'none'))}.</p>
<div class="tiles">{tiles_html}</div>
<h2>Flags</h2>{_table(flag_rows, 'flagcounts')}
<h2>Distributions</h2><div class="charts">{charts}</div>
<h2>Flagged columns</h2>{filt('flagged')}{flagged_tbl}
<h2>Restriction suggestions (sex / cohort-specific)</h2>
<p class="muted small">Report only. Sources: yaml = restrict_to, data = cases confined to a level or subset of levels (share = fraction of cases inside).</p>
{filt('sugg')}{sugg_tbl}
<h2>Pairwise overlap</h2>{pair_tbl}
<h2>Profiles</h2>{links}
</main><script>{_JS}</script></body></html>"""


# ── Flagged manifest and suggested filter ─────────────────────────────────────

def flagged_manifest(manifest, stats):
    """The manifest with a `flags` column (comma list; '' when unflagged)."""
    flags = stats.drop_duplicates('column_name').set_index('column_name')['flags']
    out = manifest.copy()
    out['flags'] = out['column_name'].map(flags).fillna('')
    return out


# Flags whose columns the suggested filter drops outright, in section order.
SUGGEST_DROP = ['empty', 'zero_variance', 'single_level', 'low_cases', 'high_missing']


def _q(v):
    """A YAML-safe double-quoted scalar."""
    return json.dumps(str(v))


def suggest_filter_yaml(stats, suggestions, pairs, eda_cfg):
    """
    Text of an editable filter config (see filter_phenotypes.py).

    Active entries: drops for SUGGEST_DROP flags, and restrict_to restrictions.
    Everything that needs a human call is written commented out: a
    restrict_separated rule at the EDA's separation_threshold (with the columns it
    would restrict listed), pairwise duplicates and missing-by-covariate columns.
    """
    settings = eda_cfg['settings']
    dictionaries = eda_cfg.get('covariate_dictionaries') or {}
    st = stats[~stats['is_covariate']].drop_duplicates('column_name')
    flag_sets = {r.column_name: set(filter(None, str(r.flags).split(',')))
                 for r in st.itertuples()}
    info = st.set_index('column_name')

    def note(col):
        r = info.loc[col]
        bits = []
        if pd.notna(r.get('n_cases')):
            bits.append(f"{int(r['n_cases'])} cases")
        if pd.notna(r.get('missing_frac')):
            bits.append(f"missing {float(r['missing_frac']):.2f}")
        return ', '.join(bits)

    sugg = suggestions.copy()
    sugg['active'] = sugg['sources'].str.contains('yaml')

    lines = [
        '# Suggested phenotype filter — generated by eda_report.py from this run\'s EDA.',
        '# Edit freely, then run the filter step. Raw outputs are never changed; the',
        '# filter writes new files under {output_dir}/filtered/:',
        '#   nextflow run main.nf <same args> --run_gather --run_eda --filter_config this.yaml -resume',
        '#   filter_phenotypes.py --config this.yaml --input phenotype_table.wide.tsv \\',
        '#       --covariates eda/eda_covariates.tsv --summary eda/eda_summary.tsv --output_dir filtered',
        '#',
        '# drop_columns     : columns removed from the filtered table',
        '# keep_columns     : always kept — wins over drop_columns and drop_flags',
        '# drop_flags       : also drop every column carrying any of these EDA flags',
        '#                    (rule-based alternative to the explicit lists below)',
        '# restrictions     : column (glob ok) -> {covariate, level}; samples in any other',
        '#                    covariate level are set to NA, cases included',
        '# restrict_separated: {min_share, covariates} — restrict every column whose cases',
        '#                    are confined to a level / subset of levels (eda_summary.tsv)',
        '# restrict_unknown : na | keep — samples with no value for that covariate',
        f"# Thresholds used: min_cases {settings['min_cases']}, max_missing {settings['max_missing']}.",
    ]
    n_noctrl = sum('no_controls' in f for f in flag_sets.values())
    if n_noctrl:
        lines.append(f'# NOTE: {n_noctrl} binary columns have no controls (no --sample_list, or '
                     'missing_as_control: false) — rerun with a roster rather than dropping them.')
    lines += ['', 'restrict_unknown: na', 'drop_flags: []', 'keep_columns: []', '', 'drop_columns:']

    placed = set()
    for flag in SUGGEST_DROP:
        cols = sorted(c for c, f in flag_sets.items() if flag in f and c not in placed)
        if not cols:
            continue
        lines.append(f'  # {flag} — {FLAG_HELP[flag]}')
        lines += [f'  - {_q(c)}  # {note(c)}' for c in cols]
        placed.update(cols)

    review = []
    if pairs is not None and len(pairs):
        dup = pairs[pairs['relation'] == 'duplicate']
        if len(dup):
            review.append('  # review: pairwise duplicates — keep one of each pair')
            review += [f'  # - {_q(b)}  # duplicate of {a} (jaccard {j:.3f})'
                       for a, b, j in zip(dup['column_a'], dup['column_b'], dup['jaccard'])]
    miss_cols = sorted(c for c, f in flag_sets.items()
                       if 'missing_by_covariate' in f and c not in placed)
    if miss_cols:
        review.append('  # review: missingness differs by covariate (cohort-specific source?)')
        review += [f'  # - {_q(c)}  # {info.loc[c, "miss_gap_covariate"]}: '
                   f'{info.loc[c, "miss_gap_detail"]}' for c in miss_cols]
    if not placed and not review:
        lines[-1] = 'drop_columns: []'
    lines += review

    lines += ['', 'restrictions:']
    act = sugg[sugg['active']].drop_duplicates('column_name')
    for r in act.itertuples():
        lines.append(f'  {_q(r.column_name)}: {{covariate: {_q(r.covariate)}, level: {_q(r.level)}}}'
                     f'  # {r.level_label} · restrict_to')
    if not len(act):
        lines[-1] = 'restrictions: {}'

    # Data-driven: one rule instead of per-column entries; review, then uncomment.
    # (the filter skips a separated covariate a column already has a restriction on)
    explicit = set(zip(act['column_name'], act['covariate']))
    data = sugg[sugg['sources'].str.contains('data')
                & np.array([(c, v) not in explicit
                            for c, v in zip(sugg['column_name'], sugg['covariate'])], dtype=bool)]
    threshold = settings['separation_threshold']
    covs = (list(settings.get('covariates') or []) + list(settings.get('covariates_file_cols') or [])
            + [INTERACTION_SEP.join(c) for c in settings.get('separation_interactions') or []])
    lines += ['', '# restrict_separated: restricts every column whose cases are confined to a',
              '# level (or subset of levels) of one of `covariates` with share >= min_share,',
              '# read from eda_summary.tsv. Review — separation can be a real association.',
              f'# At the EDA threshold ({threshold}) it would restrict {data["column_name"].nunique()} column(s):']
    for r in data.itertuples():
        lines.append(f'#   {_q(r.column_name)}  {r.covariate} = {r.level} ({r.level_label})'
                     f'  share {float(r.share):.3f}')
    lines.append(f'# restrict_separated: {{min_share: {threshold}, '
                 f'covariates: [{", ".join(_q(c) for c in covs)}]}}')
    return '\n'.join(lines) + '\n'


# ── Driver ────────────────────────────────────────────────────────────────────

def run_eda_report(stats_path, manifest_path, eda_config_path, output_dir, case_paths=None,
                   manifest_output=None):
    with open(eda_config_path) as f:
        eda_cfg = json.load(f)
    settings = eda_cfg['settings']
    os.makedirs(output_dir, exist_ok=True)

    stats = pd.read_csv(stats_path, sep='\t', dtype={'column_name': str, 'sep_level': str})
    manifest = None
    if manifest_path and os.path.isfile(manifest_path) and os.path.getsize(manifest_path):
        manifest = pd.read_csv(manifest_path, sep='\t', dtype=str, keep_default_na=False)
        mf = manifest.drop_duplicates('column_name')
        stats = stats.drop(columns=['code'], errors='ignore').merge(
            mf[['column_name', 'phenotype', 'code']].rename(columns={'phenotype': '_mp'}),
            on='column_name', how='left')
        stats['phenotype'] = stats['_mp'].where(stats['_mp'].notna() & (stats['_mp'] != ''),
                                                stats['phenotype'])
        stats = stats.drop(columns=['_mp'])
    else:
        stats['code'] = ''
    stats['is_covariate'] = stats['is_covariate'].astype(str).str.lower() == 'true'

    flags = compute_flags(stats, settings)
    pairs = None
    if settings.get('pairwise'):
        pairs = pairwise_overlap(case_paths, stats, settings)
        flags = pd.concat([flags, pairwise_flags(pairs, stats)], ignore_index=True)
        pairs.to_csv(os.path.join(output_dir, 'pairwise_overlap.tsv'), sep='\t', index=False)

    stats = attach_flags(stats, flags)
    selection = select_profiles(stats, eda_cfg)
    suggestions = restriction_suggestions(stats, eda_cfg)

    stats.sort_values(['severity', 'column_name'], ascending=[False, True]).to_csv(
        os.path.join(output_dir, 'eda_summary.tsv'), sep='\t', index=False)
    flags.to_csv(os.path.join(output_dir, 'eda_flags.tsv'), sep='\t', index=False)
    selection.to_csv(os.path.join(output_dir, 'profile_selection.tsv'), sep='\t', index=False)
    suggestions.to_csv(os.path.join(output_dir, 'restriction_suggestions.tsv'), sep='\t', index=False)
    with open(os.path.join(output_dir, 'eda_report.html'), 'w') as f:
        f.write(render_html(stats, flags, suggestions, pairs, selection, settings))
    with open(os.path.join(output_dir, 'filter_suggested.yaml'), 'w') as f:
        f.write(suggest_filter_yaml(stats, suggestions, pairs, eda_cfg))
    if manifest is not None and manifest_output:
        flagged_manifest(manifest, stats).to_csv(manifest_output, sep='\t', index=False)

    print(f'EDA: {len(stats)} columns, {int((stats["severity"] > 0).sum())} flagged, '
          f'{len(selection)} selected for profiling, {len(suggestions)} restriction '
          f'suggestions → {output_dir}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stats', required=True, help='Concatenated eda_stats.py output')
    p.add_argument('--manifest', default=None, help='phenotype_manifest.tsv')
    p.add_argument('--eda_config', required=True)
    p.add_argument('--output_dir', default='.')
    p.add_argument('--cases', nargs='*', default=[], help='*.cases.tsv from eda_stats.py')
    p.add_argument('--manifest_output', default=None,
                   help='Write the manifest with a `flags` column here')
    args = p.parse_args()
    run_eda_report(args.stats, args.manifest, args.eda_config, args.output_dir, args.cases,
                   args.manifest_output)


if __name__ == '__main__':
    main()
