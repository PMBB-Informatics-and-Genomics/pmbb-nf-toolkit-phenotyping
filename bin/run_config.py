#!/usr/bin/env python3
"""
run_config.py — Run configuration page for the EDA report.

Usage:
    run_config.py --run_info run_info.json --biobank_config config.yaml \
                  --eda_config eda_config.json --output_dir eda

Writes to --output_dir:

    run_config.html      run parameters + resolved global / data_type / eda settings
    run_config.json      the same, machine-readable
    biobank_config.yaml  verbatim copy of the biobank YAML this run used
    biobank_config.html  the same YAML, viewable in a browser

run_info.json comes from main.nf: {section: {label: value}} as logged at start,
plus start_epoch_ms. Settings tables show the effective value of every setting
and where it came from (default / global / data_types / yaml). Phenotype-level
overrides are left out on purpose: they would grow with the phenotype count,
and configs/<data_type>/<phenotype>.json already holds them.
"""

import argparse
import hashlib
import html
import json
import os
import re
import shutil
import time
from collections import Counter

import pandas as pd
import yaml

from eda_report import _CSS, _JS, _table
from resolve_config import (ADVANCED_CODED_DEFAULTS, CAT_DEFAULTS, CODED_DEFAULTS,
                            GLOBAL_DEFAULTS, QUANT_DEFAULTS)

TYPE_DEFAULTS = {
    'quantitative': QUANT_DEFAULTS,
    'categorical': CAT_DEFAULTS,
    'coded': CODED_DEFAULTS,
    'advanced_coded': ADVANCED_CODED_DEFAULTS,
}


def _fmt(v):
    if isinstance(v, str):
        return v
    return json.dumps(v, default=str)


def _elapsed(start_ms, now_ms):
    secs = max(int((now_ms - start_ms) / 1000), 0)
    h, rem = divmod(secs, 3600)
    m, s = divmod(rem, 60)
    return f'{h}h {m:02d}m {s:02d}s' if h else f'{m}m {s:02d}s'


def global_rows(global_block):
    """Built-in global defaults plus any extra key set under global:."""
    keys = list(GLOBAL_DEFAULTS) + [k for k in global_block if k not in GLOBAL_DEFAULTS]
    return [{'setting': k,
             'value': _fmt(global_block.get(k, GLOBAL_DEFAULTS.get(k))),
             'source': 'global' if k in global_block else 'default'}
            for k in keys]


def data_type_rows(data_type, global_block, type_block):
    """
    Effective value of each data_type setting before phenotype overrides, in
    resolve_config.merge_params order: built-in → global → data_types.<type>.
    Global-only keys (sample_id_col, …) are in the global table instead.
    """
    defaults = TYPE_DEFAULTS[data_type]
    keys = list(defaults) + [k for k in type_block if k not in defaults]
    rows = []
    for k in keys:
        if k in type_block:
            value, source = type_block[k], 'data_types'
        elif k in global_block:
            value, source = global_block[k], 'global'
        else:
            value, source = defaults.get(k), 'default'
        rows.append({'setting': k, 'value': _fmt(value), 'source': source})
    return rows


def eda_rows(eda_block, settings):
    """Resolved eda settings (paths absolute, 'auto' expanded) with their source."""
    return [{'setting': k, 'value': _fmt(v), 'source': 'yaml' if k in eda_block else 'default'}
            for k, v in settings.items()]


def build_run_config(run_info, config, eda_config, now_ms=None):
    run_info = dict(run_info)
    start_ms = run_info.pop('start_epoch_ms', None)
    now_ms = int(time.time() * 1000) if now_ms is None else now_ms
    sections = {name: dict(rows) for name, rows in run_info.items()}
    if start_ms is not None:
        sections.setdefault('RUN', {})['elapsed at report build'] = _elapsed(start_ms, now_ms)

    global_block = config.get('global') or {}
    types_block = config.get('data_types') or {}
    counts = Counter((p or {}).get('data_type') or global_block.get('data_type')
                     for p in (config.get('phenotypes') or {}).values())
    data_types = {
        t: {'n_phenotypes': counts[t],
            'settings': data_type_rows(t, global_block, types_block.get(t) or {})}
        for t in TYPE_DEFAULTS if counts[t]
    }
    return {
        'run': sections,
        'global': global_rows(global_block),
        'data_types': data_types,
        'eda': eda_rows(config.get('eda') or {}, eda_config.get('settings') or {}),
    }


def _kv_table(rows, table_id):
    return _table(pd.DataFrame(list(rows.items()), columns=['parameter', 'value']), table_id)


def _split_comment(line):
    """Split a YAML line at its first # that starts a comment (not inside quotes)."""
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in '\'"':
            quote = ch
        elif ch == '#' and (i == 0 or line[i - 1].isspace()):
            return line[:i], line[i:]
    return line, ''


_KEY = re.compile(r'^(\s*(?:-\s+)?)([^\s:#\'"][^:#]*?|"[^"]*"|\'[^\']*\')(:)(?=\s|$)')


def _yaml_line_html(line):
    code, comment = _split_comment(line)
    m = _KEY.match(code)
    if m:
        out = (html.escape(m.group(1)) + f'<span class="k">{html.escape(m.group(2))}</span>:'
               + html.escape(code[m.end():]))
    else:
        out = html.escape(code)
    if comment:
        out += f'<span class="c">{html.escape(comment)}</span>'
    return out


def render_yaml_html(text, source_path, sha256):
    lines = text.splitlines()
    # Jump links to top-level blocks (global:, data_types:, phenotypes:, eda:, …)
    nav = [(m.group(1), i) for i, l in enumerate(lines, 1)
           if (m := re.match(r'^([A-Za-z_][\w-]*):', l))]
    nav_html = ' · '.join(f'<a href="#L{i}">{html.escape(k)}</a>' for k, i in nav)
    body = ''.join(f'<span class="l" id="L{i}"><a class="n" href="#L{i}">{i}</a>'
                   f'{_yaml_line_html(l)}\n</span>' for i, l in enumerate(lines, 1))
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Biobank Config</title><style>{_CSS}
pre.yaml{{background:var(--surface-2);border:1px solid var(--border);border-radius:8px;padding:12px 0;
overflow-x:auto;font:12.5px/1.55 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;margin:12px 0}}
pre.yaml .l{{display:block;padding:0 12px 0 0}}pre.yaml .l:target{{background:var(--border)}}
pre.yaml .n{{display:inline-block;width:4.5em;padding-right:1em;text-align:right;color:var(--muted);
text-decoration:none;user-select:none}}
pre.yaml .k{{color:var(--accent)}}pre.yaml .c{{color:var(--muted);font-style:italic}}</style>
</head><body><main>
<p class="small"><a href="eda_report.html">← EDA report</a> · <a href="run_config.html">Run configuration</a></p>
<h1>Biobank config</h1>
<p class="muted">Verbatim copy of the YAML this run used (<code>{html.escape(source_path)}</code>,
sha256 <code>{sha256[:12]}</code>). Raw file: <a href="biobank_config.yaml">biobank_config.yaml</a>.</p>
<p class="small">{nav_html}</p>
<pre class="yaml">{body}</pre>
</main></body></html>"""


def render_html(rc):
    run = rc['run']
    parts = []
    for i, (name, rows) in enumerate(run.items()):
        parts.append(f'<h2>{html.escape(name.capitalize())}</h2>{_kv_table(rows, f"run{i}")}')
    parts.append('<h2>Global settings</h2>'
                 + _table(pd.DataFrame(rc['global']), 'global'))
    parts.append('<h2>Data type settings</h2><p class="muted small">Before phenotype-level '
                 'overrides; see configs/&lt;data_type&gt;/&lt;phenotype&gt;.json for the '
                 'settings each phenotype ran with.</p>')
    if not rc['data_types']:
        parts.append('<p class="muted">No phenotypes.</p>')
    for t, d in rc['data_types'].items():
        n = d['n_phenotypes']
        parts.append(f'<h3>{html.escape(t)} <span class="muted small">({n:,} phenotype'
                     f'{"s" if n != 1 else ""})</span></h3>'
                     + _table(pd.DataFrame(d['settings']), f'dt_{t}'))
    parts.append('<h2>EDA settings</h2>' + _table(pd.DataFrame(rc['eda']), 'eda'))

    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>EDA Run Configuration</title><style>{_CSS}
h3{{font-size:15px;margin:20px 0 6px}}</style></head><body><main>
<p class="small"><a href="eda_report.html">← EDA report</a> · <a href="biobank_config.html">Biobank config (YAML)</a></p>
<h1>Run configuration</h1>
<p class="muted">Parameters and resolved biobank YAML settings for this run. Source: where
the effective value came from (default = built-in). Full timings and per-task resource use
are in pipeline_report.html and pipeline_timeline.html in the pipeline output directory.</p>
{''.join(parts)}
</main><script>{_JS}</script></body></html>"""


def run_config(run_info_path, config_path, eda_config_path, output_dir, now_ms=None):
    with open(run_info_path) as f:
        run_info = json.load(f)
    with open(config_path) as f:
        config = yaml.safe_load(f) or {}
    with open(eda_config_path) as f:
        eda_config = json.load(f)
    rc = build_run_config(run_info, config, eda_config, now_ms)
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'run_config.json'), 'w') as f:
        json.dump(rc, f, indent=2, default=str)
    with open(os.path.join(output_dir, 'run_config.html'), 'w') as f:
        f.write(render_html(rc))

    # copyfile follows the staged symlink, so the zip holds the real file
    shutil.copyfile(config_path, os.path.join(output_dir, 'biobank_config.yaml'))
    with open(config_path, 'rb') as f:
        raw = f.read()
    source = (run_info.get('WORKFLOW PARAMETERS') or {}).get('--biobank_config') or config_path
    with open(os.path.join(output_dir, 'biobank_config.html'), 'w') as f:
        f.write(render_yaml_html(raw.decode('utf-8', errors='replace'), str(source),
                                 hashlib.sha256(raw).hexdigest()))
    return rc


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run_info', required=True, help='Run parameters JSON from main.nf')
    p.add_argument('--biobank_config', required=True, help='Biobank YAML as passed to the run')
    p.add_argument('--eda_config', required=True, help='eda_config.json from RESOLVE_CONFIG')
    p.add_argument('--output_dir', default='.')
    args = p.parse_args()
    run_config(args.run_info, args.biobank_config, args.eda_config, args.output_dir)


if __name__ == '__main__':
    main()
