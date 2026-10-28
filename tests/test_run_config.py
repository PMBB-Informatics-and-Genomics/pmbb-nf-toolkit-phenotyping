import json
import os
import sys

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from resolve_config import EDA_DEFAULTS, resolve
from run_config import _yaml_line_html, build_run_config, run_config

CONFIG = {
    'global': {'sample_id_col': 'id', 'min_occurrences': 2, 'outlier_method': 'zscore'},
    'data_types': {'quantitative': {'stats': ['mean'], 'outlier_method': 'none'}},
    'phenotypes': {
        'sex': {'table': 't.tsv', 'data_type': 'categorical', 'value_col': 'v',
                'dictionary': {'Female': 0, 'Male': 1}},
        'ldl': {'table': 't.tsv', 'data_type': 'quantitative', 'value_col': 'v'},
        'hdl': {'table': 't.tsv', 'data_type': 'quantitative', 'value_col': 'v',
                'stats': ['median']},
    },
    'eda': {'covariates': ['sex'], 'min_cases': 50},
}
RUN_INFO = {
    'RUN': {'profile': 'test', 'resumed': False},
    'WORKFLOW PARAMETERS': {'--run_eda': 'ENABLED'},
    'start_epoch_ms': 1_000_000,
}


def _rows(rows):
    return {r['setting']: (r['value'], r['source']) for r in rows}


def _eda_config(tmp_path):
    (tmp_path / 't.tsv').write_text('id\tv\nA\t1\n')
    p = tmp_path / 'c.yaml'
    p.write_text(yaml.safe_dump(CONFIG))
    resolve(str(p), str(tmp_path / 'out'))
    return p, tmp_path / 'out' / 'eda_config.json'


def test_sources_follow_merge_order(tmp_path):
    _, eda_path = _eda_config(tmp_path)
    rc = build_run_config(RUN_INFO, CONFIG, json.loads(eda_path.read_text()),
                          now_ms=1_000_000 + 3_725_000)

    assert rc['run']['RUN']['elapsed at report build'] == '1h 02m 05s'
    assert rc['run']['WORKFLOW PARAMETERS'] == {'--run_eda': 'ENABLED'}

    g = _rows(rc['global'])
    assert g['sample_id_col'] == ('id', 'global')
    assert g['reference_date'] == ('today', 'default')
    assert g['min_occurrences'] == ('2', 'global')  # extra global key is listed

    # Only data types that have phenotypes; phenotype overrides never shown
    assert set(rc['data_types']) == {'quantitative', 'categorical'}
    assert rc['data_types']['quantitative']['n_phenotypes'] == 2
    q = _rows(rc['data_types']['quantitative']['settings'])
    assert q['stats'] == ('["mean"]', 'data_types')          # not hdl's ["median"]
    assert q['outlier_method'] == ('none', 'data_types')     # data_types beats global
    assert q['min_occurrences'] == ('2', 'global')
    assert q['scale'] == ('false', 'default')
    assert 'sample_id_col' not in q

    e = _rows(rc['eda'])
    assert e['min_cases'] == ('50', 'yaml')
    assert e['covariates'] == ('["sex"]', 'yaml')
    assert e['min_level_n'] == (str(EDA_DEFAULTS['min_level_n']), 'default')


def test_writes_html_and_json(tmp_path):
    cfg_path, eda_path = _eda_config(tmp_path)
    info = tmp_path / 'run_info.json'
    info.write_text(json.dumps(RUN_INFO))
    out = tmp_path / 'eda'
    run_config(str(info), str(cfg_path), str(eda_path), str(out))

    page = (out / 'run_config.html').read_text()
    assert 'Global settings' in page and 'EDA settings' in page
    assert 'href="eda_report.html"' in page
    assert json.loads((out / 'run_config.json').read_text())['eda']

    # Verbatim copy + viewer, with jump links to the top-level blocks
    assert (out / 'biobank_config.yaml').read_text() == cfg_path.read_text()
    view = (out / 'biobank_config.html').read_text()
    assert 'href="biobank_config.yaml"' in view
    for block in ('global', 'data_types', 'phenotypes', 'eda'):
        assert f'>{block}</a>' in view


def test_yaml_line_highlighting():
    assert _yaml_line_html('  min_cases: 5  # low') == (
        '  <span class="k">min_cases</span>: 5  <span class="c"># low</span>')
    # # inside quotes or a value without a leading space is not a comment
    assert 'class="c"' not in _yaml_line_html('  value_col: "A # B"')
    assert 'class="c"' not in _yaml_line_html('  code: E11#x')
    assert _yaml_line_html('  - E11*') == '  - E11*'
    assert _yaml_line_html('<b>: x').startswith('<span class="k">&lt;b&gt;</span>')
