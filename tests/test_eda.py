import json
import os
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from eda_stats import run_eda_stats, separation, bin_covariate, parse_sep_detail
from eda_report import run_eda_report, select_profiles, compute_flags, attach_flags
from eda_profile import downsample
from resolve_config import EDA_DEFAULTS, resolve

N = 400
IDS = [f'P{i:04d}' for i in range(N)]
SEX = np.array([0, 1] * (N // 2))  # 0 = female, 1 = male


def _settings(**over):
    s = dict(EDA_DEFAULTS)
    s['covariates'] = ['SEX']
    s['min_cases'] = 20
    s['min_level_n'] = 20
    s.update(over)
    return s


def _eda_config(tmp_path, phenotypes=None, **over):
    cfg = {
        'settings': _settings(**over),
        'sample_id_col': 'sample_ID',
        'covariate_dictionaries': {'SEX': {'female': 0, 'male': 1}},
        'phenotypes': phenotypes or {
            'SEX': {'data_type': 'categorical', 'eda': True, 'restrict_to': None},
            'BMI': {'data_type': 'quantitative', 'eda': True, 'restrict_to': None},
            'ICD10': {'data_type': 'coded', 'eda': ['E11*'], 'restrict_to': None},
            'PCOS': {'data_type': 'advanced_coded', 'eda': False,
                     'restrict_to': {'SEX': 'female'}},
        },
    }
    p = tmp_path / 'eda_config.json'
    p.write_text(json.dumps(cfg))
    return str(p)


def _write(tmp_path, name, df):
    p = tmp_path / name
    df.to_csv(p, sep='\t', index=False)
    return str(p)


@pytest.fixture
def run_dir(tmp_path):
    rng = np.random.default_rng(0)
    files = {}
    files['SEX'] = _write(tmp_path, 'SEX.cat.tsv', pd.DataFrame({'sample_ID': IDS, 'SEX': SEX}))

    # N80.0: endometriosis — cases only among females
    n80 = np.where((SEX == 0) & (rng.random(N) < 0.3), 1, 0)
    files['N80.0'] = _write(tmp_path, 'N80.0.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'N80.0': n80}))
    # E11.9: common, sex-balanced
    e11 = (rng.random(N) < 0.3).astype(int)
    files['E11.9'] = _write(tmp_path, 'E11.9.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'E11.9': e11}))
    # E11.8: duplicate of E11.9 (pairwise)
    files['E11.8'] = _write(tmp_path, 'E11.8.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'E11.8': e11}))
    # I10: rare
    i10 = np.zeros(N, dtype=int)
    i10[:5] = 1
    files['I10'] = _write(tmp_path, 'I10.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'I10': i10}))
    # BMI: measured only in males (half the samples)
    male_ids = [s for s, x in zip(IDS, SEX) if x == 1][:180]
    files['BMI'] = _write(tmp_path, 'BMI.quant.tsv', pd.DataFrame(
        {'sample_ID': male_ids, 'BMI_mean': rng.normal(28, 4, len(male_ids))}))
    # PCOS advanced coded with restrict_to (and, here, also separated by SEX)
    pcos = np.where((SEX == 0) & (rng.random(N) < 0.2), 1, 0)
    files['PCOS'] = _write(tmp_path, 'PCOS.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'PCOS': pcos}))

    manifest = pd.DataFrame([
        ['SEX', 'SEX', 'categorical', ''],
        ['N80.0', 'ICD10', 'coded', 'N80.0'],
        ['E11.9', 'ICD10', 'coded', 'E11.9'],
        ['E11.8', 'ICD10', 'coded', 'E11.8'],
        ['I10', 'ICD10', 'coded', 'I10'],
        ['BMI_mean', 'BMI', 'quantitative', ''],
        ['PCOS', 'PCOS', 'advanced_coded', ''],
    ], columns=['column_name', 'phenotype', 'data_type', 'code'])
    manifest_path = _write(tmp_path, 'phenotype_manifest.tsv', manifest)
    return tmp_path, files, manifest_path


def _run(run_dir, **over):
    tmp_path, files, manifest_path = run_dir
    cfg = _eda_config(tmp_path, **over)
    stats_path = str(tmp_path / 'stats.tsv')
    cases_path = str(tmp_path / 'b.cases.tsv')
    run_eda_stats(list(files.values()), cfg, stats_path, covariate_paths=[files['SEX']],
                  cases_output=cases_path)
    out = tmp_path / 'eda'
    run_eda_report(stats_path, manifest_path, cfg, str(out), [cases_path])
    summary = pd.read_csv(out / 'eda_summary.tsv', sep='\t', dtype={'sep_level': str})
    return out, summary.set_index('column_name')


# ── eda_stats ─────────────────────────────────────────────────────────────────

def test_separation_complete_and_quasi():
    s = _settings()
    # all 60 cases in level 0 of two equal 200-sample levels
    assert separation(np.array([200., 200.]), np.array([60., 0.]), s)[:3] == ('complete', [0], 'cases')
    # 96% of cases in level 0 (default separation_threshold 0.95)
    assert separation(np.array([200., 200.]), np.array([96., 4.]), s)[0] == 'quasi'
    # ...but not at a stricter threshold
    assert separation(np.array([200., 200.]), np.array([96., 4.]), _settings(separation_threshold=0.99)) is None
    # balanced → none
    assert separation(np.array([200., 200.]), np.array([30., 30.]), s) is None
    # too rare: expected cases outside level 0 below min_expected
    assert separation(np.array([200., 200.]), np.array([3., 0.]), s) is None


def test_separation_subset_of_levels():
    s = _settings()
    n = np.array([100., 100., 100., 100., 100.])   # e.g. five cohorts
    # all cases in cohorts 1 and 3
    sep_type, idx, outcome, share = separation(n, np.array([0., 30., 0., 20., 0.]), s)
    assert (sep_type, idx, outcome, share) == ('complete', [1, 3], 'cases', 1.0)
    # one level wins over a subset
    assert separation(n, np.array([0., 50., 0., 0., 0.]), s)[1] == [1]
    # limited by max_subset_levels
    assert separation(n, np.array([0., 30., 0., 20., 0.]), _settings(max_subset_levels=1)) is None
    # cases everywhere → none
    assert separation(n, np.array([10., 10., 10., 10., 10.]), s) is None


def test_separation_min_cases():
    s = _settings()
    n = np.array([200., 200.])
    # 2 cases, both female: rare → never flagged (separation_min_cases 10, min_expected 5)
    assert separation(n, np.array([2., 0.]), s) is None
    assert separation(n, np.array([2., 0.]), _settings(separation_min_cases=1, min_expected=0)) is not None
    # 12 cases all female: flagged by default, not with a higher minimum
    assert separation(n, np.array([12., 0.]), s)[0] == 'complete'
    assert separation(n, np.array([12., 0.]), _settings(separation_min_cases=20)) is None


def test_separation_interactions(tmp_path):
    cohort = np.where(np.arange(N) < N // 2, 1, 2)       # independent of SEX
    rng = np.random.default_rng(1)
    in_a_female = (cohort == 1) & (SEX == 0)
    x = np.where(in_a_female & (rng.random(N) < 0.4), 1, 0)  # only cohort 1 females
    y = np.where((cohort == 1) & (rng.random(N) < 0.4), 1, 0)  # cohort 1, both sexes
    sex_f = _write(tmp_path, 'SEX.cat.tsv', pd.DataFrame({'sample_ID': IDS, 'SEX': SEX}))
    coh_f = _write(tmp_path, 'COHORT.cat.tsv', pd.DataFrame({'sample_ID': IDS, 'COHORT': cohort}))
    x_f = _write(tmp_path, 'X.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'X': x}))
    y_f = _write(tmp_path, 'Y.diag.tsv', pd.DataFrame({'sample_ID': IDS, 'Y': y}))
    cfg = _eda_config(tmp_path, covariates=['SEX', 'COHORT'],
                      separation_interactions=[['COHORT', 'SEX']])
    stats_path, cov_out = str(tmp_path / 's.tsv'), str(tmp_path / 'cov.tsv')
    run_eda_stats([x_f, y_f], cfg, stats_path, covariate_paths=[sex_f, coh_f],
                  covariates_output=cov_out)
    s = pd.read_csv(stats_path, sep='\t', dtype=str).set_index('column_name')
    cov = pd.read_csv(cov_out, sep='\t', dtype=str)
    assert 'COHORT&SEX' in cov.columns and cov.loc[0, 'COHORT&SEX'] == '1&0'
    # X: the interaction pins it to cohort 1 females
    assert 'COHORT&SEX=1&0:complete:cases(1.000)' in s.loc['X', 'sep_detail']
    # Y: COHORT=1 alone explains it; COHORT&SEX=1&0|1&1 is dropped as redundant
    assert 'COHORT=1:complete' in s.loc['Y', 'sep_detail']
    assert 'COHORT&SEX' not in s.loc['Y', 'sep_detail']

    # the filter applies it via eda_covariates.tsv like any other covariate
    wide = _write(tmp_path, 'wide.tsv', pd.DataFrame({'sample_ID': IDS, 'X': x}))
    summary = _write(tmp_path, 'summary.tsv', s.reset_index()[['column_name', 'sep_detail']]
                     .assign(flags=''))
    (tmp_path / 'f.yaml').write_text('restrict_separated: {min_share: 0.95, covariates: [COHORT&SEX]}\n')
    table, log = run_filter(str(tmp_path / 'f.yaml'), wide, str(tmp_path / 'filtered'),
                            cov_out, summary)
    kept = table['X'] != ''
    assert set(cov.loc[kept.values, 'COHORT&SEX']) == {'1&0'}
    assert 'COHORT&SEX == 1&0' in log.set_index('column_name').loc['X', 'reason']


def test_parse_sep_detail():
    d = parse_sep_detail('SEX=0:complete:cases(1.000);COHORT=A|B:quasi:cases(0.970)')
    assert d[1] == {'covariate': 'COHORT', 'levels': ['A', 'B'], 'sep_type': 'quasi',
                    'outcome': 'cases', 'share': 0.97}
    assert parse_sep_detail(None) == [] and parse_sep_detail('') == []


def test_bin_covariate_quantiles():
    s = _settings(age_bins=4)
    binned = bin_covariate(pd.Series(np.arange(100, dtype=float)), s)
    assert binned.nunique() == 4
    few = bin_covariate(pd.Series([0.0, 1.0, 1.0, np.nan]), s)
    assert list(few[:3]) == ['0', '1', '1'] and pd.isna(few[3])


def test_stats_metrics(run_dir):
    _, summary = _run(run_dir)
    assert summary.loc['N80.0', 'sep_type'] == 'complete'
    assert summary.loc['N80.0', 'sep_level'] == '0'
    assert pd.isna(summary.loc['E11.9', 'sep_type'])
    assert summary.loc['I10', 'n_cases'] == 5
    # BMI measured in 180 of 400 universe samples, only males
    assert summary.loc['BMI_mean', 'missing_frac'] == pytest.approx(220 / 400)
    assert summary.loc['BMI_mean', 'miss_gap'] > 0.8
    assert bool(summary.loc['SEX', 'is_covariate'])


def test_no_covariates_skips_separation(run_dir):
    tmp_path, files, _ = run_dir
    cfg = _eda_config(tmp_path, covariates=[])
    out = str(tmp_path / 's.tsv')
    run_eda_stats([files['N80.0']], cfg, out)
    s = pd.read_csv(out, sep='\t')
    assert s['sep_type'].isna().all()
    assert s.loc[0, 'n_cases'] > 0


# ── eda_report ────────────────────────────────────────────────────────────────

def test_flags(run_dir):
    out, summary = _run(run_dir, pairwise=True)
    flags = summary['flags'].fillna('')
    assert 'separation_complete' in flags['N80.0']
    assert 'low_cases' in flags['I10']
    assert 'high_missing' in flags['BMI_mean']
    assert 'missing_by_covariate' in flags['BMI_mean']
    assert 'pairwise_duplicate' in flags['E11.9']
    pairs = pd.read_csv(out / 'pairwise_overlap.tsv', sep='\t')
    assert set(pairs.iloc[0][['column_a', 'column_b']]) == {'E11.9', 'E11.8'}
    assert (out / 'eda_report.html').read_text().startswith('<!doctype html>')


def test_profile_selection_optin_topk_cap(run_dir):
    out, _ = _run(run_dir, profile_flagged_top_k=1)
    sel = pd.read_csv(out / 'profile_selection.tsv', sep='\t')
    reasons = dict(zip(sel['column_name'], sel['reason']))
    # cascade: BMI/SEX opted in, ICD10 only E11* globs; top-1 flagged = N80.0
    assert reasons['BMI_mean'] == 'opt_in'
    assert reasons['E11.9'] == 'opt_in' and reasons['E11.8'] == 'opt_in'
    assert 'I10' not in reasons
    assert reasons['N80.0'] == 'flagged'
    assert set(sel.loc[sel['phenotype'] == 'ICD10', 'group']) == {'ICD10'}

    out, _ = _run(run_dir, max_profile_columns=2)
    sel = pd.read_csv(out / 'profile_selection.tsv', sep='\t')
    assert len(sel) == 2


def test_profile_group_size_splits():
    stats = pd.DataFrame({
        'column_name': [f'C{i}' for i in range(5)], 'source_file': 'x.diag.tsv',
        'phenotype': 'ICD10', 'code': [f'C{i}' for i in range(5)], 'severity': 0,
    })
    cfg = {'settings': _settings(profile_group_size=2),
           'phenotypes': {'ICD10': {'eda': True}}}
    sel = select_profiles(stats, cfg)
    assert sorted(sel['group'].unique()) == ['ICD10_001', 'ICD10_002', 'ICD10_003']


def test_restriction_suggestions(run_dir):
    out, _ = _run(run_dir)
    sugg = pd.read_csv(out / 'restriction_suggestions.tsv', sep='\t', dtype=str)
    n80 = sugg[sugg['column_name'] == 'N80.0'].iloc[0]
    assert n80['sources'] == 'data'             # data-driven only: no code/name patterns
    assert n80['level'] == '0' and n80['level_label'] == 'female'
    assert float(n80['share']) == 1.0
    pcos = sugg[sugg['column_name'] == 'PCOS'].iloc[0]
    assert pcos['sources'] == 'data,yaml'
    assert 'E11.9' not in set(sugg['column_name'])
    assert 'SEX' not in set(sugg['column_name'])


def test_no_controls_flag():
    # 10 cases, 990 NA non-cases: no_controls only, not high_missing / missing_by_covariate
    stats = pd.DataFrame([{
        'column_name': 'X', 'source_file': 'X.diag.tsv', 'kind': 'binary', 'n': 10,
        'n_cases': 10, 'n_controls': 0, 'missing_frac': 0.99, 'n_unique': 1,
        'sep_type': None, 'sep_detail': None, 'miss_gap': 0.9,
        'miss_gap_covariate': 'SEX', 'miss_gap_detail': '0=0.1;1=1.0',
    }])
    flags = compute_flags(stats, _settings())
    assert list(flags['flag']) == ['no_controls']
    assert attach_flags(stats, flags).loc[0, 'severity'] == 4


# ── eda_profile ───────────────────────────────────────────────────────────────

def test_downsample_keeps_cases():
    df = pd.DataFrame({'rare': [1] * 10 + [0] * 990}, index=[f'S{i}' for i in range(1000)])
    small = downsample(df, 100, seed=1)
    assert len(small) == 100
    assert small['rare'].sum() == 10
    assert len(downsample(df, 5000)) == 1000


# ── resolve_config ────────────────────────────────────────────────────────────

def _resolve(tmp_path, cfg):
    (tmp_path / 't.tsv').write_text('id\tv\nA\t1\n')
    p = tmp_path / 'c.yaml'
    p.write_text(yaml.safe_dump(cfg))
    out = tmp_path / 'out'
    resolve(str(p), str(out))
    return out


def _base_cfg():
    return {
        'global': {'sample_id_col': 'id'},
        'phenotypes': {
            'sex': {'table': 't.tsv', 'data_type': 'categorical', 'value_col': 'v',
                    'dictionary': {'Female': 2, 'Male': 1}},
            'ldl': {'table': 't.tsv', 'data_type': 'quantitative', 'value_col': 'v'},
            'icd': {'table': 't.tsv', 'data_type': 'coded', 'value_col': 'v'},
            'pcos': {'table': 't.tsv', 'data_type': 'coded', 'value_col': 'v',
                     'eda': ['E28*'], 'restrict_to': {'sex': 'Female'}},
        },
    }


def test_resolve_eda_cascade_and_config(tmp_path):
    cfg = _base_cfg()
    cfg['eda'] = {'covariates': ['sex'], 'min_cases': 50}
    out = _resolve(tmp_path, cfg)
    eda = json.loads((out / 'eda_config.json').read_text())
    assert eda['settings']['min_cases'] == 50
    assert eda['settings']['separation_threshold'] == 0.95
    assert 'restriction_rules' not in eda['settings']
    assert eda['covariate_dictionaries'] == {'sex': {'Female': 2, 'Male': 1}}
    ph = eda['phenotypes']
    assert ph['ldl']['eda'] is True and ph['sex']['eda'] is True
    assert ph['icd']['eda'] is False
    assert ph['pcos']['eda'] == ['E28*']
    assert ph['pcos']['restrict_to'] == {'sex': 'Female'}
    assert json.loads((out / 'coded' / 'icd.json').read_text())['eda'] is False


def test_resolve_separation_interactions_auto(tmp_path):
    cfg = _base_cfg()
    cfg['eda'] = {'covariates': ['sex', 'ldl'], 'separation_interactions': 'auto'}
    out = _resolve(tmp_path, cfg)
    s = json.loads((out / 'eda_config.json').read_text())['settings']
    assert s['separation_interactions'] == [['sex', 'ldl']]
    assert s['separation_min_cases'] == 10


def test_resolve_eda_global_override(tmp_path):
    cfg = _base_cfg()
    cfg['global']['eda'] = False
    cfg['data_types'] = {'coded': {'eda': True}}
    cfg['eda'] = {'covariates': ['sex']}
    out = _resolve(tmp_path, cfg)
    ph = json.loads((out / 'eda_config.json').read_text())['phenotypes']
    assert ph['ldl']['eda'] is False
    assert ph['icd']['eda'] is True
    assert ph['pcos']['eda'] == ['E28*']


@pytest.mark.parametrize('eda_block', [
    {'covariates': ['nope']},
    {'covariates': ['icd']},
    {'bogus_key': 1},
    {'restriction_mode': 'apply'},   # removed setting → unknown key
    {'quasi_threshold': 0.99},       # renamed → separation_threshold
    {'restriction_rules': []},       # removed (data-driven only)
    {'max_missing': 2},
    {'separation_threshold': 1.5},
    {'max_subset_levels': 0},
    {'separation_min_cases': 0},
    {'separation_interactions': [['sex', 'nope']]},
    {'separation_interactions': [['sex']]},
    {'separation_interactions': 'all'},
])
def test_resolve_eda_block_errors(tmp_path, eda_block):
    cfg = _base_cfg()
    cfg['eda'] = eda_block
    if 'covariates' not in eda_block:
        cfg['eda']['covariates'] = ['sex']
    with pytest.raises(SystemExit):
        _resolve(tmp_path, cfg)


def test_resolve_restrict_to_unknown_covariate(tmp_path):
    cfg = _base_cfg()
    cfg['eda'] = {'covariates': ['ldl']}
    with pytest.raises(SystemExit):
        _resolve(tmp_path, cfg)


def test_resolve_bad_eda_value(tmp_path):
    cfg = _base_cfg()
    cfg['phenotypes']['ldl']['eda'] = 'yes'
    with pytest.raises(SystemExit):
        _resolve(tmp_path, cfg)


# ── flagged manifest / suggested filter ───────────────────────────────────────

def test_flagged_manifest_and_suggested_yaml(run_dir):
    tmp_path, files, manifest_path = run_dir
    cfg = _eda_config(tmp_path, pairwise=True)
    stats_path = str(tmp_path / 'stats.tsv')
    cases_path = str(tmp_path / 'b.cases.tsv')
    cov_out = str(tmp_path / 'eda_covariates.tsv')
    run_eda_stats(list(files.values()), cfg, stats_path, covariate_paths=[files['SEX']],
                  cases_output=cases_path, covariates_output=cov_out)
    out = tmp_path / 'eda'
    mf_out = str(tmp_path / 'manifest_flagged.tsv')
    run_eda_report(stats_path, manifest_path, cfg, str(out), [cases_path], mf_out)

    mf = pd.read_csv(mf_out, sep='\t', dtype=str, keep_default_na=False)
    assert list(mf.columns) == ['column_name', 'phenotype', 'data_type', 'code', 'flags']
    flags = dict(zip(mf['column_name'], mf['flags']))
    assert 'low_cases' in flags['I10'] and flags['SEX'] == ''

    cov = pd.read_csv(cov_out, sep='\t', dtype=str)
    assert list(cov.columns) == ['sample_ID', 'SEX'] and len(cov) == N

    text = (out / 'filter_suggested.yaml').read_text()
    sug = yaml.safe_load(text)
    assert sug['drop_columns'] == ['I10', 'BMI_mean']
    # only restrict_to is active; data-driven separation goes through restrict_separated
    assert sug['restrictions'] == {'PCOS': {'covariate': 'SEX', 'level': '0'}}
    assert 'restrict_separated' not in sug
    assert '# restrict_separated: {min_share: 0.95, covariates: ["SEX"]}' in text
    assert '#   "N80.0"  SEX = 0 (female)  share 1.000' in text
    assert '# - "E11.8"' in text   # pairwise duplicate left for review, not active


# ── filter_phenotypes ─────────────────────────────────────────────────────────

from filter_phenotypes import run_filter


def _filter_setup(tmp_path, cfg):
    wide = pd.DataFrame({
        'sample_ID': ['A', 'B', 'C', 'D'],
        'N80.0': ['1', '0', '1', ''],      # C is male with a case → miscode
        'E11.9': ['1', '0', '0', '1'],
        'I10':   ['0', '0', '1', '0'],
        'BMI_mean': ['27.50', '', '31.0', '22.1'],
    })
    wide_path = _write(tmp_path, 'wide.tsv', wide)
    cov_path = _write(tmp_path, 'cov.tsv', pd.DataFrame(
        {'sample_ID': ['A', 'B', 'C', 'D'], 'SEX': ['0', '0', '1', ''],   # D: unknown sex
         'SITE': ['X', 'Y', 'Z', 'X']}))
    summary_path = _write(tmp_path, 'summary.tsv', pd.DataFrame({
        'column_name': ['N80.0', 'E11.9', 'I10', 'BMI_mean'],
        'flags': ['separation_complete', 'separation_quasi', 'low_cases',
                  'high_missing,missing_by_covariate'],
        'sep_detail': ['SEX=0:complete:cases(1.000)', 'SITE=X|Y:quasi:cases(0.960)', '', '']}))
    manifest_path = _write(tmp_path, 'mf.tsv', pd.DataFrame({
        'column_name': ['N80.0', 'E11.9', 'I10', 'BMI_mean'], 'phenotype': 'x',
        'data_type': 'coded', 'code': '', 'flags': ''}))
    cfg_path = tmp_path / 'filter.yaml'
    cfg_path.write_text(cfg)
    out = tmp_path / 'filtered'
    table, log = run_filter(str(cfg_path), wide_path, str(out), cov_path, summary_path,
                            manifest_path)
    return wide_path, out, table.set_index('sample_ID'), log.set_index('column_name')


def test_filter_drop_keep_restrict(tmp_path):
    cfg = """
drop_flags: [high_missing]
drop_columns: ["I1*"]
keep_columns: []
restrictions:
  "N80.0": {covariate: SEX, level: "0"}
"""
    wide_path, out, table, log = _filter_setup(tmp_path, cfg)
    assert list(table.columns) == ['N80.0', 'E11.9']
    # male C (a case) and unknown-sex D set to NA; females untouched
    assert list(table['N80.0']) == ['1', '0', '', '']
    assert log.loc['N80.0', 'n_set_na'] == 1 and log.loc['N80.0', 'n_cases_set_na'] == 1
    assert log.loc['I10', 'reason'] == 'drop_columns ~ I1*'
    assert log.loc['BMI_mean', 'reason'] == 'drop_flags: high_missing'
    # raw input untouched; kept values byte-identical
    raw = pd.read_csv(wide_path, sep='\t', dtype=str, keep_default_na=False)
    assert list(raw['N80.0']) == ['1', '0', '1', '']
    mf = pd.read_csv(out / 'phenotype_manifest.filtered.tsv', sep='\t')
    assert list(mf['column_name']) == ['N80.0', 'E11.9']


def test_filter_keep_overrides_and_unknown_keep(tmp_path):
    cfg = """
restrict_unknown: keep
drop_columns: [BMI_mean]
keep_columns: [BMI_mean]
restrictions:
  "N80.0": {covariate: SEX, level: "0"}
"""
    _, _, table, log = _filter_setup(tmp_path, cfg)
    assert 'BMI_mean' in table.columns and table.loc['A', 'BMI_mean'] == '27.50'
    assert log.loc['BMI_mean', 'action'] == 'kept'
    assert list(table['N80.0']) == ['1', '0', '', '']   # D was already NA


def test_filter_commented_out_sections(tmp_path):
    # suggested files can leave a section with only commented entries → None
    cfg = "drop_columns:\n  # - \"I10\"\nrestrictions:\n"
    _, _, table, log = _filter_setup(tmp_path, cfg)
    assert list(table.columns) == ['N80.0', 'E11.9', 'I10', 'BMI_mean']
    assert log.empty


@pytest.mark.parametrize('cfg', [
    'bogus: 1\n',
    'restrict_unknown: maybe\n',
    'restrictions:\n  X: {covariate: SEX}\n',
    'restrictions:\n  "N80.0": {covariate: AGE, level: "1"}\n',
    'restrictions:\n  "N80.0": {covariate: SEX, level: []}\n',
    'restrict_separated: {covariates: [SEX]}\n',
    'restrict_separated: {min_share: 2}\n',
    'restrict_separated: {min_share: 0.9, bogus: 1}\n',
])
def test_filter_config_errors(tmp_path, cfg):
    with pytest.raises(ValueError):
        _filter_setup(tmp_path, cfg)


def test_filter_level_list(tmp_path):
    cfg = 'restrictions:\n  "E11.9": {covariate: SITE, level: [X, Y]}\n'
    _, _, table, log = _filter_setup(tmp_path, cfg)
    assert list(table['E11.9']) == ['1', '0', '', '1']    # C is SITE Z
    assert log.loc['E11.9', 'reason'].startswith('SITE == X|Y')


def test_filter_restrict_separated(tmp_path):
    # both separations ≥ 0.95, all covariates
    _, _, table, log = _filter_setup(tmp_path, 'restrict_separated: {min_share: 0.95}\n')
    assert list(table['N80.0']) == ['1', '0', '', '']     # SEX == 0; D unknown → NA
    assert list(table['E11.9']) == ['1', '0', '', '1']    # SITE in X|Y
    assert 'separated, share 1.000' in log.loc['N80.0', 'reason']
    assert log.loc['N80.0', 'n_cases_set_na'] == 1

    # threshold above E11.9's share, and covariate filter
    _, _, table, log = _filter_setup(tmp_path, 'restrict_separated: {min_share: 0.99}\n')
    assert 'E11.9' not in log.index
    _, _, table, log = _filter_setup(
        tmp_path, 'restrict_separated: {min_share: 0.9, covariates: [SITE]}\n')
    assert set(log.index) == {'E11.9'}


def test_filter_explicit_restriction_wins_over_separated(tmp_path):
    cfg = """
restrictions:
  "N80.0": {covariate: SEX, level: ["0", "1"]}
restrict_separated: {min_share: 0.95}
"""
    _, _, table, log = _filter_setup(tmp_path, cfg)
    assert list(table['N80.0']) == ['1', '0', '1', '']    # explicit rule keeps male C
    assert 'separated' not in log.loc['N80.0', 'reason']
