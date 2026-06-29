import json
import os
import sys
import pytest
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from resolve_config import resolve, merge_params, validate


# ── merge_params ──────────────────────────────────────────────────────────────

def test_merge_params_global_to_phenotype():
    """Global params cascade down when not overridden."""
    global_p = {'sample_id_col': 'PMBB_ID', 'outlier_method': 'iqr'}
    data_type_p = {'stats': ['mean']}
    pheno_p = {'table': '/tmp/bmi.tsv', 'value_col': 'BMI_VALUE'}
    result = merge_params(global_p, data_type_p, pheno_p, 'quantitative')
    assert result['sample_id_col'] == 'PMBB_ID'
    assert result['outlier_method'] == 'iqr'
    assert result['stats'] == ['mean']
    assert result['table'] == '/tmp/bmi.tsv'


def test_merge_params_phenotype_overrides_data_type():
    """Phenotype-level stats override data_type-level stats."""
    global_p = {}
    data_type_p = {'stats': ['mean', 'std']}
    pheno_p = {'stats': ['median']}
    result = merge_params(global_p, data_type_p, pheno_p, 'quantitative')
    assert result['stats'] == ['median']


def test_merge_params_data_type_overrides_global():
    """DataType-level outlier_method overrides global."""
    global_p = {'outlier_method': 'iqr'}
    data_type_p = {'outlier_method': 'none'}
    pheno_p = {}
    result = merge_params(global_p, data_type_p, pheno_p, 'quantitative')
    assert result['outlier_method'] == 'none'


# ── validate ──────────────────────────────────────────────────────────────────

def test_validate_missing_table_no_preprocessed_path():
    errors = validate('BMI', {'data_type': 'quantitative', 'value_col': 'X'})
    assert any('table' in e for e in errors)


def test_validate_preprocessed_path_skips_table_check():
    errors = validate('BMI', {
        'data_type': 'quantitative',
        'value_col': 'X',
        'preprocessed_path': '/some/file.tsv',
    })
    assert errors == []


def test_validate_bad_data_type():
    errors = validate('BMI', {
        'data_type': 'UNKNOWN',
        'table': '/x.tsv',
        'value_col': 'x',
    })
    assert any('data_type' in e for e in errors)


# ── resolve (integration) ─────────────────────────────────────────────────────

def test_resolve_writes_json_per_phenotype(tmp_path, pmbb_test_config):
    resolve(pmbb_test_config, str(tmp_path))
    # quantitative
    bmi_json = tmp_path / 'quantitative' / 'BMI.json'
    assert bmi_json.exists()
    data = json.loads(bmi_json.read_text())
    assert data['data_type'] == 'quantitative'
    assert data['phenotype_name'] == 'BMI'
    assert data['sample_id_col'] == 'PMBB_ID'          # inherited from global
    assert data['outlier_method'] == 'iqr'              # inherited from global
    assert data['stats'] == ['mean', 'std']             # from data_types.quantitative

def test_resolve_categorical_json(tmp_path, pmbb_test_config):
    resolve(pmbb_test_config, str(tmp_path))
    sex_json = tmp_path / 'categorical' / 'SEX.json'
    assert sex_json.exists()
    data = json.loads(sex_json.read_text())
    assert data['data_type'] == 'categorical'
    assert data['dictionary'] == {'female': 0, 'male': 1, 'unknown': 'NA'}

def test_resolve_coded_json(tmp_path, pmbb_test_config):
    resolve(pmbb_test_config, str(tmp_path))
    icd_json = tmp_path / 'coded' / 'ICD10.json'
    assert icd_json.exists()
    data = json.loads(icd_json.read_text())
    assert data['data_type'] == 'coded'
    assert data['subsample'] == 'all'
    assert data['filter'] == {'CODE_TYPE': ['ICD10']}

def test_resolve_output_name_defaults_to_phenotype_name(tmp_path, pmbb_test_config):
    resolve(pmbb_test_config, str(tmp_path))
    bmi_data = json.loads((tmp_path / 'quantitative' / 'BMI.json').read_text())
    assert bmi_data['output_name'] == 'BMI'

def test_resolve_output_name_override(tmp_path, pmbb_test_config):
    resolve(pmbb_test_config, str(tmp_path))
    hdlc_data = json.loads((tmp_path / 'quantitative' / 'HDLC.json').read_text())
    assert hdlc_data['output_name'] == 'HDLC'

def test_resolve_missing_data_type_exits(tmp_path, tmp_path_factory):
    bad_config = tmp_path_factory.mktemp('cfg') / 'bad.yaml'
    bad_config.write_text('phenotypes:\n  BMI:\n    table: x.tsv\n    value_col: y\n')
    with pytest.raises(SystemExit):
        resolve(str(bad_config), str(tmp_path))

def test_resolve_table_path_made_absolute(tmp_path, pmbb_test_config):
    """Relative table paths should be resolved relative to the YAML file's directory."""
    resolve(pmbb_test_config, str(tmp_path))
    bmi_data = json.loads((tmp_path / 'quantitative' / 'BMI.json').read_text())
    assert os.path.isabs(bmi_data['table'])


# ── coded chunking ───────────────────────────────────────────────────────

def test_resolve_coded_chunked_emits_multiple_jsons(tmp_path, pmbb_test_config):
    """chunk_size=2 with 3 unique ICD10 codes should produce 2 chunk JSON files."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=2)
    chunk_files = list((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    assert len(chunk_files) == 2


def test_resolve_coded_chunk_subsample_is_list(tmp_path, pmbb_test_config):
    """Each chunk JSON should have subsample as a list of specific codes."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=2)
    chunk_files = sorted((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    for f in chunk_files:
        data = json.loads(f.read_text())
        assert isinstance(data['subsample'], list)
        assert len(data['subsample']) > 0


def test_resolve_coded_chunk_covers_all_codes(tmp_path, pmbb_test_config):
    """Union of all chunk subsample lists equals the full set of unique filtered codes."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=2)
    chunk_files = sorted((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    all_codes = []
    for f in chunk_files:
        all_codes.extend(json.loads(f.read_text())['subsample'])
    # icd.tsv has 3 ICD10 codes after CODE_TYPE filter: N80.0, E11.9, I11.9
    assert set(all_codes) == {'N80.0', 'E11.9', 'I11.9'}


def test_resolve_coded_chunk_no_overlap(tmp_path, pmbb_test_config):
    """No code should appear in more than one chunk."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=2)
    chunk_files = sorted((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    seen = []
    for f in chunk_files:
        codes = json.loads(f.read_text())['subsample']
        for code in codes:
            assert code not in seen, f"Duplicate code {code} across chunks"
            seen.append(code)


def test_resolve_coded_chunk_zero_disables_chunking(tmp_path, pmbb_test_config):
    """chunk_size=0 should produce the original single ICD10.json with subsample='all'."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=0)
    assert (tmp_path / 'coded' / 'ICD10.json').exists()
    chunk_files = list((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    assert len(chunk_files) == 0
    data = json.loads((tmp_path / 'coded' / 'ICD10.json').read_text())
    assert data['subsample'] == 'all'


def test_resolve_coded_chunk_respects_filter(tmp_path, pmbb_test_config):
    """Chunking should honour the phenotype's filter — ICD9 code '617' must be excluded."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=10)
    chunk_files = list((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    all_codes = []
    for f in chunk_files:
        all_codes.extend(json.loads(f.read_text())['subsample'])
    assert '617' not in all_codes


def test_resolve_coded_chunk_preserves_other_fields(tmp_path, pmbb_test_config):
    """Chunk JSONs must preserve data_type, table, value_col, and other config fields."""
    resolve(pmbb_test_config, str(tmp_path), coded_chunk_size=2)
    chunk_files = sorted((tmp_path / 'coded').glob('ICD10_chunk_*.json'))
    for f in chunk_files:
        data = json.loads(f.read_text())
        assert data['data_type'] == 'coded'
        assert data['phenotype_name'] == 'ICD10'
        assert 'table' in data
        assert 'value_col' in data


# ── sep inheritance ───────────────────────────────────────────────────────────

def test_sep_inherits_from_global(tmp_path):
    """sep set at global level flows into resolved phenotype JSON."""
    config = {
        'global': {'sample_id_col': 'PID', 'sep': 'csv'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'BMI': {'data_type': 'quantitative', 'table': '/fake/bmi.csv', 'value_col': 'VAL'}
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    out_dir = tmp_path / 'out'
    resolve(str(cfg_file), str(out_dir))
    resolved = json.loads((out_dir / 'quantitative' / 'BMI.json').read_text())
    assert resolved['sep'] == 'csv'


def test_sep_phenotype_overrides_global(tmp_path):
    """Phenotype-level sep overrides global sep."""
    config = {
        'global': {'sample_id_col': 'PID', 'sep': 'csv'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'BMI': {'data_type': 'quantitative', 'table': '/fake/bmi.tsv',
                    'value_col': 'VAL', 'sep': 'tsv'}
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    out_dir = tmp_path / 'out'
    resolve(str(cfg_file), str(out_dir))
    resolved = json.loads((out_dir / 'quantitative' / 'BMI.json').read_text())
    assert resolved['sep'] == 'tsv'


def test_sep_absent_resolves_to_none(tmp_path):
    """When sep is not set anywhere, resolved JSON has sep: null."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'BMI': {'data_type': 'quantitative', 'table': '/fake/bmi.tsv', 'value_col': 'VAL'}
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    out_dir = tmp_path / 'out'
    resolve(str(cfg_file), str(out_dir))
    resolved = json.loads((out_dir / 'quantitative' / 'BMI.json').read_text())
    assert resolved['sep'] is None


# ── event_mask ───────────────────────────────────────────────────────────────

def test_event_mask_defaults_to_empty_list(tmp_path):
    """Phenotype without event_mask in YAML gets event_mask: [] in output JSON."""
    config = {
        'global': {'sample_id_col': 'PMBB_ID'},
        'phenotypes': {
            'HbA1c': {
                'data_type': 'quantitative',
                'table': '/fake/table.tsv',
                'value_col': 'VALUE',
            }
        }
    }
    cfg_path = tmp_path / 'config.yaml'
    cfg_path.write_text(yaml.dump(config))
    resolve(str(cfg_path), str(tmp_path))
    out_json = tmp_path / 'quantitative' / 'HbA1c.json'
    data = json.loads(out_json.read_text())
    assert 'event_mask' in data
    assert data['event_mask'] == []


def test_event_mask_preserved_from_yaml(tmp_path):
    """Phenotype with event_mask in YAML has that list in output JSON."""
    config = {
        'global': {'sample_id_col': 'PMBB_ID'},
        'phenotypes': {
            'HbA1c': {
                'data_type': 'quantitative',
                'table': '/fake/table.tsv',
                'value_col': 'VALUE',
                'event_mask': ['third_trimester_pregnancy', 'first_trimester_loss'],
            }
        }
    }
    cfg_path = tmp_path / 'config.yaml'
    cfg_path.write_text(yaml.dump(config))
    resolve(str(cfg_path), str(tmp_path))
    out_json = tmp_path / 'quantitative' / 'HbA1c.json'
    data = json.loads(out_json.read_text())
    assert data['event_mask'] == ['third_trimester_pregnancy', 'first_trimester_loss']


# ── risk_classification resolution ───────────────────────────────────────────

def test_risk_data_type_default_emits_rule(tmp_path):
    """DataType-level risk_classification produces a rule in risk_classification_config.json."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {
            'quantitative': {
                'risk_classification': {'method': 'percentile', 'high': 90, 'low': 10}
            }
        },
        'phenotypes': {
            'BMI': {'data_type': 'quantitative', 'table': '/fake/bmi.tsv', 'value_col': 'VAL'}
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    risk_cfg = json.loads((tmp_path / 'risk_classification_config.json').read_text())
    assert len(risk_cfg['rules']) == 1
    rule = risk_cfg['rules'][0]
    assert rule['phenotype'] == 'BMI'
    assert rule['method'] == 'percentile'
    assert rule['high'] == 90
    assert rule['low'] == 10
    assert rule['source_col'] == 'BMI_mean'
    assert rule['output_col'] == 'BMI_risk'


def test_risk_pheno_override_replaces_data_type(tmp_path):
    """Per-phenotype risk_classification completely replaces the data_type-level rule."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {
            'quantitative': {
                'risk_classification': {'method': 'percentile', 'high': 90, 'low': 10}
            }
        },
        'phenotypes': {
            'BMI': {
                'data_type': 'quantitative', 'table': '/fake/bmi.tsv', 'value_col': 'VAL',
                'risk_classification': {'method': 'quantile_bin', 'n': 4, 'high': 1, 'low': 1},
            }
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    risk_cfg = json.loads((tmp_path / 'risk_classification_config.json').read_text())
    rule = risk_cfg['rules'][0]
    assert rule['method'] == 'quantile_bin'
    assert rule['n'] == 4
    assert rule['high'] == 1
    assert rule['low'] == 1


def test_risk_opt_out_excludes_phenotype(tmp_path):
    """risk_classification: false excludes phenotype from the rule list."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {
            'quantitative': {
                'risk_classification': {'method': 'percentile', 'high': 90}
            }
        },
        'phenotypes': {
            'BMI': {
                'data_type': 'quantitative', 'table': '/fake/bmi.tsv', 'value_col': 'VAL',
                'risk_classification': False,
            },
            'HDL': {
                'data_type': 'quantitative', 'table': '/fake/hdl.tsv', 'value_col': 'VAL',
            },
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    risk_cfg = json.loads((tmp_path / 'risk_classification_config.json').read_text())
    phenotypes = [r['phenotype'] for r in risk_cfg['rules']]
    assert 'BMI' not in phenotypes
    assert 'HDL' in phenotypes


def test_risk_no_rules_no_file(tmp_path):
    """When no phenotype has a risk rule, risk_classification_config.json is not written."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'BMI': {'data_type': 'quantitative', 'table': '/fake/bmi.tsv', 'value_col': 'VAL'}
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    assert not (tmp_path / 'risk_classification_config.json').exists()


def test_risk_custom_stat_in_source_col(tmp_path):
    """Custom stat in risk_classification is reflected in source_col."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'BP': {
                'data_type': 'quantitative', 'table': '/fake/bp.tsv', 'value_col': 'VAL',
                'risk_classification': {'method': 'threshold', 'high': 130, 'stat': 'median'},
            }
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    risk_cfg = json.loads((tmp_path / 'risk_classification_config.json').read_text())
    rule = risk_cfg['rules'][0]
    assert rule['stat'] == 'median'
    assert rule['source_col'] == 'BP_median'


def test_risk_output_name_used_in_cols(tmp_path):
    """When phenotype has output_name, source_col and output_col use output_name."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'HDLC': {
                'data_type': 'quantitative', 'table': '/fake/hdl.tsv', 'value_col': 'VAL',
                'output_name': 'HDL',
                'risk_classification': {'method': 'percentile', 'high': 90},
            }
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    risk_cfg = json.loads((tmp_path / 'risk_classification_config.json').read_text())
    rule = risk_cfg['rules'][0]
    assert rule['source_col'] == 'HDL_mean'
    assert rule['output_col'] == 'HDL_risk'


def test_risk_cutoffs_and_labels_preserved(tmp_path):
    """cutoffs list and labels are preserved verbatim in the emitted rule."""
    config = {
        'global': {'sample_id_col': 'PID'},
        'data_types': {'quantitative': {}},
        'phenotypes': {
            'BP': {
                'data_type': 'quantitative', 'table': '/fake/bp.tsv', 'value_col': 'VAL',
                'risk_classification': {
                    'method': 'percentile', 'cutoffs': [33, 67],
                    'labels': ['normal', 'elevated', 'high'],
                },
            }
        },
    }
    cfg_file = tmp_path / 'cfg.yaml'
    cfg_file.write_text(yaml.dump(config))
    resolve(str(cfg_file), str(tmp_path))
    risk_cfg = json.loads((tmp_path / 'risk_classification_config.json').read_text())
    rule = risk_cfg['rules'][0]
    assert rule['cutoffs'] == [33, 67]
    assert rule['labels'] == ['normal', 'elevated', 'high']


# ── prefix_matching defaults ──────────────────────────────────────────────────

def test_diag_defaults_include_prefix_matching_true():
    from resolve_config import CODED_DEFAULTS
    assert CODED_DEFAULTS.get('prefix_matching') is True


def test_advanced_diag_defaults_include_prefix_matching_true():
    from resolve_config import ADVANCED_CODED_DEFAULTS
    assert ADVANCED_CODED_DEFAULTS.get('prefix_matching') is True


def test_merge_params_coded_prefix_matching_default_true():
    result = merge_params({}, {}, {'table': '/tmp/x.tsv', 'value_col': 'CODE'}, 'coded')
    assert result['prefix_matching'] is True


def test_merge_params_coded_prefix_matching_phenotype_override():
    result = merge_params(
        {}, {},
        {'table': '/tmp/x.tsv', 'value_col': 'CODE', 'prefix_matching': False},
        'coded',
    )
    assert result['prefix_matching'] is False


def test_merge_params_coded_prefix_matching_data_type_override():
    result = merge_params(
        {}, {'prefix_matching': False},
        {'table': '/tmp/x.tsv', 'value_col': 'CODE'},
        'coded',
    )
    assert result['prefix_matching'] is False


def test_merge_params_advanced_coded_prefix_matching_default_true():
    result = merge_params({}, {}, {}, 'advanced_coded')
    assert result['prefix_matching'] is True


# ── prefix_matching validation ────────────────────────────────────────────────

def test_validate_coded_prefix_matching_non_bool_rejected():
    errors = validate('CPT', {
        'data_type': 'coded',
        'table': '/tmp/cpt.tsv',
        'value_col': 'CPT_CODE',
        'prefix_matching': 'yes',
    })
    assert any("'prefix_matching' must be a bool" in e for e in errors)


def test_validate_coded_prefix_matching_false_accepted():
    errors = validate('CPT', {
        'data_type': 'coded',
        'table': '/tmp/cpt.tsv',
        'value_col': 'CPT_CODE',
        'prefix_matching': False,
    })
    assert errors == []


def test_validate_advanced_coded_prefix_matching_non_bool_rejected():
    errors = validate('T2Diab', {
        'data_type': 'advanced_coded',
        'prefix_matching': 'true',
        'sources': {
            'ICD10': {'from_phenotype': 'ICD10', 'case_codes': ['E11']},
        },
    })
    assert any("'prefix_matching' must be a bool" in e for e in errors)


def test_validate_source_prefix_matching_non_bool_rejected():
    errors = validate('T2Diab', {
        'data_type': 'advanced_coded',
        'sources': {
            'CPT': {
                'from_phenotype': 'CPT_codes',
                'case_codes': ['99213'],
                'prefix_matching': 1,
            },
        },
    })
    assert any("source 'CPT'" in e and "'prefix_matching' must be a bool" in e for e in errors)


def test_validate_source_prefix_matching_false_accepted():
    errors = validate('T2Diab', {
        'data_type': 'advanced_coded',
        'sources': {
            'CPT': {
                'from_phenotype': 'CPT_codes',
                'case_codes': ['99213'],
                'prefix_matching': False,
            },
        },
    })
    assert errors == []
