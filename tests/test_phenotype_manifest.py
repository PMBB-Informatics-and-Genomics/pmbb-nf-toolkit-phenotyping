import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from phenotype_manifest import run_phenotype_manifest, MANIFEST_COLUMNS


def _write(tmp_path, name, columns, cfg):
    result = tmp_path / f'{name}.tsv'
    pd.DataFrame(columns=['sample_ID'] + columns).to_csv(result, sep='\t', index=False)
    cfg_path = tmp_path / f'{name}.json'
    cfg_path.write_text(json.dumps(cfg))
    out = tmp_path / f'{name}.manifest.tsv'
    run_phenotype_manifest(str(result), str(cfg_path), str(out))
    return pd.read_csv(out, sep='\t', dtype=str, keep_default_na=False)


def test_header_columns():
    assert MANIFEST_COLUMNS == ['column_name', 'phenotype', 'data_type', 'code']


def test_quantitative_one_row_per_stat_column(tmp_path):
    m = _write(tmp_path, 'BMI', ['BMI_mean', 'BMI_std'],
               {'phenotype_name': 'BMI', 'data_type': 'quantitative'})
    assert m.values.tolist() == [
        ['BMI_mean', 'BMI', 'quantitative', ''],
        ['BMI_std', 'BMI', 'quantitative', ''],
    ]


def test_categorical_one_hot_columns(tmp_path):
    m = _write(tmp_path, 'ANC', ['ANCESTRY_AFR', 'ANCESTRY_EUR'],
               {'phenotype_name': 'ANCESTRY', 'data_type': 'categorical'})
    assert m['column_name'].tolist() == ['ANCESTRY_AFR', 'ANCESTRY_EUR']
    assert set(m['data_type']) == {'categorical'}


def test_coded_code_is_raw_without_prefix(tmp_path):
    m = _write(tmp_path, 'phe', ['phe_401.1'],
               {'phenotype_name': 'PheCode', 'data_type': 'coded', 'column_prefix': 'phe_'})
    assert m.values.tolist() == [['phe_401.1', 'PheCode', 'coded', '401.1']]


def test_coded_without_prefix_code_equals_column(tmp_path):
    m = _write(tmp_path, 'icd', ['N80.0'],
               {'phenotype_name': 'ICD10', 'data_type': 'coded', 'column_prefix': ''})
    assert m.values.tolist() == [['N80.0', 'ICD10', 'coded', 'N80.0']]


def test_advanced_coded_has_no_code(tmp_path):
    m = _write(tmp_path, 'T2D', ['T2Diab'],
               {'phenotype_name': 'T2Diab', 'data_type': 'advanced_coded'})
    assert m.values.tolist() == [['T2Diab', 'T2Diab', 'advanced_coded', '']]


def test_chunked_coded_uses_yaml_key(tmp_path):
    # chunk JSONs keep the parent phenotype_name
    m = _write(tmp_path, 'ICD10_chunk_001', ['E11.9'],
               {'phenotype_name': 'ICD10', 'data_type': 'coded'})
    assert m['phenotype'].tolist() == ['ICD10']
