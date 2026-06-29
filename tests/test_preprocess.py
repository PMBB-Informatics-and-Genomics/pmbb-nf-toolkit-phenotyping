import json
import os
import sys
from pathlib import Path
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from preprocess import (
    apply_filter, get_matching_codes, run_preprocess,
    filter_outliers_iqr, filter_outliers_zscore,
)


# ── apply_filter ──────────────────────────────────────────────────────────────

def test_apply_filter_exact_match():
    df = pd.DataFrame({'TYPE': ['HDL', 'LDL', 'TG'], 'VAL': [1, 2, 3]})
    result = apply_filter(df, {'TYPE': ['HDL']})
    assert list(result['TYPE']) == ['HDL']


def test_apply_filter_wildcard():
    df = pd.DataFrame({'TYPE': ['HDL', 'LDL', 'VLDL', 'TG'], 'VAL': [1, 2, 3, 4]})
    result = apply_filter(df, {'TYPE': ['*DL']})
    assert set(result['TYPE']) == {'HDL', 'LDL', 'VLDL'}


def test_apply_filter_multiple_columns():
    df = pd.DataFrame({
        'TYPE': ['HDL', 'HDL', 'LDL'],
        'SITE': ['A', 'B', 'A'],
        'VAL': [1, 2, 3],
    })
    result = apply_filter(df, {'TYPE': ['HDL'], 'SITE': ['A']})
    assert len(result) == 1


def test_apply_filter_empty_dict_returns_all():
    df = pd.DataFrame({'TYPE': ['HDL', 'LDL'], 'VAL': [1, 2]})
    result = apply_filter(df, {})
    assert len(result) == 2


def test_apply_filter_missing_column_skipped():
    df = pd.DataFrame({'VAL': [1, 2]})
    result = apply_filter(df, {'NONEXISTENT': ['X']})
    assert len(result) == 2


# ── get_matching_codes ────────────────────────────────────────────────────────

def test_get_matching_codes_all():
    series = pd.Series(['N80.0', 'E11.9', 'I11.9'])
    codes = get_matching_codes(series, 'all')
    assert set(codes) == {'N80.0', 'E11.9', 'I11.9'}


def test_get_matching_codes_exact_list():
    series = pd.Series(['N80.0', 'E11.9', 'I11.9', 'N80.1'])
    codes = get_matching_codes(series, ['N80.0', 'E11.9'])
    assert set(codes) == {'N80.0', 'E11.9'}


def test_get_matching_codes_wildcard():
    series = pd.Series(['N80.0', 'N80.1', 'E11.9', 'E11.0'])
    codes = get_matching_codes(series, ['N80.*', 'E11.9'])
    assert set(codes) == {'N80.0', 'N80.1', 'E11.9'}


def test_get_matching_codes_no_match_returns_empty():
    series = pd.Series(['N80.0', 'E11.9'])
    codes = get_matching_codes(series, ['Z99.*'])
    assert codes == []


# ── run_preprocess (integration) ──────────────────────────────────────────────

def _make_config(tmp_path, overrides):
    base = {
        'phenotype_name': 'TEST',
        'data_type': 'quantitative',
        'sample_id_col': 'PMBB_ID',
        'date_col': 'VISIT_DATE',
        'value_col': 'VALUE',
        'filter': {},
        'subsample': None,
        'preprocessed_path': None,
        'output_name': 'TEST',
        'min_occurrences': 1,
        'outlier_method': 'none',
        'outlier_mode': 'cap',
        'outlier_iqr_multiplier': 1.5,
        'outlier_zscore_sd': 3.0,
        'scale': False,
        'scale_method': 'zscore',
        'stats': ['mean'],
        'qcut_bins': 0,
        'reference_date': 'today',
    }
    base.update(overrides)
    p = tmp_path / 'config.json'
    p.write_text(json.dumps(base))
    return str(p)


def test_run_preprocess_quantitative(tmp_path, pmbb_vitals_bmi):
    cfg = _make_config(tmp_path, {
        'phenotype_name': 'BMI',
        'table': pmbb_vitals_bmi,
        'value_col': 'BMI_VALUE',
        'date_col': 'VISIT_DATE',
        'output_name': 'BMI',
    })
    run_preprocess(cfg, str(tmp_path))
    out = tmp_path / 'BMI.long.tsv'
    assert out.exists()
    df = pd.read_csv(out, sep='\t')
    assert set(df.columns) >= {'sample_ID', 'data_type', 'phenotype', 'value', 'occurrence_date'}
    assert all(df['phenotype'] == 'BMI')
    assert all(df['data_type'] == 'quantitative')


def test_run_preprocess_filter(tmp_path, pmbb_lipids):
    cfg = _make_config(tmp_path, {
        'phenotype_name': 'HDLC',
        'table': pmbb_lipids,
        'value_col': 'RESULT_VALUE_NUM',
        'date_col': 'RESULT_DATE',
        'filter': {'MEASUREMENT_TYPE': ['HDL']},
        'output_name': 'HDLC',
    })
    run_preprocess(cfg, str(tmp_path))
    df = pd.read_csv(tmp_path / 'HDLC.long.tsv', sep='\t')
    assert len(df) == 4   # 4 HDL rows in labs_lipids.tsv (P001, P002, P003, P004)
    assert all(df['phenotype'] == 'HDLC')


def test_run_preprocess_filter_wildcard(tmp_path, pmbb_lipids):
    cfg = _make_config(tmp_path, {
        'phenotype_name': 'DL_cholesterol',
        'table': pmbb_lipids,
        'value_col': 'RESULT_VALUE_NUM',
        'date_col': 'RESULT_DATE',
        'filter': {'MEASUREMENT_TYPE': ['*DL']},
        'output_name': 'DL_cholesterol',
    })
    run_preprocess(cfg, str(tmp_path))
    df = pd.read_csv(tmp_path / 'DL_cholesterol.long.tsv', sep='\t')
    assert len(df) == 8   # 4 HDL + 4 LDL rows (TG does not match *DL)


def test_run_preprocess_coded_all(tmp_path, pmbb_icd):
    cfg = _make_config(tmp_path, {
        'phenotype_name': 'ICD10',
        'data_type': 'coded',
        'table': pmbb_icd,
        'value_col': 'ICD_CODE',
        'date_col': 'ENCOUNTER_DATE',
        'filter': {'CODE_TYPE': ['ICD10']},
        'subsample': 'all',
        'output_name': None,
    })
    run_preprocess(cfg, str(tmp_path))
    tsv_files = list(tmp_path.glob('*.long.tsv'))
    code_names = {f.stem.replace('.long', '') for f in tsv_files}
    assert code_names == {'N80.0', 'E11.9', 'I11.9'}


def test_run_preprocess_coded_specific_codes(tmp_path, pmbb_icd):
    cfg = _make_config(tmp_path, {
        'phenotype_name': 'ICD10_subset',
        'data_type': 'coded',
        'table': pmbb_icd,
        'value_col': 'ICD_CODE',
        'date_col': 'ENCOUNTER_DATE',
        'filter': {'CODE_TYPE': ['ICD10']},
        'subsample': ['N80.*'],
        'output_name': None,
    })
    run_preprocess(cfg, str(tmp_path))
    tsv_files = list(tmp_path.glob('*.long.tsv'))
    code_names = {f.stem.replace('.long', '') for f in tsv_files}
    assert code_names == {'N80.0'}


def test_run_preprocess_preprocessed_path(tmp_path):
    already_done = tmp_path / 'already.tsv'
    already_done.write_text('sample_ID\tphenotype\tvalue\nP001\tBMI\t28.4\n')
    cfg = _make_config(tmp_path, {
        'phenotype_name': 'BMI',
        'preprocessed_path': str(already_done),
        'table': None,
        'value_col': None,
    })
    run_preprocess(cfg, str(tmp_path))
    out = tmp_path / 'BMI.long.tsv'
    assert out.exists()
    assert out.read_text() == already_done.read_text()


def test_run_preprocess_outlier_iqr(tmp_path):
    """IQR outlier capping is applied to quantitative data."""
    data = 'PMBB_ID\tVALUE\tVISIT_DATE\n'
    data += 'P001\t10.0\t2020-01-01\n'
    data += 'P002\t12.0\t2020-01-01\n'
    data += 'P003\t11.0\t2020-01-01\n'
    data += 'P004\t10.5\t2020-01-01\n'
    data += 'P005\t9999.0\t2020-01-01\n'
    table = tmp_path / 'quant.tsv'
    table.write_text(data)
    cfg = _make_config(tmp_path, {
        'table': str(table),
        'outlier_method': 'iqr',
        'outlier_mode': 'cap',
        'outlier_iqr_multiplier': 1.5,
    })
    run_preprocess(cfg, str(tmp_path))
    df = pd.read_csv(tmp_path / 'TEST.long.tsv', sep='\t')
    values = pd.to_numeric(df['value'])
    assert values.max() < 9999.0


# ── _read_table sep_config ────────────────────────────────────────────────────

def test_read_table_csv_with_sep_config(tmp_path):
    """sep_config='csv' reads a comma-separated file with a .tsv extension."""
    f = tmp_path / 'data.tsv'   # extension would normally trigger tab detection
    f.write_text('ID,VAL\nP001,1.2\nP002,3.4\n')
    from preprocess import _read_table
    df = _read_table(str(f), sep_config='csv')
    assert list(df.columns) == ['ID', 'VAL']
    assert df['VAL'].tolist() == ['1.2', '3.4']

def test_read_table_pipe_with_sep_config(tmp_path):
    """sep_config='pipe' reads a pipe-delimited file."""
    f = tmp_path / 'data.txt'
    f.write_text('ID|VAL\nP001|1.2\nP002|3.4\n')
    from preprocess import _read_table
    df = _read_table(str(f), sep_config='pipe')
    assert list(df.columns) == ['ID', 'VAL']
    assert df['VAL'].tolist() == ['1.2', '3.4']

def test_read_table_autodetect_tsv(tmp_path):
    """No sep_config: .tsv extension still auto-detects as tab."""
    f = tmp_path / 'data.tsv'
    f.write_text('ID\tVAL\nP001\t1.2\n')
    from preprocess import _read_table
    df = _read_table(str(f))
    assert list(df.columns) == ['ID', 'VAL']

def test_read_table_autodetect_csv(tmp_path):
    """No sep_config: .csv extension auto-detects as comma."""
    f = tmp_path / 'data.csv'
    f.write_text('ID,VAL\nP001,1.2\n')
    from preprocess import _read_table
    df = _read_table(str(f))
    assert list(df.columns) == ['ID', 'VAL']


# ── sample_list filtering ──────────────────────────────────────────────────────

def test_run_preprocess_sample_list_filters_quantitative(tmp_path, pmbb_test_dir):
    """Only samples in sample_list appear in output."""
    table = Path(pmbb_test_dir) / 'vitals_bmi.tsv'
    sample_list = tmp_path / 'samples.txt'
    sample_list.write_text('P001\nP002\n')
    cfg = _make_config(tmp_path, {
        'table': str(table),
        'data_type': 'quantitative',
        'sample_id_col': 'PMBB_ID',
        'date_col': 'VISIT_DATE',
        'value_col': 'BMI_VALUE',
        'sep': 'tsv',
    })
    run_preprocess(cfg, str(tmp_path), sample_list_path=str(sample_list))
    result = pd.read_csv(tmp_path / 'TEST.long.tsv', sep='\t')
    assert set(result['sample_ID']) <= {'P001', 'P002'}
    assert 'P003' not in result['sample_ID'].values
    assert 'P004' not in result['sample_ID'].values


def test_run_preprocess_no_sample_list_keeps_all(tmp_path, pmbb_test_dir):
    """Without sample_list, all samples appear in output."""
    table = Path(pmbb_test_dir) / 'vitals_bmi.tsv'
    cfg = _make_config(tmp_path, {
        'table': str(table),
        'data_type': 'quantitative',
        'sample_id_col': 'PMBB_ID',
        'date_col': 'VISIT_DATE',
        'value_col': 'BMI_VALUE',
        'sep': 'tsv',
    })
    run_preprocess(cfg, str(tmp_path))
    result = pd.read_csv(tmp_path / 'TEST.long.tsv', sep='\t')
    assert set(result['sample_ID']) == {'P001', 'P002', 'P003', 'P004'}


def test_run_preprocess_sample_list_filters_coded(tmp_path, pmbb_test_dir):
    """sample_list is applied before coded code explosion."""
    table = Path(pmbb_test_dir) / 'icd.tsv'
    sample_list = tmp_path / 'samples.txt'
    # Only keep P001
    sample_list.write_text('P001\n')
    cfg = _make_config(tmp_path, {
        'table': str(table),
        'data_type': 'coded',
        'sample_id_col': 'PMBB_ID',
        'value_col': 'ICD_CODE',
        'date_col': 'ENCOUNTER_DATE',
        'sep': 'tsv',
        'filter': {'CODE_TYPE': ['ICD10']},
        'subsample': 'all',
    })
    run_preprocess(cfg, str(tmp_path), sample_list_path=str(sample_list))
    output_files = list(tmp_path.glob('*.long.tsv'))
    assert len(output_files) > 0, "Expected at least one .long.tsv output file"
    for f in output_files:
        df = pd.read_csv(f, sep='\t')
        assert set(df['sample_ID']) == {'P001'}, f"Found non-P001 samples in {f.name}"


def test_run_preprocess_sample_list_empty_file(tmp_path, pmbb_test_dir):
    """Empty sample_list produces empty output (no rows)."""
    table = Path(pmbb_test_dir) / 'vitals_bmi.tsv'
    sample_list = tmp_path / 'samples.txt'
    sample_list.write_text('')
    cfg = _make_config(tmp_path, {
        'table': str(table),
        'data_type': 'quantitative',
        'sample_id_col': 'PMBB_ID',
        'date_col': 'VISIT_DATE',
        'value_col': 'BMI_VALUE',
        'sep': 'tsv',
    })
    run_preprocess(cfg, str(tmp_path), sample_list_path=str(sample_list))
    result = pd.read_csv(tmp_path / 'TEST.long.tsv', sep='\t')
    assert len(result) == 0
