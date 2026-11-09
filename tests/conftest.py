"""
Shared pytest fixtures for pmbb-nf-toolkit-phenotyping tests.
"""

import sys
import os
import pytest
import pandas as pd

# Ensure bin/ is importable regardless of cwd
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_data')
sys.path.insert(0, os.path.abspath(BIN_DIR))


@pytest.fixture(scope='session')
def test_data_dir():
    return os.path.abspath(TEST_DATA_DIR)



@pytest.fixture
def minimal_long_df():
    """Minimal canonical long-format DataFrame for unit testing."""
    return pd.DataFrame({
        'sample_ID': ['P001', 'P001', 'P002', 'P003'],
        'data_type':   ['ICD10', 'quantitative', 'ICD10', 'quantitative'],
        'phenotype': ['N80.0', 'BMI', 'N80.1', 'BMI'],
        'value':     ['N80.0', '28.4', 'N80.1', '31.7'],
        'occurrence_date': ['2020-01-01', '2020-01-01', '2021-03-15', '2021-03-15'],
    })


@pytest.fixture
def quant_df():
    """Small quantitative DataFrame with a clear outlier for filtering tests."""
    return pd.DataFrame({
        'sample_ID': [f'P{i:03d}' for i in range(1, 11)],
        'data_type':   ['quantitative'] * 10,
        'phenotype': ['ALT'] * 10,
        'value':     ['20.0', '25.0', '22.0', '18.0', '30.0',
                      '26.0', '24.0', '19.0', '21.0', '2500.0'],
        'occurrence_date': ['2020-01-01'] * 10,
    })


PMBB_TEST_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_data', 'PMBB')


@pytest.fixture(scope='session')
def pmbb_test_dir():
    return os.path.abspath(PMBB_TEST_DIR)


@pytest.fixture(scope='session')
def pmbb_test_config(pmbb_test_dir):
    return os.path.join(pmbb_test_dir, 'pmbb_test_config.yaml')


@pytest.fixture(scope='session')
def pmbb_vitals_bmi(pmbb_test_dir):
    return os.path.join(pmbb_test_dir, 'vitals_bmi.tsv')


@pytest.fixture(scope='session')
def pmbb_lipids(pmbb_test_dir):
    return os.path.join(pmbb_test_dir, 'labs_lipids.tsv')


@pytest.fixture(scope='session')
def pmbb_covariates(pmbb_test_dir):
    return os.path.join(pmbb_test_dir, 'covariates.tsv')


@pytest.fixture(scope='session')
def pmbb_icd(pmbb_test_dir):
    return os.path.join(pmbb_test_dir, 'icd.tsv')
