"""
Unit tests for bin/icd_utils.py

Coverage:
  - normalize_data_type()
  - code_matches_pattern()
  - build_prefix_index()
  - sample_has_any_pattern()
  - detect_data_type_from_values()
  - detect_data_type_from_column_name()
"""

import pytest
from icd_utils import (
    normalize_data_type,
    code_matches_pattern,
    build_prefix_index,
    sample_has_any_pattern,
    detect_data_type_from_values,
    detect_data_type_from_column_name,
    normalize_patterns,
)


# ---------------------------------------------------------------------------
# normalize_data_type
# ---------------------------------------------------------------------------

class TestNormalizeDataType:

    @pytest.mark.parametrize('raw,expected', [
        # ICD9 variants
        ('ICD09',     'ICD9'),
        ('ICD-9',     'ICD9'),
        ('ICD9CM',    'ICD9'),
        ('ICD-9-CM',  'ICD9'),
        ('icd9',      'ICD9'),
        ('icd09',     'ICD9'),
        ('icd-9',     'ICD9'),
        ('ICD_9',     'ICD9'),
        # ICD10 variants
        ('ICD-10',    'ICD10'),
        ('ICD10CM',   'ICD10'),
        ('ICD-10-CM', 'ICD10'),
        ('icd10',     'ICD10'),
        ('icd-10',    'ICD10'),
        ('ICD_10',    'ICD10'),
        # Already canonical — pass through unchanged
        ('ICD9',      'ICD9'),
        ('ICD10',     'ICD10'),
        # Other data_type aliases
        ('Quantitative', 'quantitative'),
        ('QUANT',        'quantitative'),
        ('Lab',          'measurements'),
        ('LAB',          'measurements'),
        ('Drug',         'drugs'),
        ('MEDICATION',   'drugs'),
        ('SNOMED_CT',    'SNOMED'),
        ('Demo',         'demographic'),
        ('DEMO',         'demographic'),
    ])
    def test_known_aliases(self, raw, expected):
        assert normalize_data_type(raw) == expected

    def test_unknown_data_type_passthrough(self):
        assert normalize_data_type('CUSTOM_CODE') == 'CUSTOM_CODE'

    def test_leading_trailing_whitespace_stripped(self):
        assert normalize_data_type('  ICD09  ') == 'ICD9'

    def test_non_string_integer(self):
        # Non-strings should not raise; returned as string
        result = normalize_data_type(42)
        assert isinstance(result, str)

    def test_empty_string_passthrough(self):
        # Empty string is not in aliases → returned as-is
        assert normalize_data_type('') == ''


# ---------------------------------------------------------------------------
# code_matches_pattern
# ---------------------------------------------------------------------------

class TestCodeMatchesPattern:

    # -- Exact match (no wildcards) -----------------------------------------

    def test_exact_match_true(self):
        assert code_matches_pattern('N80.0', 'N80.0') is True

    def test_exact_match_false_different_subcode(self):
        assert code_matches_pattern('N80.1', 'N80.0') is False

    def test_exact_match_false_prefix_only(self):
        # 'N80' without wildcard is exact — 'N80.0' should NOT match
        assert code_matches_pattern('N80.0', 'N80') is False

    def test_exact_match_icd9(self):
        assert code_matches_pattern('617', '617') is True

    def test_exact_match_icd9_subcoded(self):
        assert code_matches_pattern('617.1', '617.1') is True

    # -- Single wildcard (N80.*) --------------------------------------------

    def test_wildcard_matches_direct_child(self):
        assert code_matches_pattern('N80.0', 'N80.*') is True

    def test_wildcard_matches_deeper_child(self):
        assert code_matches_pattern('N80.11', 'N80.*') is True

    def test_wildcard_does_not_match_no_dot(self):
        # 'N80.*' prefix is 'N80.' — 'N80X0' has no dot after N80, so no match
        assert code_matches_pattern('N80X0', 'N80.*') is False

    def test_wildcard_does_not_match_adjacent_code(self):
        # N80.* should NOT match N81.0
        assert code_matches_pattern('N81.0', 'N80.*') is False

    def test_wildcard_icd9(self):
        assert code_matches_pattern('617.1', '617.*') is True
        assert code_matches_pattern('617.11', '617.*') is True

    def test_wildcard_icd9_no_match_different_code(self):
        assert code_matches_pattern('618.0', '617.*') is False

    # -- Double wildcard (I70.0**) ------------------------------------------

    def test_double_wildcard_matches_subcode(self):
        assert code_matches_pattern('I70.00', 'I70.0**') is True

    def test_double_wildcard_matches_multi_digit(self):
        assert code_matches_pattern('I70.09', 'I70.0**') is True

    def test_double_wildcard_does_not_match_sibling(self):
        # I70.0** prefix is 'I70.0' — I70.10 does NOT start with I70.0
        assert code_matches_pattern('I70.10', 'I70.0**') is False

    # -- Trailing dot pattern -----------------------------------------------

    def test_trailing_dot_matches_child(self):
        assert code_matches_pattern('I83.0', 'I83.') is True

    def test_trailing_dot_does_not_match_no_dot(self):
        # 'I83.' requires the dot — 'I830' must not match
        assert code_matches_pattern('I830', 'I83.') is False

    def test_trailing_dot_does_not_match_exact(self):
        # 'I83' itself does not start with 'I83.' (missing the dot)
        assert code_matches_pattern('I83', 'I83.') is False

    # -- Fringe / edge cases ------------------------------------------------

    def test_empty_pattern_returns_false(self):
        assert code_matches_pattern('N80.0', '') is False

    def test_empty_code_returns_false(self):
        assert code_matches_pattern('', 'N80.*') is False

    def test_both_empty_returns_false(self):
        assert code_matches_pattern('', '') is False

    def test_whitespace_stripped_from_inputs(self):
        assert code_matches_pattern('  N80.0  ', '  N80.0  ') is True

    def test_wildcard_only_pattern(self):
        # '.*' → prefix is '.', matches any code starting with '.'
        # This is a pathological pattern; verify it doesn't crash
        result = code_matches_pattern('N80.0', '.*')
        assert isinstance(result, bool)

    def test_ukbb_no_dot_icd10_exact(self):
        # UK Biobank encodes ICD10 without a dot: 'I251' not 'I25.1'
        assert code_matches_pattern('I251', 'I251') is True
        assert code_matches_pattern('I252', 'I251') is False

    def test_ukbb_no_dot_wildcard(self):
        # 'I25*' strip trailing '*' → prefix 'I25'; 'I251' starts with 'I25'
        assert code_matches_pattern('I251', 'I25*') is True
        assert code_matches_pattern('I251', 'I26*') is False

    def test_pattern_with_leading_wildcard_not_treated_as_prefix(self):
        # '*80.0' has '*' in it; rstrip('*') → '*80.0' (star not at end stays)
        # This documents current behavior — the leading '*' is preserved in prefix
        pat = '*80.0'
        result = code_matches_pattern('N80.0', pat)
        assert isinstance(result, bool)  # should not crash

    def test_z_code_icd10_wildcard(self):
        assert code_matches_pattern('Z34.10', 'Z34.*') is True
        assert code_matches_pattern('Z35.0',  'Z34.*') is False


# ---------------------------------------------------------------------------
# build_prefix_index
# ---------------------------------------------------------------------------

class TestBuildPrefixIndex:

    def test_basic_grouping(self):
        idx = build_prefix_index(['N80.0', 'N80.*', 'I70.0**'])
        assert 'N8' in idx
        assert 'I7' in idx
        assert 'N80.0' in idx['N8']
        assert 'N80.*' in idx['N8']
        assert 'I70.0**' in idx['I7']

    def test_empty_input(self):
        assert build_prefix_index([]) == {}

    def test_empty_strings_skipped(self):
        idx = build_prefix_index(['', '  ', 'N80.*'])
        assert '' not in idx
        assert '  ' not in idx
        assert 'N8' in idx

    def test_none_values_skipped(self):
        # None in the list should not crash
        idx = build_prefix_index([None, 'N80.0'])
        assert 'N8' in idx

    def test_single_char_pattern(self):
        # Pattern shorter than 2 chars — key is the pattern itself
        idx = build_prefix_index(['N'])
        assert 'N' in idx


# ---------------------------------------------------------------------------
# sample_has_any_pattern
# ---------------------------------------------------------------------------

class TestSampleHasAnyPattern:

    def test_match_found(self):
        assert sample_has_any_pattern({'N80.0', 'I10'}, ['N80.*']) is True

    def test_no_match(self):
        assert sample_has_any_pattern({'N81.0', 'I10'}, ['N80.*']) is False

    def test_empty_code_set(self):
        assert sample_has_any_pattern(set(), ['N80.*']) is False

    def test_empty_patterns(self):
        assert sample_has_any_pattern({'N80.0'}, []) is False

    def test_all_empty_string_patterns_skipped(self):
        assert sample_has_any_pattern({'N80.0'}, ['', ' ', None]) is False

    def test_multiple_patterns_one_matches(self):
        assert sample_has_any_pattern({'I70.01'}, ['N80.*', 'I70.0**']) is True

    def test_with_prefix_index(self):
        patterns = ['N80.*', 'I70.0**']
        idx = build_prefix_index(patterns)
        assert sample_has_any_pattern({'N80.5'}, patterns, prefix_index=idx) is True
        assert sample_has_any_pattern({'Z34.0'}, patterns, prefix_index=idx) is False

    def test_icd9_patterns(self):
        assert sample_has_any_pattern({'617', '617.1'}, ['617.*']) is True
        assert sample_has_any_pattern({'618.0'}, ['617.*']) is False

    def test_none_code_in_set_skipped(self):
        # None in set should not raise
        assert sample_has_any_pattern({None, 'N80.0'}, ['N80.*']) is True

    def test_exact_match_in_large_set(self):
        codes = {f'N{i:02d}.0' for i in range(1, 99)}
        codes.add('N80.1')
        assert sample_has_any_pattern(codes, ['N80.*']) is True


# ---------------------------------------------------------------------------
# detect_data_type_from_values
# ---------------------------------------------------------------------------

class TestDetectDataTypeFromValues:

    def test_icd10_values(self):
        values = ['N80.0', 'N80.1', 'I10', 'Z34.10']
        assert detect_data_type_from_values(values) == 'ICD10'

    def test_icd9_values(self):
        values = ['617', '617.1', '250.00', '714.0']
        assert detect_data_type_from_values(values) == 'ICD9'

    def test_quantitative_float_values(self):
        # All parseable as numbers → quantitative
        values = ['27.8', '28.4', '29.1', '31.7']
        assert detect_data_type_from_values(values) == 'quantitative'

    def test_quantitative_integer_values(self):
        values = ['1', '0', '1', '0', '1']
        assert detect_data_type_from_values(values) == 'quantitative'

    def test_empty_values_returns_unknown(self):
        assert detect_data_type_from_values([]) == 'unknown'

    def test_all_nan_returns_unknown(self):
        import pandas as pd
        assert detect_data_type_from_values([None, None, float('nan')]) == 'unknown'

    def test_categorical_values(self):
        values = ['male', 'female', 'unknown', 'male', 'female']
        result = detect_data_type_from_values(values)
        assert result == 'categorical'

    def test_majority_icd10(self):
        # > 70% ICD10 → ICD10, even with a stray string
        values = ['N80.0', 'N80.1', 'I10', 'N80.2', 'N80.3', 'N80.4', 'I25.1', 'UNKNOWN', 'N80.5']
        result = detect_data_type_from_values(values)
        assert result == 'ICD10'


# ---------------------------------------------------------------------------
# detect_data_type_from_column_name
# ---------------------------------------------------------------------------

class TestDetectDataTypeFromColumnName:

    @pytest.mark.parametrize('col,expected', [
        ('N80.0',      'ICD10'),   # starts with letter+digits
        ('I70.0',      'ICD10'),
        ('617',        'ICD9'),    # starts with digits
        ('250.00',     'ICD9'),
        ('BMI',        'quantitative'),
        ('HEIGHT',     'quantitative'),
        ('WEIGHT',     'quantitative'),
        ('SEX',        'demographic'),
        ('GENDER',     'demographic'),
        ('ANCESTRY',   'demographic'),
        ('RACE',       'demographic'),
        ('ICD10_col',  'ICD10'),
        ('icd10_N80',  'ICD10'),
        ('ICD9_code',  'ICD9'),
        ('ICD09_code', 'ICD9'),
    ])
    def test_known_column_names(self, col, expected):
        assert detect_data_type_from_column_name(col) == expected

    def test_unknown_column_returns_none(self):
        assert detect_data_type_from_column_name('SOME_RANDOM_COL') is None

    def test_empty_string_returns_none(self):
        assert detect_data_type_from_column_name('') is None

    def test_case_insensitive_keyword_matching(self):
        assert detect_data_type_from_column_name('bmi') == 'quantitative'
        assert detect_data_type_from_column_name('sex') == 'demographic'


# ---------------------------------------------------------------------------
# normalize_patterns
# ---------------------------------------------------------------------------

class TestNormalizePatterns:

    def test_bare_code_gets_wildcard(self):
        assert normalize_patterns(['E11']) == ['E11.*']

    def test_dotted_code_unchanged(self):
        assert normalize_patterns(['E10.5']) == ['E10.5']

    def test_explicit_wildcard_unchanged(self):
        assert normalize_patterns(['E11.*']) == ['E11.*']

    def test_empty_list_returns_empty(self):
        assert normalize_patterns([]) == []

    def test_none_list_returns_empty(self):
        assert normalize_patterns(None) == []

    def test_icd9_bare_code_gets_wildcard(self):
        assert normalize_patterns(['250']) == ['250.*']

    def test_empty_string_elements_filtered(self):
        assert normalize_patterns(['', ' ', 'E11']) == ['E11.*']

    def test_prefix_matching_false_bare_code_not_expanded(self):
        assert normalize_patterns(['E11'], prefix_matching=False) == ['E11']

    def test_prefix_matching_false_icd9_bare_code_not_expanded(self):
        assert normalize_patterns(['250'], prefix_matching=False) == ['250']

    def test_prefix_matching_false_explicit_wildcard_unchanged(self):
        assert normalize_patterns(['E11.*'], prefix_matching=False) == ['E11.*']

    def test_prefix_matching_false_dotted_code_unchanged(self):
        assert normalize_patterns(['E11.9'], prefix_matching=False) == ['E11.9']

    def test_prefix_matching_false_mixed_list(self):
        assert normalize_patterns(['E11', 'E11.*', 'E10.5'], prefix_matching=False) == ['E11', 'E11.*', 'E10.5']

    def test_prefix_matching_true_explicit_same_as_default(self):
        assert normalize_patterns(['E11'], prefix_matching=True) == ['E11.*']
