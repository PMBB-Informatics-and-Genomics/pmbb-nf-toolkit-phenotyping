#!/usr/bin/env python3
"""
Shared utilities for ICD code matching, data_type normalization, and auto-detection.
Used by preprocess.py and all downstream phenotyping modules.
"""

import re

# ---------------------------------------------------------------------------
# DataType normalization
# ---------------------------------------------------------------------------

DATA_TYPE_ALIASES = {
    # ICD9 variants
    'ICD09': 'ICD9',
    'ICD-9': 'ICD9',
    'ICD9CM': 'ICD9',
    'ICD-9-CM': 'ICD9',
    'icd9': 'ICD9',
    'icd09': 'ICD9',
    'icd-9': 'ICD9',
    'ICD_9': 'ICD9',
    # ICD10 variants
    'ICD-10': 'ICD10',
    'ICD10CM': 'ICD10',
    'ICD-10-CM': 'ICD10',
    'icd10': 'ICD10',
    'icd-10': 'ICD10',
    'ICD_10': 'ICD10',
    # Other standardizations
    'SNOMED_CT': 'SNOMED',
    'snomedct': 'SNOMED',
    'snomed': 'SNOMED',
    'Quantitative': 'quantitative',
    'QUANT': 'quantitative',
    'quant': 'quantitative',
    'Categorical': 'categorical',
    'CAT': 'categorical',
    'Binary': 'binary',
    'BINARY': 'binary',
    'Demo': 'demographic',
    'DEMO': 'demographic',
    'Demographics': 'demographic',
    'Vital': 'vitals',
    'VITAL': 'vitals',
    'Drug': 'drugs',
    'DRUG': 'drugs',
    'Medication': 'drugs',
    'MEDICATION': 'drugs',
    'Lab': 'measurements',
    'LAB': 'measurements',
    'Laboratory': 'measurements',
    'Measurement': 'measurements',
}

def normalize_data_type(data_type: str) -> str:
    """
    Normalize data_type labels to canonical form.

    Examples:
        'ICD09'   -> 'ICD9'
        'icd10'   -> 'ICD10'
        'ICD-10'  -> 'ICD10'
        'ICD10CM' -> 'ICD10'

    Unknown data_types are returned as-is (uppercased).
    """
    if not isinstance(data_type, str):
        return str(data_type)
    data_type = data_type.strip()
    return DATA_TYPE_ALIASES.get(data_type, data_type)


# ---------------------------------------------------------------------------
# ICD pattern matching
# ---------------------------------------------------------------------------
# NOTE: We use prefix/exact matching (str.startswith) rather than regex.
# The old Snakemake approach of `code.replace('*', '.')` was a bug: in regex,
# '.' matches ANY character, so 'N80.*' would incorrectly match 'N80X0'.
# Our approach: strip trailing wildcards/dots → prefix match.

def code_matches_pattern(patient_code: str, pattern: str) -> bool:
    """
    Test whether a patient ICD code matches a phenotype pattern.

    Wildcard rules:
        'N80.0'   -> exact match only
        'N80.*'   -> prefix match on 'N80.' (startswith)
        'N80.0**' -> prefix match on 'N80.0' (double-star same as single)
        'I70.0**' -> prefix match on 'I70.0'
        'I83.'    -> prefix match on 'I83' (trailing dot stripped)
        'N80'     -> prefix match on 'N80' (no dot = any child code)
    """
    if not pattern or not patient_code:
        return False
    pattern = pattern.strip()
    patient_code = patient_code.strip()

    if '*' in pattern:
        # Strip only trailing wildcards — the dot is part of the prefix (ICD decimal separator).
        # e.g. 'N80.*'   → prefix 'N80.'  (matches N80.0, N80.1, but NOT N80X0)
        #      'I70.0**' → prefix 'I70.0' (matches I70.00, I70.01, but NOT I70.10)
        prefix = pattern.rstrip('*')
        return patient_code.startswith(prefix)
    elif pattern.endswith('.'):
        # Trailing dot → prefix match including the dot, so 'N80.' matches 'N80.0'
        # but NOT 'N800' (the dot is the ICD decimal separator and is required).
        return patient_code.startswith(pattern)
    else:
        return patient_code == pattern


def build_prefix_index(patterns: list) -> dict:
    """
    Build a dict mapping 2-char prefix → list of patterns with that prefix.
    Used to skip irrelevant patterns quickly before full matching.

        {'N8': ['N80.0', 'N80.*', 'N80.1'], 'I7': ['I70.0**'], ...}
    """
    idx: dict = {}
    for p in patterns:
        if p and p.strip():
            key = p.strip()[:2]
            idx.setdefault(key, []).append(p.strip())
    return idx


def sample_has_any_pattern(
    code_set: set,
    patterns: list,
    prefix_index: dict | None = None,
) -> bool:
    """
    Return True if any code in code_set matches any pattern in patterns.

    Uses prefix_index for fast candidate filtering when provided.
    """
    if not patterns or not code_set:
        return False
    # Filter out empty/None patterns
    patterns = [p for p in patterns if p and str(p).strip()]
    if not patterns:
        return False

    if prefix_index is None:
        prefix_index = build_prefix_index(patterns)

    for code in code_set:
        if not code:
            continue
        key = str(code).strip()[:2]
        if key in prefix_index:
            for pat in prefix_index[key]:
                if code_matches_pattern(str(code).strip(), pat):
                    return True
    return False


def normalize_patterns(patterns: list | None, prefix_matching: bool = True) -> list:
    """
    Normalize ICD pattern lists for prefix matching.

    Bare codes (no dot, no wildcard) get '.*' appended when prefix_matching=True
    so they prefix-match all subcodes: 'E11' → 'E11.*' matches 'E11.9', 'E11.65'.
    Dotted codes ('E10.5') and patterns with existing wildcards are left as-is.
    When prefix_matching=False, bare codes are kept as exact-match patterns.
    """
    if not patterns:
        return []
    normalized = []
    for p in patterns:
        p = str(p).strip()
        if p and prefix_matching and '*' not in p and not p.endswith('.') and '.' not in p:
            p = p + '.*'
        if p:
            normalized.append(p)
    return normalized


# ---------------------------------------------------------------------------
# DataType auto-detection
# ---------------------------------------------------------------------------

# ICD10: letter followed by 2+ digits (optionally with decimal)
_ICD10_VALUE_RE = re.compile(r'^[A-Z]\d{2}(\.\d*)?$', re.IGNORECASE)
# ICD9: 3-5 digits optionally with decimal (e.g. 274.0, 617.1, 714)
_ICD9_VALUE_RE = re.compile(r'^\d{3,5}(\.\d+)?$')
# ICD10 column name pattern (e.g. "N80.0", "ICD10_N80.0", columns starting with letter+digits)
_ICD10_COL_RE = re.compile(r'^[A-Z]\d{2}', re.IGNORECASE)
_ICD9_COL_RE = re.compile(r'^\d{3}')


def detect_data_type_from_values(values, sample_size: int = 50) -> str:
    """
    Heuristically detect data_type type from a sample of values.

    Returns one of: 'ICD10', 'ICD9', 'quantitative', 'categorical', 'unknown'
    """
    import pandas as pd
    sample = pd.Series(values).dropna().astype(str).head(sample_size)
    if len(sample) == 0:
        return 'unknown'

    icd10_matches = sample.str.match(r'^[A-Z]\d{2}', na=False).sum()
    icd9_matches = sample.str.match(r'^\d{3,5}(\.\d+)?$', na=False).sum()

    frac = len(sample)
    if icd10_matches / frac > 0.7:
        return 'ICD10'
    if icd9_matches / frac > 0.7:
        return 'ICD9'

    # Try numeric (quantitative) — after ICD checks so numeric ICD9 codes
    # (e.g. 617, 250.00) are not misclassified as quantitative.
    try:
        pd.to_numeric(sample, errors='raise')
        return 'quantitative'
    except (ValueError, TypeError):
        pass

    return 'categorical'


def detect_data_type_from_column_name(col_name: str) -> str | None:
    """
    Heuristically detect data_type type from a column name.
    Returns None if no confident detection.
    """
    col_upper = col_name.upper()
    if 'ICD9' in col_upper or 'ICD09' in col_upper:
        return 'ICD9'
    if 'ICD10' in col_upper or 'ICD-10' in col_upper:
        return 'ICD10'
    if 'SNOMED' in col_upper:
        return 'SNOMED'
    if 'BMI' in col_upper or 'HEIGHT' in col_upper or 'WEIGHT' in col_upper:
        return 'quantitative'
    if 'SEX' in col_upper or 'GENDER' in col_upper or 'ANCESTRY' in col_upper or 'RACE' in col_upper:
        return 'demographic'
    if _ICD10_COL_RE.match(col_name):
        return 'ICD10'
    if _ICD9_COL_RE.match(col_name):
        return 'ICD9'
    return None
