# pmbb-nf-toolkit-phenotyping

Modular Nextflow DSL2 phenotyping pipeline for biobanks (Penn Medicine Biobank, All of Us, UK Biobank, and others).

## Overview

The pipeline takes a single biobank YAML config describing all phenotypes and their source tables, resolves it into per-phenotype JSON configs, then preprocesses and phenotypes each one in parallel. An optional final GATHER step joins all outputs into a wide or long table.

The canonical intermediate format is a long-format TSV:

```
sample_ID | data_type | phenotype | value | occurrence_date
```

### Modules

| Module | Status | Description |
|--------|--------|-------------|
| `RESOLVE_CONFIG` | ✅ Implemented | Reads biobank YAML, resolves inheritance, writes per-phenotype JSONs |
| `PREPROCESS` | ✅ Implemented | Per-phenotype: read table, apply filters, outlier filtering, scaling; coded values explode into individual TSVs (basic) or one combined TSV (advanced) |
| `BASIC_CODED` | ✅ Implemented | Binary presence/absence per code value; works with any vocabulary (ICD, CPT, SNOMED, LOINC, phecodes, etc.); min_occurrences threshold; optional `column_prefix` |
| `ADVANCED_CODED` | ✅ Implemented | Named case/control phenotype from code patterns across any vocabulary; case_codes, case_exclude, control_exclude |
| `BASIC_QUANTITATIVE` | ✅ Implemented | Per-sample summary statistics (mean, median, std, etc.), quantile binning, min_occurrences |
| `BASIC_CATEGORICAL` | ✅ Implemented | Auto (default), dictionary, one-hot, and binary encoding; min_occurrences |
| `GATHER` | ✅ Implemented | Outer-join all per-phenotype TSVs into wide or long table |
| Phenotype manifest | ✅ Implemented | `phenotype_manifest.tsv`: one row per output column (column, phenotype, data_type, raw code) for subsetting the gathered table — see [Outputs](#outputs) |
| `EVENT_MASK` | ✅ Implemented | Compute per-sample date-range windows from named clinical events (ICD codes + offsets) |
| `APPLY_MASK` | ✅ Implemented | Filter opted-in phenotype long.tsvs against mask windows before downstream calculation |
| `RISK_CLASSIFICATION` | ✅ Implemented | Percentile/threshold/quantile_bin risk categorization on gathered wide output |
| `EDA_STATS` / `EDA_REPORT` / `EDA_PROFILE` | ✅ Implemented | Optional QC report (`--run_eda`): per-column metrics and flags for every phenotype (separation, missingness, low cases), data-driven restriction suggestions against any covariate, bounded ydata profiles — see [EDA Report](#eda-report---run_eda) |
| `FILTER_PHENOTYPES` | ✅ Implemented | Optional post-processing (`--filter_config`): drops / sex- or cohort-restricts columns into new `filtered/` files; raw outputs untouched — see [Filtering](#filtering---filter_config) |

---

## Quick Start

### Requirements

- [Nextflow](https://nextflow.io/) ≥ 23.04
- Python ≥ 3.10 with pandas, numpy, PyYAML (or use a conda env / container)

### Run with a biobank config

```bash
nextflow run main.nf \
  --biobank_config biobanks/PMBB/pmbb_config.yaml \
  --output_dir results/
```

### Run with gather (produces a single wide-format table)

```bash
nextflow run main.nf \
  --biobank_config biobanks/PMBB/pmbb_config.yaml \
  --output_dir results/ \
  --run_gather true \
  --gather_format wide
```

### Run test profile

```bash
nextflow run main.nf -profile test
```

This uses `test_data/PMBB/pmbb_test_config.yaml` with small synthetic fixtures.

To also exercise the EDA report and the filter step:

```bash
nextflow run main.nf -profile test --run_eda \
    --sample_list test_data/PMBB/filter/roster.txt \
    --filter_config test_data/PMBB/pmbb_test_filter.yaml
```

`test_data/PMBB/filter/` holds frozen EDA outputs from that run, so
`bin/filter_phenotypes.py` can also be run standalone (command in the header of
`pmbb_test_filter.yaml`).

---

## Pipeline Topology

```
biobank_config.yaml
       │
       ▼
  RESOLVE_CONFIG (runs once)
  ├── resolves global → data_type → phenotype inheritance
  ├── validates required fields
  └── writes one JSON per phenotype → output_dir/configs/
       │
       ├─ quantitative JSONs ─┐
       ├─ categorical JSONs   ├─ all fan out to PREPROCESS in parallel
       └─ coded JSONs  ──┘
       │
       ▼
  PREPROCESS (per phenotype — parallel)
  ├── if preprocessed_path → pass file through (skip processing)
  ├── reads raw table, applies row-level filter (glob wildcards)
  ├── quantitative: outlier filtering + optional scaling → single .long.tsv
  ├── categorical: → single .long.tsv
  ├── coded (basic):    scans value_col, matches subsample → one .long.tsv per code
  └── coded (advanced): case_codes present → one combined .long.tsv for all codes
       │
       ├──────────────────┬──────────────────┬──────────────────┐
       ▼                  ▼                  ▼                  ▼
  BASIC_CODED  ADVANCED_CODED  BASIC_QUANT   BASIC_CATEGORICAL
  (per code)        (per named phenotype) (per phenotype) (per phenotype)
  └── column = {column_prefix}{code}
       │                  │                  │                  │
       │   each also writes manifest rows (phenotype_manifest.py)
       └──────────────────┴──────────────────┴──────────────────┘
                          │
                          ├──► phenotype_manifest.tsv (always; one row per output column)
                          ▼
                    GATHER (optional)
                    ├── outer-join all TSVs on sample_ID
                    └── output: phenotype_table.wide.tsv / .long.tsv
                          │
                          ▼
                    RISK_CLASSIFICATION (optional, requires GATHER wide output)
                    ├── reads risk rules from risk_classification_config.json (emitted by RESOLVE_CONFIG)
                    ├── applies percentile / threshold / quantile_bin rules per phenotype
                    ├── output: risk_classification/risk_matrix.tsv (sample_ID + one risk column per rule)
                    ├── output: risk_classification/risk_report.txt (computed cutoffs and class counts)
                    └── NOT merged into phenotype_table.wide.tsv; NOT listed in phenotype_manifest.tsv

per-phenotype results ──► EDA_STATS (optional, --run_eda; batches of --eda_batch_size files)
                          ├── per-column metrics vs. eda.covariates (separation, missingness)
                          ▼
                    EDA_REPORT ──► flags, eda_report.html, restriction suggestions,
                          │        profile_selection.tsv (bounded subset)
                          ▼
                    EDA_PROFILE (one ydata report per selected group, parallel)

GATHER wide + EDA tables ──► FILTER_PHENOTYPES (optional, --filter_config)
                             └── filtered/phenotype_table.filtered.wide.tsv + filter_log.tsv
```

---

## Outputs

```
{output_dir}/
├── phenotype_table.wide.tsv      # GATHER (--run_gather): sample_ID + every phenotype column
├── phenotype_table.long.tsv      # GATHER with --gather_format long|both
├── phenotype_manifest.tsv        # always: one row per phenotype output column (see below)
├── phenotypes/                   # one TSV per phenotype (*.quant.tsv, *.cat.tsv, *.diag.tsv) + reports;
│                                 #   *.codes.tsv = value→code map for auto-encoded categoricals
├── long_format/                  # PREPROCESS long-format TSVs (raw codes) + resolved JSONs
├── configs/                      # RESOLVE_CONFIG per-phenotype JSONs — the effective settings each job used
├── risk_classification/          # only with --run_risk_classification (separate — see below)
│   ├── risk_matrix.tsv
│   └── risk_report.txt
├── eda_report.zip                # only with --run_eda: the whole eda/ folder, for sharing (overwritten each run)
├── eda/                          # only with --run_eda (see EDA Report below)
│   ├── eda_report.html
│   ├── run_config.html / .json   #   run parameters + resolved global / data_type / eda settings
│   ├── biobank_config.yaml / .html #  verbatim copy of the biobank YAML + browser view
│   ├── eda_summary.tsv / eda_flags.tsv / restriction_suggestions.tsv / profile_selection.tsv
│   ├── filter_suggested.yaml     #   editable starting point for --filter_config
│   ├── eda_covariates.tsv        #   covariate levels used (read by the filter)
│   ├── pairwise_overlap.tsv      #   only with eda.pairwise: true
│   └── profiles/<group>.html
├── filtered/                     # only with --filter_config (raw outputs are untouched)
│   ├── phenotype_table.filtered.wide.tsv
│   ├── phenotype_manifest.filtered.tsv
│   ├── filter_log.tsv
│   └── filter_config_used.yaml
├── pipeline_report.html
└── pipeline_timeline.html
```

### `phenotype_manifest.tsv`

A tab-separated index of **every column the pipeline produced** — the same columns as the header of `phenotype_table.wide.tsv` (minus `sample_ID`). It is meant for the "run once, build one giant table, subset later" workflow: pull out all ICD10 codes, all PheCodes, all quantitative phenotypes, etc. without guessing from column names.

| Column | Contents |
|---|---|
| `column_name` | Exact column header in `phenotype_table.wide.tsv` / the per-phenotype TSV (includes any `column_prefix`, stat suffix, or one-hot suffix). |
| `phenotype` | The phenotype's key in the YAML config (e.g. `ICD10`, `PheCode`, `BMI`, `T2Diab`). For chunked coded runs this is still the original key, not `ICD10_chunk_001`. Not affected by `output_name`. |
| `data_type` | `quantitative`, `categorical`, `coded`, or `advanced_coded`. |
| `code` | `coded` rows only: the **raw code** as it appears in the input table — `column_name` with `column_prefix` removed. Blank for other data types. |
| `warnings` | `advanced_coded` rows only: config problems found while phenotyping, `; `-separated — patterns that matched no code in the input data (with a hint when similar codes exist) and `0 cases`. Blank when there are none. See [Code matching rules](#code-matching-rules). Find them with e.g. `manifest[manifest.warnings != '']`. |
| `flags` | Only with `--run_eda`: comma-separated EDA flags for the column (e.g. `low_cases,high_missing`), blank if none. See [EDA Report](#eda-report---run_eda). Filter with e.g. `manifest[manifest.flags == '']`. |

Example:

```
column_name     phenotype  data_type       code    warnings
T2Diab          T2Diab     advanced_coded
Cirrhosis       Cirrhosis  advanced_coded          case_codes pattern 'K74.6' matched no codes (data has K74.60, K74.69); 0 cases
ANCESTRY_AFR    ANCESTRY   categorical
SEX             SEX        categorical
icd10_E11.9     ICD10      coded           E11.9
phe_401.1       PheCode    coded           401.1
BMI_mean        BMI        quantitative
BMI_std         BMI        quantitative
```

Rows are sorted by `data_type`, then `phenotype`, then column; the columns of one phenotype keep their output order (e.g. `BMI_mean` before `BMI_std`). The manifest is written whether or not `--run_gather` is on — without GATHER it indexes the per-phenotype files in `phenotypes/`.

**Subsetting the gathered table:**

```bash
# column names for all ICD10 codes / all quantitative phenotypes
awk -F'\t' '$2=="ICD10"{print $1}'        phenotype_manifest.tsv > icd10_cols.txt
awk -F'\t' '$3=="quantitative"{print $1}' phenotype_manifest.tsv > quant_cols.txt

# detect colliding coded columns (same column_name from two sources)
tail -n +2 phenotype_manifest.tsv | cut -f1 | sort | uniq -d
```

```python
import pandas as pd
manifest = pd.read_csv('phenotype_manifest.tsv', sep='\t', dtype=str, keep_default_na=False)
cols = manifest.loc[manifest['phenotype'] == 'PheCode', 'column_name'].tolist()
phecodes = pd.read_csv('phenotype_table.wide.tsv', sep='\t', usecols=['sample_ID'] + cols)
```

**Not included:** risk classification columns (e.g. `BMI_risk`). They live only in `risk_classification/risk_matrix.tsv`, which is **not** merged into `phenotype_table.wide.tsv`; join it on `sample_ID` yourself if you want them together. See [Risk Classification → Outputs](#outputs-1).

---

## Biobank Config YAML

Each biobank run uses one YAML config. Parameters cascade: `global:` → `data_types.<type>:` → individual phenotype block. Later levels override earlier ones.

See [`examples/basic_config_example.yaml`](examples/basic_config_example.yaml) for a full annotated example covering quantitative, categorical, and coded phenotypes.

### Field Reference

| Field | Applies to | Description |
|---|---|---|
| `table` | all | Path to raw source table (absolute, or relative to the YAML file) |
| `data_type` | all | `quantitative`, `categorical`, `coded`, or `advanced_coded` |
| `value_col` | all | Column containing the measurement or code value |
| `sample_id_col` | all | Sample ID column (inheritable from global) |
| `date_col` | all | Date column (inheritable from global). Required for the `first` / `last` stats. |
| `sep` | all | Delimiter for reading the raw table. Named aliases: `tsv`, `tab`, `csv`, `comma`, `pipe`, `space`, `whitespace`. Any single character also accepted (quote in YAML: `sep: "\|"`). Auto-detected from extension when omitted: `.tsv`→tab, `.csv`→comma, `.txt`→tab, other→comma. |
| `filter` | all | Dict of `column: [values]`; row-level AND filter; glob wildcards in values |
| `subsample` | coded (basic) | `all` or list of code values (wildcards ok); each match → separate phenotype column |
| `column_prefix` | coded (basic) | String prepended to each exploded output column and its `.diag.tsv` filename (e.g. `"phe_"` → `phe_250.2`). Default `""` (raw code). **Output naming only** — does not change long-format files, `advanced_coded` matching, or event masking. Not allowed on `advanced_coded`. See [Coded Column Prefixes](#coded-column-prefixes-and-multiple-vocabularies). |
| `case_codes` | coded | Code patterns that make a patient a **case** (1). Works with any vocabulary — ICD, CPT, SNOMED, LOINC, phecodes, etc. Bare codes match the code itself and its subcodes by default (`E11` matches `E11`, `E11.9`, `E11.65`); disable with `prefix_matching: false` for exact-match systems like CPT. Wildcards (`E11.**`) also work. See [Code matching rules](#code-matching-rules). |
| `case_exclude` | coded | Code patterns that **exclude a patient (→ NA)** — overrides both case and control |
| `control_exclude` | coded | Code patterns that **exclude a non-case from being a control** (→ NA) |
| `prefix_matching` | coded, advanced_coded | `true` (default): bare codes match themselves and all sub-codes (e.g. `E11` → `E11`, `E11.9`, `E11.65`). `false`: bare codes require an exact match — use for non-hierarchical systems like CPT or phecodes. Settable at phenotype level or per source in `advanced_coded` (source-level wins). |
| `sources` | advanced_coded | Dict of named sources. Keys are arbitrary labels (e.g. `ICD10`, `CPT`); each entry specifies which `coded` phenotype to read and what code patterns to apply. Sources can mix vocabularies. Patterns from all sources are pooled before case/control assignment. |
| `from_phenotype` | advanced_coded (per source) | Name of a `coded` phenotype defined in this YAML whose preprocessed output feeds this source. Required. |
| `case_codes` | advanced_coded (per source) | Code patterns that make a patient a **case**. Required. Bare codes match themselves and their subcodes by default; set `prefix_matching: false` on the source for exact matching (e.g. CPT). See [Code matching rules](#code-matching-rules). |
| `case_exclude` | advanced_coded (per source) | Code patterns that **exclude a patient (→ NA)**, overriding case and control. |
| `control_exclude` | advanced_coded (per source) | Code patterns that **exclude a non-case from being a control** (→ NA). |
| `preprocessed_path` | all | Skip PREPROCESS; pass this long-format TSV directly downstream |
| `output_name` | all | Override the phenotype name in output columns |
| `min_occurrences` | all | Minimum observations per sample; below threshold → NA at phenotyping step. For `advanced_coded` it counts raw rows per exact code and also applies to exclusion codes — see [`min_occurrences` behavior](#min_occurrences-behavior). |
| `stats` | quantitative | Summary stats, one `{output_name}_{stat}` column each: `mean`, `median`, `std`, `min`, `max`, `count`, `squared` (the **mean**, squared), `first` / `last` (value at the earliest / latest date — requires `date_col`; see [`first` / `last` stats](#first--last-stats)). An unknown name (e.g. a typo) stops the run at RESOLVE_CONFIG with an error listing the valid stats. |
| `qcut_bins` | quantitative | Number of quantile bins (0 = continuous) |
| `scale` | quantitative | Whether to scale values |
| `scale_method` | quantitative | `zscore` or `minmax` |
| `outlier_method` | quantitative | `iqr`, `zscore`, or `none` |
| `outlier_mode` | quantitative | `cap` or `remove` |
| `outlier_iqr_multiplier` | quantitative | IQR multiplier (default 1.5) |
| `outlier_zscore_sd` | quantitative | Z-score SD threshold (default 3.0) |
| `dictionary` | categorical | Map raw string values to numeric codes; string `"NA"` → `NaN`. **Encoding precedence 1** — wins over `one_hot` and `binary`. See [Categorical Encoding](#categorical-encoding). |
| `one_hot` | categorical | `true` = expand categories into individual 0/1 columns. **Encoding precedence 2** — ignored if a `dictionary` is set. |
| `binary` | categorical | `true` = 1 if the sample has ≥ `min_occurrences` rows, regardless of value. **Encoding precedence 3** — ignored if `dictionary` or `one_hot: true` is set. |
| `categories` | categorical | Ordered list of valid categories. `one_hot`: which columns to create. Auto encoding: sets the code order (1, 2, …); values not listed → NA. |
| `missing_as_control` | categorical, coded | `true` = absent sample → 0; `false` = absent → NA. For coded/advanced_coded, "absent" means in `--sample_list` but without the code — without `--sample_list` there is no roster and no 0s are written. |
| `eda` | all | Whether this phenotype's columns get a ydata profile under `--run_eda`: `true`, `false`, or a list of code / column globs (`["E11*", I10]`). Cascades global → data_types → phenotype. Built-in defaults: `true` for quantitative/categorical, `false` for coded/advanced_coded. Metrics and flags are computed for **every** column regardless. See [EDA Report](#eda-report---run_eda). |
| `restrict_to` | all | A priori restriction `{covariate: level}`, e.g. `{sex: female}`; the covariate must be in `eda.covariates` (or `covariates_file`). Never applied to the pipeline outputs: it is listed in `eda/restriction_suggestions.tsv` and written as an active entry in `eda/filter_suggested.yaml`. |
| `risk_classification` | quantitative | Risk classification rule for this phenotype. Set at `data_types.quantitative` level as a default; override per phenotype or opt out with `false`. See [Risk Classification](#risk-classification) below. |

### Categorical Encoding

Each categorical phenotype uses exactly **one** encoding. When more than one is set, the **first match wins**:

| Precedence | Setting | Encoding | Output |
|---|---|---|---|
| 1 | `dictionary: {...}` | dictionary | One column; each value mapped through the dictionary (unlisted values → NA) |
| 2 | `one_hot: true` | one_hot | One 0/1 column per category (`{output_name}_{category}`) |
| 3 | `binary: true` | binary | One 0/1 column: does the sample have any row? The value itself is ignored. |
| 4 | *(none of the above)* | **auto** (default) | One column. If every value is numeric, values pass through unchanged. Otherwise each distinct value gets an integer code 1..N — natural-sorted (`Freeze 2` before `Freeze 10`) or in `categories` order — and the mapping is written to `phenotypes/{name}.codes.tsv`. |

For all encodings, a sample with several rows takes its most common value (ties → the lowest code), and a sample absent from the table is NA (or 0 with `missing_as_control: true`).

> **Check your global settings.** Options cascade from `global:` and `data_types.categorical:` into every categorical phenotype. A `one_hot: true` or `binary: true` set at those levels applies to **every** categorical phenotype that doesn't have a `dictionary` — override it per phenotype with `one_hot: false` / `binary: false`.
>
> To see what was actually applied, open `{output_dir}/configs/categorical/{name}.json`. Each file holds the fully merged settings for that phenotype, plus an `encoding` field (`dictionary`, `one_hot`, `binary`, or `auto`) recording which encoding ran.

### Advanced Coded Phenotypes

`data_type: advanced_coded` produces a single named case/control output column per phenotype (e.g. `T2Diab: 1/0/NA`) built from code patterns across one or more source tables. Use it instead of `data_type: coded` when you want a named phenotype rather than one column per code value. Sources can mix vocabularies — e.g. combine an ICD10 source and a CPT source in a single phenotype definition.

Each phenotype must define a `sources` dict. Keys are arbitrary source labels; each source specifies:
- `from_phenotype` — name of a `coded` phenotype (defined elsewhere in this YAML) whose preprocessed output feeds this source. Required.
- `case_codes` — code patterns that make a patient a case. Required. Works with any vocabulary.
- `case_exclude` — code patterns that set a patient to NA, whether or not they have a case code (overrides case and control).
- `control_exclude` — code patterns that mark a non-case as NA rather than 0 (ambiguous controls).

When multiple sources are defined, case/exclude lists are pooled across all sources before classification.

#### Code matching rules

How a pattern in `case_codes`, `case_exclude` or `control_exclude` (and in event-mask `codes`) matches the codes stored in your data. Codes are compared as text and are expected to be stored **with the ICD dot** (`K74.60`, not `K7460`).

| Pattern | Rule | Matches | Does not match |
|---|---|---|---|
| `N80` (bare) | The code itself **and** its subcodes | `N80`, `N80.0`, `N80.11` | `N800`, `N80X0`, `N81.0` |
| `C61` (bare, billable with no subcodes) | Same rule — the bare code itself matches | `C61` | `C610`, `C62` |
| `N80.*` | Subcodes only — **not** the bare code | `N80.0`, `N80.11` | `N80`, `N80X0` |
| `N80.0` (dotted, no `*`) | Exact match only | `N80.0` | `N80.00`, `N80.01`, `N80` |
| `K74.6` (dotted header code) | Exact match only — header codes are never billed, so this usually matches **no one** | — | `K74.60`, `K74.69` |
| `K74.6*` | Prefix match on `K74.6` | `K74.6`, `K74.60`, `K74.69` | `K74.5`, `K74.0` |
| `I70.0**` | Prefix match on `I70.0` (`**` same as `*`) | `I70.0`, `I70.00`, `I70.09` | `I70.10` |
| `I83.` (trailing dot) | Prefix match on `I83.` | `I83.0`, `I83.90` | `I83`, `I830` |
| `250` (ICD-9, bare) | The code itself and its subcodes | `250`, `250.00`, `250.13` | `2500`, `251.0` |
| `99213` with `prefix_matching: false` | Bare codes are exact | `99213` | `99213.1` |

Rules of thumb:

- **To include subcodes of a dotted code, add `*`.** Papers usually list header codes (`K74.6`, `K70.3`, `I85.0`) meaning "and its subcodes" — write `K74.6*`.
- **Write a bare code (`N80`) when you want the category, `N80.*` only when you specifically want to leave out a stored `N80`.** Many 3-character ICD-10 codes are billable with no subcodes (`C61` prostate cancer, `I10` hypertension), so `C61.*` finds no one.
- For non-hierarchical vocabularies like CPT or phecodes, set `prefix_matching: false` (at phenotype or per-source level) so bare codes are exact. Explicit wildcards still prefix-match.

**Patterns that match nothing are flagged in the manifest.** For each `advanced_coded` phenotype, any `case_codes` / `case_exclude` / `control_exclude` pattern that matches no code in the phenotype's input data is recorded in the `warnings` column of [`phenotype_manifest.tsv`](#phenotype_manifesttsv). When codes starting with the pattern exist in the data, up to five are listed as a hint — e.g. `case_codes pattern 'K74.6' matched no codes (data has K74.60, K74.69)`. No hint means the code is absent from the source entirely (not in this release, or withheld). The check runs on all input rows, before `min_occurrences`, and only against your own data — no external code lists. A phenotype that ends up with `0 cases` is flagged too. At the end of the run, the log prints one line saying how many advanced phenotypes have warnings; the full messages are also in each ADVANCED_CODED task's `.command.err`.

```yaml
T2Diab:
  data_type: advanced_coded
  missing_as_control: false
  sources:
    ICD10:
      from_phenotype: ICD10      # references the 'ICD10' coded phenotype above
      case_codes:
        - E11       # prefix match: E11.9, E11.65, etc.
        - O24.1     # exact match (has a dot)
      case_exclude:
        - E10       # type 1 diabetes → NA (excluded from cases and controls)
        - E12
        - E13
```

See [`examples/advanced_coded_example.yaml`](examples/advanced_coded_example.yaml) for more examples including `control_exclude` and multi-source configs.

### Coded Column Prefixes and Multiple Vocabularies

Basic `coded` phenotypes explode into one output column per code, named after the code itself. When a config has more than one coded source whose code formats overlap — most commonly **ICD-9 and PheCodes**, which are both numeric (`401.1` exists in both) — the exploded outputs collide:

- both write `phenotypes/401.1.diag.tsv` (one overwrites the other),
- both write `long_format/401.1.long.tsv`,
- GATHER sees two `401.1` columns.

#### Recommended: one run per coded vocabulary

**It is best to run the pipeline once per coded phenotype/vocabulary** (e.g. one run for ICD10, one for ICD9, one for PheCodes), each with its own `--output_dir`. Every output directory (`long_format/`, `phenotypes/`, the gathered table) then contains only one vocabulary, so nothing can collide, and each run is smaller and easier to rerun. Keep any `advanced_coded` phenotypes in the run(s) that contain the coded phenotypes they reference via `from_phenotype`.

#### Alternative: `column_prefix` in a single run

If you need several coded vocabularies in one run, give each coded phenotype a `column_prefix`:

```yaml
phenotypes:
  ICD9:
    data_type: coded
    table: icd9.tsv
    value_col: code
    column_prefix: "icd9_"     # → icd9_401.1
  PheCode:
    data_type: coded
    table: phecodes.tsv
    value_col: phecode
    column_prefix: "phe_"      # → phe_401.1
```

`column_prefix` can also be set once under `data_types.coded:` and overridden per phenotype. Note that `long_format/` can still collide in a single run (`401.1.long.tsv` from both sources), because long-format files are intentionally left unprefixed (see below) — another reason to prefer separate runs.

#### Where the prefix is applied (and where it is not)

The prefix is applied at the **very last step only** — when BASIC_CODED names its output column. Everything upstream and every consumer of the long-format files sees raw codes.

| Stage / file | Prefixed? | Details |
|---|---|---|
| `bin/resolve_config.py` | — | Accepts `column_prefix` (default `""`), validates it is a string, and rejects it on `advanced_coded` phenotypes. Passes it through in the coded JSON. |
| `bin/preprocess.py` → `long_format/{code}.long.tsv` | **No** | Filenames and the `phenotype`/`value` columns keep the raw code exactly as it appears in the input table. |
| `bin/advanced_coded.py` | **No** | Matches `case_codes` / `case_exclude` / `control_exclude` against the raw `value` column of the long files. Output column is the advanced phenotype's own name (e.g. `T2Diab`). |
| `bin/event_mask.py` | **No** | Derives codes from the raw long-file names. |
| `bin/basic_coded.py` | **Yes** | Output column = `{column_prefix}{code}`; output file = `{column_prefix}{code}.diag.tsv`. |
| `modules/basic_coded/main.nf` | **Yes** | Publishes `phenotypes/{column_prefix}{code}.diag.tsv`. |
| `bin/phenotype_manifest.py` → `phenotype_manifest.tsv` | Both | `column_name` is prefixed; the `code` column is the raw code with the prefix stripped. |
| `phenotype_table.wide.tsv` (GATHER) | **Yes** | Contains the prefixed basic columns (and unprefixed advanced phenotype columns). |

> **Writing your own `advanced_coded` phenotypes:** always use the **raw codes as they appear in your input tables** (e.g. `E11`, `401.1`, `250.2`) in `case_codes`, `case_exclude`, `control_exclude`, and event-mask `codes`. **Do not** write patterns with a `column_prefix` in mind (e.g. `phe_250.2`) — they will never match, because advanced_coded never sees prefixed codes. This is what keeps the shared/community advanced phenotype definitions portable: the same definition works whether or not a given biobank config uses prefixes.

#### Finding coded columns after the run

[`phenotype_manifest.tsv`](#phenotype_manifesttsv) lists every exploded coded column with its source phenotype and raw code, so you can select e.g. all `PheCode` columns from the gathered table regardless of prefix. Two rows with the same `column_name` indicate a collision.

### `first` / `last` stats

`first` and `last` report each sample's value at its **earliest** / **latest** `occurrence_date` — e.g. baseline and most recent BMI.

- **`date_col` is required.** Without one, PREPROCESS would fill every row with the same `reference_date` and "first" would just mean "first row in the file", so RESOLVE_CONFIG stops with an error instead. (Not required with `preprocessed_path`, whose `occurrence_date` column is used as-is.)
- **The date column must exist in the table.** A `date_col` that isn't a column in the table — a typo, or a `global:` `date_col` inherited by a table that doesn't have it — stops PREPROCESS with an error. (For other stats a missing date column silently falls back to `reference_date`, as before.)
- **Rows with a blank or unparseable date are skipped** by `first` / `last` — PREPROCESS leaves their `occurrence_date` blank instead of filling `reference_date` (which defaults to today and would otherwise make an undated row look like the most recent measurement), and logs how many rows had no date. Those rows still count toward the other stats (`mean`, `count`, …) and `min_occurrences`. A sample whose rows are all undated gets NA.
- **Same-date ties** keep the order the rows appear in the source table.
- `first` / `last` work with `min_occurrences` and `qcut_bins` like any other stat.

### `min_occurrences` behavior

`min_occurrences` is applied during phenotyping, **not** during PREPROCESS. All rows are preserved in long-format output so downstream steps have full longitudinal data.

- Quantitative: sample with < N measurements → all stat columns are `NA`
- Categorical: sample with < N occurrences → treated as missing (NA or 0 per `missing_as_control`)
- Diagnostic (basic): sample with < N rows for a given code → that code column is NA
- Diagnostic (advanced): a code qualifies for pattern matching only if it appears ≥ N times **for that specific code**. For example, with `min_occurrences: 2`, a patient with one `E11.9` and one `E11.65` visit is **not** a case — each individual code falls short, even though they have two E11-family encounters total. If you want any qualifying encounter to count, set `min_occurrences: 1` (the default). Two further consequences:
  - **It counts raw rows, not distinct dates.** Two rows of the same code on the same day count as 2.
  - **It applies to exclusion codes too.** The threshold is applied to every code before any pattern matching, so under `min_occurrences: 2` a patient with a single `case_exclude` or `control_exclude` code is **not** excluded. Raising `min_occurrences` can therefore *increase* a phenotype's case count — e.g. an HCC phenotype with a non-HCC-cancer exclusion had 191 cases at `min_occurrences: 1` and 220 at `2`, because once-recorded exclusion codes stopped excluding.

---

## Event Masking (`--event_mask_config`)

Event masking excludes observations within date windows triggered by clinical events
(e.g., pregnancy) from phenotype calculations. This runs **after preprocessing** and
**before** phenotype classification.

### How it works

1. Define named events in an `event_mask.yaml` file (see `examples/event_mask_example.yaml`).
   Each event has ICD codes, a source coded phenotype, and day offsets:

   ```yaml
   events:
     third_trimester_pregnancy:
       source: ICD10          # references an existing coded phenotype
       codes: [O80, O82, O83, Z37]   # bare codes = prefix match
       offset_before: -270   # days before event date
       offset_after: 90      # days after event date
   ```

2. Opt a phenotype into masking by adding `event_mask` to its config:

   ```yaml
   phenotypes:
     HbA1c:
       data_type: quantitative
       date_col: ENCOUNTER_DATE
       event_mask: [third_trimester_pregnancy]
       ...
   ```

3. Run the pipeline with `--event_mask_config path/to/event_mask.yaml`.

### What gets masked

For each opted-in phenotype, any observation row whose `occurrence_date` falls within a
mask window for that sample is **excluded** before downstream calculation. This means:

- **`basic_coded`**: Code occurrences within windows don't count toward `min_occurrences`.
  A sample with all code occurrences masked may become control/NA.
- **`basic_quantitative`**: Measurements within windows are excluded from summary stats and binning.
- **`basic_categorical`**: Category values within windows are excluded from frequency/encoding.

Overlapping windows (from multiple events or multiple event occurrences) are merged per sample.

### Requirements and limitations

> **Phenotypes without a real `date_col` cannot be event-masked.** When `date_col` is
> not set, `occurrence_date` is filled from `reference_date` (a synthetic date shared
> by all observations). Masking on synthetic dates is meaningless — `apply_mask` skips
> masking for these phenotypes and logs a WARNING.

### Design note (Solution 1 — raw-data window computation)

Event windows are **always** computed from raw (unmasked) long.tsvs. If a phenotype is
both the source for an event and opts into masking by that same event, the self-masking
behavior is valid and intentional: event code occurrences that triggered a window are
also filtered from that phenotype's data during the window period. This is logged as
a warning but is not an error.

---

## Risk Classification

Risk classification runs after `GATHER` and classifies each sample's aggregated quantitative value into a binary (high-risk / not) or multi-level (low / moderate / high) category. Rules are defined in the biobank YAML — no separate params file needed.

### Enabling

```bash
nextflow run main.nf \
  --biobank_config biobanks/PMBB/pmbb_config.yaml \
  --run_gather true \
  --run_risk_classification true \
  --output_dir results/
```

`--run_gather` is required because risk classification reads the gathered wide TSV.

### Config

Rules are defined under `data_types.quantitative.risk_classification` (applies to all quantitative phenotypes) and/or per phenotype. Per-phenotype config completely replaces the data_type-level default. Opt out with `risk_classification: false`.

```yaml
data_types:
  quantitative:
    risk_classification:
      method: percentile   # percentile | threshold | quantile_bin
      stat: mean           # aggregated stat column to classify on (default: mean)
      high: 90             # upper tail cutoff
      low: 10              # lower tail cutoff (omit for one-sided)

phenotypes:
  BP_systolic:
    risk_classification:   # override: use fixed clinical threshold
      method: threshold
      high: 130
      low: 90

  HDLC:
    risk_classification: false   # opt out
```

### Methods

| Method | Cutoff scale | Description |
|---|---|---|
| `percentile` | 0–100 | Cutoffs computed from data via `np.nanpercentile`. `high: 90` → value > p90 = case. |
| `threshold` | raw value | Fixed absolute cutoffs, no data-driven computation. `high: 130` → value > 130 = case. |
| `quantile_bin` | bin counts | Equal-frequency bins via `pd.qcut(q=N)`. `high: k` → top k bins = case; `low: k` → bottom k bins = case. |

### Binary vs multi-level output

**Binary** — use `high` and/or `low`. Both tails collapse to 1; all other values → 0.

```yaml
risk_classification:
  method: percentile
  high: 90    # value > p90 → 1
  low: 10     # value < p10 → 1 (both tails flagged)
```

**Multi-level (percentile / threshold)** — use `cutoffs` list instead of `high`/`low`. N cutoffs → N+1 integer levels (0 = lowest). Optional `labels` replaces integers with strings.

```yaml
risk_classification:
  method: threshold
  cutoffs: [200, 240]
  labels: [normal, borderline, high]
```

**Multi-level (quantile_bin)** — omit `high`/`low`. All N bin indices become output levels.

```yaml
risk_classification:
  method: quantile_bin
  n: 5
  labels: [Q1, Q2, Q3, Q4, Q5]
```

### Output naming

The risk output column is always `{output_name}_risk` where `output_name` is the same value used for all other output columns (phenotype key, or the explicit `output_name` field if set). The source column read from the wide TSV is `{output_name}_{stat}`.

### Outputs

Written to `{output_dir}/risk_classification/`:

- **`risk_matrix.tsv`** — `sample_ID` + one column per rule. Binary columns are float (0.0 / 1.0 / NA); multi-level columns are integer or string.
- **`risk_report.txt`** — computed cutoffs and per-class counts for each rule, for reproducibility.

> **Risk columns are kept separate.** RISK_CLASSIFICATION reads `phenotype_table.wide.tsv` but never writes back to it, so `*_risk` columns are **not** in the gathered table and **not** listed in `phenotype_manifest.tsv`. To combine them, join `risk_matrix.tsv` to `phenotype_table.wide.tsv` on `sample_ID`.

### Missing data

Samples with `NA` in the source column are excluded from cutoff computation and receive `NA` in the risk output column.

---

## EDA Report (`--run_eda`)

An optional QC pass over every phenotype output column, designed to stay cheap
at thousands of phenotypes × tens of thousands of samples: **metrics first,
profiles second.** Every column gets cheap metrics and flags; the expensive
[fg-data-profiling](https://pypi.org/project/fg-data-profiling/) (successor of
ydata-profiling) HTML profiles run only on a bounded subset. It reads the
per-phenotype results directly, so it does **not** need `--run_gather`.

```bash
nextflow run main.nf --biobank_config my_config.yaml --sample_list samples.txt \
    --output_dir results --run_eda
```

In `eda/eda_summary.tsv`, `data_type` is `advanced_coded` for advanced coded phenotype columns (taken from the phenotype manifest) and `coded` for exploded per-code columns; both are treated as binary.

The whole `eda/` folder (report, profiles, tables) is also zipped to
`{output_dir}/eda_report.zip` for easy download. `eda/run_config.html` (linked
from the top of the report) records the run: start time, profile, pipeline
revision, the parameters logged at startup, and the effective global,
per-data_type and `eda:` settings, each marked with where it came from
(`default`, `global`, `data_types` or `yaml`). Phenotype-level overrides are not
listed; see `configs/<data_type>/<phenotype>.json`. `eda/biobank_config.yaml` is a
verbatim copy of the YAML the run used (so the zip still documents the run if
the original is edited later), and `eda/biobank_config.html` shows it in the
browser with line numbers and jump links; both pages are linked from the report. It has a fixed name and is
**overwritten on every run**, like the rest of `{output_dir}`: copy it
elsewhere before re-running if you want to keep it.

Pass `--sample_list`: it is the denominator for missingness (quantitative
results only list measured samples) and, for coded phenotypes, the roster that
makes controls exist at all.

### What gets flagged

| Flag | Rule (defaults) |
|---|---|
| `separation_complete` | All cases (or controls) of a binary column fall in one covariate level — or in a subset of up to `max_subset_levels` (3) levels — with ≥ `separation_min_cases` (10) cases in total and ≥ `min_expected` (5) expected outside it |
| `separation_quasi` | ≥ `separation_threshold` (0.95) of them do, while that level / subset holds < 95% of samples |
| `low_cases` | Binary, fewer than `min_cases` (100) cases |
| `no_controls` | Binary with cases but no controls (no `--sample_list` / `missing_as_control: false`) — suppresses the two missingness flags |
| `high_missing` | Missing fraction > `max_missing` (0.5) |
| `missing_by_covariate` | Missing-fraction gap between covariate levels > `missing_gap` (0.8), e.g. a lab measured almost only in one sex |
| `empty` / `zero_variance` / `single_level` | No data / a constant column |
| `pairwise_duplicate` / `pairwise_nested` | Opt-in (`pairwise: true`): binary case sets with Jaccard ≥ 0.95, or ≥ 99% contained in another column |

Covariate levels with fewer than `min_level_n` (50) samples are ignored.
Numeric covariates with more than `max_covariate_levels` (10) values are
binned into `age_bins` (5) quantiles.

### Covariates

Set in the top-level `eda:` block. Preferred: name phenotypes already defined
in the same YAML (`covariates: [sex, batch, age]`) — their pipeline outputs are
reused, so the coding comes from their `dictionary`. Quantitative covariates use
their first stat column. Alternative (or in addition): `covariates_file`, an
external table keyed by `global.sample_id_col`, with optional
`covariates_file_cols`.

### Keeping profiles bounded

Profiles = columns opted in via the `eda` cascade **plus** the
`profile_flagged_top_k` (50) most severe flagged columns, capped at
`max_profile_columns` (200), split into one report per phenotype (chunks of
`profile_group_size`, 50) that run in parallel. Each profile is downsampled to
`profile_max_samples` (20000), keeping cases of rare binary columns, in minimal
mode with correlations off (`profile_correlations: true` turns on `auto`). If
the profiler isn't installed, the `EDA_PROFILE` tasks exit with status 3 and
are ignored, so everything else still runs. Because they aren't cached as
successful, `-resume` re-runs them once the profiler is installed.

### Sex- / cohort-specific phenotypes (separation)

Detection is **data-driven only**, against every covariate you list (sex,
cohort, site, binned age, …). For each binary column and covariate, the EDA
finds the smallest set of levels — one level, or up to `max_subset_levels` for
multi-level covariates like cohort — holding ≥ `separation_threshold` of the
cases. It is reported per column in `eda/eda_summary.tsv`:

- **Rare columns are never flagged.** It needs ≥ `separation_min_cases` (10)
  cases, and ≥ `min_expected` (5) cases *expected* outside the level / subset
  if cases were spread like samples. 1–2 cases that happen to be female fail
  both, and so does any column with ≤ 5 cases at the defaults.
- **Combinations of covariates** (`separation_interactions`): each listed
  combination becomes one more covariate, e.g. `COHORT&SEX` with levels like
  `1&0` (cohort 1, female); `auto` adds every pair. Used for separation only.
  An interaction result that a single covariate already explains (cases in
  cohort 1, both sexes → `1&0|1&1`) is dropped.


| Column | Meaning |
|---|---|
| `sep_type` | `complete` (all cases inside) or `quasi` (≥ `separation_threshold`) |
| `sep_covariate` / `sep_level` | Where the cases are, e.g. `COHORT` / `1\|2` or `COHORT&SEX` / `1&0` (level codes from `eda_covariates.tsv`; a subset is joined with `\|`, an interaction with `&`) |
| `sep_outcome` | `cases` (or `controls`, if the controls are the confined group) |
| `sep_share` | Fraction of cases (controls) inside that level / subset, e.g. `0.993` |
| `sep_detail` | Every covariate that showed separation, not just the top one |

`eda/restriction_suggestions.tsv` lists, per (column, covariate, level), the
`data` separations and any `restrict_to` you declared (`yaml`), with `share`
and readable `level_label`s. These are suggestions only — the EDA never changes
phenotype data. To act on them, use the [filter step](#filtering---filter_config)
(`restrict_separated`, or explicit `restrictions`).

Known sex-specific codes (ICD-9-CM, ICD-10-CM, phecode 1.2, phecodeX) are not
built into the pipeline; ready-made filter `restrictions:` blocks are in
`community/sex_specific/`.

### Config

```yaml
eda:
  covariates: [sex, batch, age]     # phenotypes in this YAML
  # covariates_file: /path/covariates.tsv
  # covariates_file_cols: [ANCESTRY, SITE]
  min_cases: 100
  max_missing: 0.5
  separation_threshold: 0.95        # share of cases inside a level / subset → flagged
  separation_min_cases: 10          # fewer cases → never tested for separation
  max_subset_levels: 3              # largest level subset tested (multi-level covariates)
  separation_interactions: [[cohort, sex]]   # or: auto (every pair); default []
  min_level_n: 50
  min_expected: 5
  missing_gap: 0.8
  profile_flagged_top_k: 50
  max_profile_columns: 200
  profile_group_size: 50
  profile_max_samples: 20000
  pairwise: false                   # pairwise_jaccard / pairwise_containment / pairwise_top_k / pairwise_scope

global:
  eda: false                        # optional: opt everything out, then opt in below
data_types:
  quantitative: {eda: true}
phenotypes:
  ICD10:  {eda: ["E11*", I10]}      # profile only these codes of an exploded coded phenotype
  PCOS:   {restrict_to: {sex: Female}}
```

Unknown `eda:` keys stop the run at RESOLVE_CONFIG. The full list with defaults is
`EDA_DEFAULTS` in `bin/resolve_config.py`; every key is annotated in
`examples/eda_example.yaml`, and `community/eda.yaml` has recommended production settings.

---

## Filtering (`--filter_config`)

The pipeline's outputs are raw: nothing in `phenotypes/` or
`phenotype_table.*.tsv` is ever filtered. Filtering is a separate, opt-in
post-processing step that writes **new** files under `{output_dir}/filtered/`.

1. Run with `--run_eda`. It writes `eda/filter_suggested.yaml`, which contains:
   - drops for `empty`, `zero_variance`, `single_level`, `low_cases` and
     `high_missing` columns
   - active restrictions for every `restrict_to`
   - a commented-out `restrict_separated` rule at the EDA's
     `separation_threshold`, listing the columns it would restrict — separation
     can be a real association, so this needs a judgement call
   - commented-out "review" entries for pairwise duplicates and
     missing-by-covariate columns
2. Copy it, edit it, and rerun with the filter step. `-resume` reuses everything
   else:

```bash
cp results/eda/filter_suggested.yaml my_filter.yaml   # edit
nextflow run main.nf <same args> --run_gather --run_eda \
    --filter_config my_filter.yaml -resume
```

Or run it standalone on a finished output directory:

```bash
filter_phenotypes.py --config my_filter.yaml \
    --input results/phenotype_table.wide.tsv \
    --covariates results/eda/eda_covariates.tsv \
    --summary results/eda/eda_summary.tsv \
    --manifest results/phenotype_manifest.tsv \
    --output_dir results/filtered
```

| Key | Meaning |
|---|---|
| `drop_columns` | Column names or globs to remove |
| `keep_columns` | Names or globs that are always kept — overrides `drop_columns` and `drop_flags` |
| `drop_flags` | Drop every column carrying any of these EDA flags, e.g. `[low_cases, empty]`. A rule-based alternative to explicit lists |
| `restrictions` | `{column or glob: {covariate: C, level: L}}`. Samples whose `C` isn't `L` are set to NA, cases included; those cases are counted in the log as likely miscodes. `L` is the level code shown in `restriction_suggestions.tsv` (from `eda_covariates.tsv`), or a list of codes (`level: ["1", "2"]`) |
| `restrict_separated` | `{min_share: S, covariates: [C, ...]}`. Rule-based: restricts every column whose cases the EDA found confined to a level (or subset of levels) of `C` with `sep_share` ≥ `S` (from `eda_summary.tsv`). `C` may be an interaction, e.g. `COHORT&SEX`. Omit `covariates` for all. A column's explicit `restrictions` entry wins for that covariate. Only separations at or above the EDA's `separation_threshold` are recorded, so `S` below it adds nothing |
| `restrict_unknown` | `na` (default) or `keep`: what happens to samples with no value for `C` |

Annotated example: `examples/filter_example.yaml`. Recommended rule-based starting point: `community/filter.yaml`.
Sex-specific code lists to paste into `restrictions:`: `community/sex_specific/`.

`filter_log.tsv` lists every column that was dropped, kept by override, or
restricted, with the reason and how many values and cases were set to NA. Kept
cells are byte-identical to the input. Requires `--run_gather` with a wide table.

---

## Standalone Script Usage

All Python scripts can be run directly without Nextflow.

### resolve_config.py

```bash
python bin/resolve_config.py \
  --config biobanks/PMBB/pmbb_config.yaml \
  --output_dir /tmp/configs

# outputs: /tmp/configs/quantitative/BMI.json, /tmp/configs/coded/ICD10.json, etc.
```

### preprocess.py

```bash
python bin/preprocess.py \
  --config /tmp/configs/quantitative/BMI.json \
  --output_dir /tmp/long_format

# outputs: /tmp/long_format/BMI.long.tsv
```

### basic_coded.py

```bash
python bin/basic_coded.py \
  --input /tmp/long_format/N80.0.long.tsv \
  --config /tmp/configs/coded/ICD10.json \
  --output /tmp/N80.0.diag.tsv

# Or let the script name the file from column_prefix + code:
python bin/basic_coded.py \
  --input /tmp/long_format/401.1.long.tsv \
  --config /tmp/configs/coded/PheCode.json \
  --output_dir /tmp/
# outputs: /tmp/phe_401.1.diag.tsv
```

### phenotype_manifest.py

Lists the output columns of one per-phenotype result (any data type). The pipeline runs it after every phenotyping step and concatenates the results into `phenotype_manifest.tsv`.

```bash
python bin/phenotype_manifest.py \
  --result /tmp/phe_401.1.diag.tsv \
  --config /tmp/configs/coded/PheCode.json \
  --output /tmp/phe_401.1.manifest.tsv
```

### advanced_coded.py

Accepts one or more long.tsvs from basic coded preprocessing; all are pooled before case/control classification.

```bash
python bin/advanced_coded.py \
  --input /tmp/long_format/E11.9.long.tsv /tmp/long_format/E11.65.long.tsv \
  --config /tmp/configs/advanced_coded/T2Diab.json \
  --output /tmp/T2Diab.diag.tsv
```

### basic_quantitative.py

```bash
python bin/basic_quantitative.py \
  --input /tmp/long_format/BMI.long.tsv \
  --config /tmp/configs/quantitative/BMI.json \
  --output /tmp/BMI.quant.tsv \
  --report /tmp/BMI_report.txt
```

### basic_categorical.py

```bash
python bin/basic_categorical.py \
  --input /tmp/long_format/SEX.long.tsv \
  --config /tmp/configs/categorical/SEX.json \
  --output /tmp/SEX.cat.tsv \
  --codes  /tmp/SEX.codes.tsv   # optional; written only when auto encoding assigns codes
```

### gather.py

```bash
python bin/gather.py \
  --inputs /tmp/*.quant.tsv /tmp/*.cat.tsv /tmp/*.diag.tsv \
  --output_prefix /tmp/phenotype_table \
  --format wide   # wide | long | both
```

### risk_classification.py

```bash
python bin/risk_classification.py \
  --wide_tsv /tmp/phenotype_table.wide.tsv \
  --risk_config /tmp/configs/risk_classification_config.json \
  --output /tmp/risk_matrix.tsv \
  --report /tmp/risk_report.txt
```

`risk_classification_config.json` is emitted by `resolve_config.py` when at least one quantitative phenotype has a `risk_classification` rule defined in the YAML.

---

## Running Tests

```bash
pytest tests/ -v
```

| Test file | Coverage |
|-----------|----------|
| `tests/test_resolve_config.py` | `merge_params` inheritance, `validate` errors (incl. `column_prefix`), `resolve` integration, path absolutization |
| `tests/test_preprocess.py` | `apply_filter` (glob wildcards, AND/OR logic), `get_matching_codes`, `run_preprocess` (quant/cat/diag/preprocessed_path) |
| `tests/test_basic_coded.py` | `run_basic_coded` (present→1, absent→NA, missing_as_control, --all_samples roster, min_occurrences, column_prefix, long.tsv untouched) |
| `tests/test_phenotype_manifest.py` | `run_phenotype_manifest` (row per column for each data_type, raw `code` with prefix stripped, chunked configs use YAML key) |
| `tests/test_advanced_coded.py` | `run_advanced_coded` (case→1, case_exclude→NA (cases and controls), control_exclude→NA, prefix matching incl. bare codes stored as-is, min_occurrences per-code semantics, unmatched-pattern / 0-case warnings) |
| `tests/test_basic_quantitative.py` | `run_basic_quantitative` (mean/std, squared, first/last ordering + undated rows, unknown-stat error, min_occurrences→NA, qcut_bins, output_name) |
| `tests/test_basic_categorical.py` | `run_basic_categorical` (auto, binary, dictionary, one-hot, encoding precedence, missing_as_control, min_occurrences) |
| `tests/test_gather.py` | `run_gather` (wide join, outer join NA fill, long format, both formats) |
| `tests/test_icd_utils.py` | ICD code matching utilities, including `normalize_patterns` |
| `tests/test_event_mask.py` | `run_event_mask` (window computation, prefix matching, merging, source filtering) |
| `tests/test_apply_mask.py` | `run_apply_mask` (filtering, boundaries, passthrough for no date_col, event type filtering) |
| `tests/test_risk_classification.py` | `apply_rule` (percentile/threshold/quantile_bin binary and multi-level, NA propagation, all-NA guard), `run_risk_classification` integration |
| `tests/test_split_samples.py` | Sample chunking utilities |
| `tests/test_merge_wide.py` | Wide-format merge utilities |
| `tests/test_run_config.py` | `build_run_config` (setting sources follow the merge order, elapsed time, only data types in use), HTML/JSON output |

---

## Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `biobank_config` | — | **Required.** Path to biobank YAML config |
| `output_dir` | `results` | Output directory |
| `run_gather` | `false` | Join all per-phenotype TSVs into a final table |
| `gather_format` | `wide` | `wide`, `long`, or `both` |
| `run_risk_classification` | `false` | Enable risk classification (requires --run_gather; config defined in pmbb_config.yaml) |
| `coded_chunk_size` | `500` | Split coded phenotypes into chunks of N codes, each running as a separate parallel job. Set to `0` to disable (single job). Recommended for biobanks with thousands of unique ICD codes. |
| `chunk_size` | `0` | Split the long-format input into sample partitions of this size for parallel processing. `0` = disabled. See module comments for per-module caveats. |
| `event_mask_config` | `null` | Path to an `event_mask.yaml` config. When set, EVENT_MASK computes windows from event ICD codes and APPLY_MASK filters opted-in phenotypes before downstream calculation. |
| `run_eda` | `false` | EDA / QC report; settings in the YAML's `eda:` block. See [EDA Report](#eda-report---run_eda). |
| `eda_batch_size` | `200` | Result files per `EDA_STATS` task. |
| `filter_config` | `null` | Filter YAML for the post-processing filter step (needs `--run_gather` and `--run_eda`). |
| `sample_list` | `null` | Path to a file with one sample ID per line. When provided, only those samples are processed; all others are dropped at the PREPROCESS step. It is also the roster for coded/advanced_coded: samples in the list without a case code get 0 when `missing_as_control: true`. Without it, coded columns have only 1/NA. |

These are the only pipeline parameters. Phenotype settings (outliers, stats,
encodings, …) go in the biobank YAML. Older parameters such as
`--outlier_method`, `--quant_stats` or `--run_quantitative` had no effect and
have been removed; passing one now logs a warning.

---

## Project Structure

```
pmbb-nf-toolkit-phenotyping/
├── main.nf                            # Scatter-gather workflow
├── nextflow.config                    # Parameters and profiles
├── bin/
│   ├── resolve_config.py              # Biobank YAML → per-phenotype JSONs
│   ├── preprocess.py                  # Per-phenotype preprocessing
│   ├── basic_coded.py            # Binary presence/absence per code
│   ├── advanced_coded.py         # Named case/control phenotype from ICD patterns
│   ├── basic_quantitative.py          # Per-sample quantitative stats
│   ├── basic_categorical.py           # Categorical encoding
│   ├── gather.py                      # Outer-join per-phenotype TSVs
│   ├── phenotype_manifest.py          # Per-phenotype rows for phenotype_manifest.tsv
│   ├── event_mask.py                  # Per-source window computation from ICD long.tsvs
│   ├── apply_mask.py                  # Filter a long.tsv against mask_intervals.tsv
│   ├── risk_classification.py         # Apply risk rules to gathered wide TSV
│   ├── split_samples.py               # Legacy: chunk long-format by sample_ID
│   ├── merge_wide.py                  # Legacy: merge wide-format chunks
│   └── icd_utils.py                   # ICD/data_type utilities (prefix matching, normalize_patterns)
├── modules/
│   ├── resolve_config/main.nf
│   ├── preprocess/main.nf
│   ├── basic_coded/main.nf
│   ├── advanced_coded/main.nf
│   ├── basic_quantitative/main.nf
│   ├── basic_categorical/main.nf
│   ├── gather/main.nf
│   ├── event_mask/main.nf
│   ├── apply_mask/main.nf
│   ├── risk_classification/main.nf
│   ├── split_samples/main.nf          # Legacy chunked mode
│   └── merge_wide/main.nf             # Legacy chunked mode
├── examples/
│   ├── basic_config_example.yaml      # Annotated example: quant/cat/coded phenotypes
│   ├── advanced_coded_example.yaml  # Advanced coded with sources dict
│   ├── event_mask_example.yaml        # Reference event_mask config with pregnancy events
│   ├── eda_example.yaml               # Every eda: key, eda cascade, restrict_to
│   └── filter_example.yaml            # Every --filter_config key
├── community/
│   ├── advanced_coded.yaml            # Production advanced_coded definitions
│   ├── event_mask.yaml                # Production obstetric event mask
│   ├── risk_classification.yaml       # Production risk classification rules
│   ├── eda.yaml                       # Recommended eda: block for production-scale runs
│   ├── filter.yaml                    # Recommended --filter_config (drop_flags, restrict_separated)
│   └── sex_specific/                  # Sex-specific ICD-9/ICD-10/phecode 1.2/phecodeX restrictions blocks
├── biobanks/
│   └── PMBB/
│       ├── pmbb_config.yaml           # PMBB phenotype/table config (update paths)
│       ├── pmbb_advanced.yaml         # Advanced coded phenotype definitions
│       └── pmbb_config.conf           # Nextflow parameter overrides for PMBB runs
├── conf/
│   ├── base.config
│   ├── test.config                    # Uses test_data/PMBB/
│   └── slurm.config
├── test_data/
│   ├── PMBB/
│   │   ├── pmbb_test_config.yaml
│   │   ├── pmbb_test_filter.yaml      # Synthetic --filter_config for the test cohort
│   │   ├── vitals_bmi.tsv
│   │   ├── labs_lipids.tsv
│   │   ├── covariates.tsv
│   │   ├── icd.tsv
│   │   └── filter/                    # Frozen EDA outputs + roster for standalone filter runs
│   └── ...                            # Legacy generic fixtures
└── tests/
    ├── conftest.py
    ├── test_resolve_config.py
    ├── test_preprocess.py
    ├── test_basic_coded.py
    ├── test_advanced_coded.py
    ├── test_basic_quantitative.py
    ├── test_basic_categorical.py
    ├── test_gather.py
    └── ...
```
