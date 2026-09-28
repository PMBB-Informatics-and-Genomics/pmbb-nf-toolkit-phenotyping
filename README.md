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

Example:

```
column_name     phenotype  data_type       code
T2Diab          T2Diab     advanced_coded
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
| `case_codes` | coded | Code patterns that make a patient a **case** (1). Works with any vocabulary — ICD, CPT, SNOMED, LOINC, phecodes, etc. Bare codes are prefix-matched by default (`E11` matches `E11.9`, `E11.65`); disable with `prefix_matching: false` for exact-match systems like CPT. Wildcards (`E11.**`) also work. |
| `case_exclude` | coded | Code patterns that **override a case to 0** — patient has a qualifying code but also one of these |
| `control_exclude` | coded | Code patterns that **exclude a non-case from being a control** (→ NA) |
| `prefix_matching` | coded, advanced_coded | `true` (default): bare codes prefix-match all sub-codes (e.g. `E11` → `E11.9`, `E11.65`). `false`: bare codes require an exact match — use for non-hierarchical systems like CPT or phecodes. Settable at phenotype level or per source in `advanced_coded` (source-level wins). |
| `sources` | advanced_coded | Dict of named sources. Keys are arbitrary labels (e.g. `ICD10`, `CPT`); each entry specifies which `coded` phenotype to read and what code patterns to apply. Sources can mix vocabularies. Patterns from all sources are pooled before case/control assignment. |
| `from_phenotype` | advanced_coded (per source) | Name of a `coded` phenotype defined in this YAML whose preprocessed output feeds this source. Required. |
| `case_codes` | advanced_coded (per source) | Code patterns that make a patient a **case**. Required. Bare codes prefix-match by default; set `prefix_matching: false` on the source for exact matching (e.g. CPT). |
| `case_exclude` | advanced_coded (per source) | Code patterns that **override a case to 0**. |
| `control_exclude` | advanced_coded (per source) | Code patterns that **exclude a non-case from being a control** (→ NA). |
| `preprocessed_path` | all | Skip PREPROCESS; pass this long-format TSV directly downstream |
| `output_name` | all | Override the phenotype name in output columns |
| `min_occurrences` | all | Minimum observations per sample; below threshold → NA at phenotyping step |
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
| `missing_as_control` | categorical, coded | `true` = absent sample → 0; `false` = absent → NA |
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
- `case_exclude` — code patterns that override a case assignment to 0 (patient qualifies but also has an exclusion code).
- `control_exclude` — code patterns that mark a non-case as NA rather than 0 (ambiguous controls).

When multiple sources are defined, case/exclude lists are pooled across all sources before classification.

**Prefix matching:** By default, bare codes with no dot automatically prefix-match all subcodes — `E11` matches `E11.9`, `E11.65`, `E11.319`, etc. Codes with a dot (e.g. `E10.5`) are exact matches. Explicit wildcards (`E11.**`) also work. For non-hierarchical vocabularies like CPT or phecodes, set `prefix_matching: false` (at phenotype or per-source level) to require exact matches.

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
        - E10       # type 1 diabetes overrides case assignment
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
- Diagnostic (advanced): a code qualifies for pattern matching only if it appears ≥ N times **for that specific code**. For example, with `min_occurrences: 2`, a patient with one `E11.9` and one `E11.65` visit is **not** a case — each individual code falls short, even though they have two E11-family encounters total. If you want any qualifying encounter to count, set `min_occurrences: 1` (the default).

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
| `tests/test_basic_coded.py` | `run_basic_coded` (present→1, absent→NA, missing_as_control, min_occurrences, column_prefix, long.tsv untouched) |
| `tests/test_phenotype_manifest.py` | `run_phenotype_manifest` (row per column for each data_type, raw `code` with prefix stripped, chunked configs use YAML key) |
| `tests/test_advanced_coded.py` | `run_advanced_coded` (case→1, case_exclude→0, control_exclude→NA, prefix matching, min_occurrences per-code semantics) |
| `tests/test_basic_quantitative.py` | `run_basic_quantitative` (mean/std, squared, first/last ordering + undated rows, unknown-stat error, min_occurrences→NA, qcut_bins, output_name) |
| `tests/test_basic_categorical.py` | `run_basic_categorical` (auto, binary, dictionary, one-hot, encoding precedence, missing_as_control, min_occurrences) |
| `tests/test_gather.py` | `run_gather` (wide join, outer join NA fill, long format, both formats) |
| `tests/test_icd_utils.py` | ICD code matching utilities, including `normalize_patterns` |
| `tests/test_event_mask.py` | `run_event_mask` (window computation, prefix matching, merging, source filtering) |
| `tests/test_apply_mask.py` | `run_apply_mask` (filtering, boundaries, passthrough for no date_col, event type filtering) |
| `tests/test_risk_classification.py` | `apply_rule` (percentile/threshold/quantile_bin binary and multi-level, NA propagation, all-NA guard), `run_risk_classification` integration |
| `tests/test_split_samples.py` | Sample chunking utilities |
| `tests/test_merge_wide.py` | Wide-format merge utilities |

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
| `run_quantitative` | `false` | Run the `BASIC_QUANTITATIVE` step. |
| `run_categorical` | `false` | Run the `BASIC_CATEGORICAL` step. |
| `event_mask_config` | `null` | Path to an `event_mask.yaml` config. When set, EVENT_MASK computes windows from event ICD codes and APPLY_MASK filters opted-in phenotypes before downstream calculation. |
| `sample_list` | `null` | Path to a file with one sample ID per line. When provided, only those samples are processed; all others are dropped at the PREPROCESS step. |

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
│   └── event_mask_example.yaml        # Reference event_mask config with pregnancy events
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
│   │   ├── vitals_bmi.tsv
│   │   ├── labs_lipids.tsv
│   │   ├── covariates.tsv
│   │   └── icd.tsv
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
