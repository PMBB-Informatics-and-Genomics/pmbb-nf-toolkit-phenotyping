Each ID (P001-P010) has 5 different ICD codes
The codes include a mix of ICD9 and ICD10
Each ID has codes representing different conditions, such as:

P001: Endometriosis, Autism, Gout
P002: Autism, Psoriasis, VTE, Urolithiasis
P003: Rheumatoid arthritis, Valvular heart disease, Sarcoidosis
P004: VTE, Aortic conditions, Gout variants
P005: Hypertension, Lupus, Varicose veins
P006: Hernia, Pulmonary embolism, Menstrual disorders
P007: Eating disorders, Cervical neoplasms, Amyloidosis
P008: Bladder cancer, Intestinal neoplasms, Breast cancer
P009: Neuropathy, Autism, Rheumatoid arthritis
P010: Valve disorders, Kidney stones, Psoriasis
P011: Has N80.0 which should be a positive match (case_include)
P012: Has N80.3 which should be excluded (case_exclude)

P013 who has:

N80.0 (which is in icd10_case_include)
N80.3 (which is in icd10_case_exclude)

P014: Has both:

Case include (N80.0)
Case exclude (N80.1)
ICD9 include (617)
ICD9 exclude (617.1)
Expected: 2 (exclusion overrides inclusion)

P015: Has only exclusion codes:

ICD9 exclude (617.1)
ICD10 case exclude (N80.2, N80.3)
Expected: 2

P016: Has a full range of ICD10 codes:

Case include (N80.0)
Multiple case excludes (N80.1-N80.4)
Expected: 2 (exclusion overrides inclusion)

P017: Has only ICD9 codes:

Case include (617)
Multiple excludes (617.1-617.4)
Expected: 2 (exclusion overrides inclusion)

P018: Mix of ICD9 and ICD10:

Both case includes (N80.0, 617)
Both case excludes (N80.1, 617.1)
Expected: 2 (exclusion overrides inclusion)

P019: Only exclusion codes:

ICD9 excludes (617.1, 617.2, 617.3)
ICD10 excludes (N80.3, N80.4)
Expected: 2

P020: Only high-numbered exclusion codes:

ICD10 excludes (N80.5, N80.6, N80.8)
ICD9 excludes (617.4, 617.5)
Expected: 2

PMBB test cohort — EDA / filter fixtures
========================================

test_data/PMBB/icd.tsv: P004 (male) carries a synthetic N80.0 (endometriosis)
row — a deliberate "miscode" so the filter's explicit SEX restriction has a
case to act on.

test_data/PMBB/covariates.tsv: synthetic COHORT column (A/B/C → codes 1/2/3 via
the COHORT dictionary in pmbb_test_config.yaml), an EDA covariate alongside SEX.
Coded phenotypes use missing_as_control: true so they have controls and are
tested for separation. Expected data-driven separation (eda_summary.tsv):
  E11.9         COHORT = 1     (A)    complete — single level
  N80.0         COHORT = 1|2   (A|B)  complete — multi-level subset
  ANCESTRY_EUR  COHORT = 3     (C)    complete — coincidence of a tiny cohort
  ANCESTRY_AFR  SEX&COHORT = 0&1|0&2 (female&A|female&B) — interaction only:
                neither SEX nor COHORT alone passes min_expected
(eda: separation_interactions: auto adds SEX&COHORT; separation_min_cases is
lowered to 2 for this tiny cohort.)

test_data/PMBB/pmbb_test_filter.yaml: synthetic --filter_config exercising every
filter key (drop_columns, keep_columns, drop_flags, restrictions,
restrict_separated, restrict_unknown). Expected filter_log.tsv is in its header.

test_data/PMBB/filter/: frozen outputs of
  nextflow run main.nf -profile test --run_eda \
      --sample_list test_data/PMBB/filter/roster.txt \
      --filter_config test_data/PMBB/pmbb_test_filter.yaml
  roster.txt                 full P001-P005 sample roster (--sample_list)
  phenotype_table.wide.tsv   GATHER wide table
  phenotype_manifest.tsv     manifest (with EDA flags column)
  eda_covariates.tsv         EDA covariate level codes (SEX: 0 female, 1 male;
                             COHORT: 1 A, 2 B, 3 C; SEX&COHORT: e.g. 0&1)
  eda_summary.tsv            per-column EDA metrics, flags, sep_* separation
Regenerate them with the command above if EDA/GATHER output formats change.
