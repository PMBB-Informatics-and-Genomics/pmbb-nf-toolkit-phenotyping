# Community Configs

Production-validated configuration files from real biobank deployments.
Copy individual phenotype entries or whole files into your own biobank config.

| File | Contents | Source |
|---|---|---|
| `advanced_diagnostics.yaml` | ~220 ICD-10 advanced diagnostic phenotypes | PMBB + BraVa (Backman et al. 2021) |
| `event_mask.yaml` | Obstetric trimester event windows (first/second/third) | Penn Medicine BioBank |
| `risk_classification.yaml` | BMI, BP, cholesterol risk classification rules | Penn Medicine BioBank |

## Difference from `examples/`

`examples/` contains **annotated how-to templates** that explain config structure and options.

`community/` contains **real production configs** you can use directly or adapt for your biobank.

## Usage

Each file has usage instructions in its header comment. The general pattern:

```yaml
# In your biobank config (e.g. my_biobank_config.yaml):
phenotypes:
  ICD10:
    concept: diagnostic
    table: /path/to/your/icd10.tsv
    value_col: condition_source_value

  # Paste entries from advanced_diagnostics.yaml:
  T2Diab:
    concept: advanced_diagnostic
    sources:
      ICD10:
        from_phenotype: ICD10
        case_codes: [E11, O24.1]
        case_exclude: [E10, E12, E13]
```
