# biobanks/

Place your biobank configuration directory here. Each biobank gets its own
subdirectory containing a config YAML and any supporting files.

```
biobanks/
└── YOUR_BIOBANK/
    ├── config.yaml          # main phenotyping config (see examples/)
    ├── event_mask.yaml      # optional: obstetric/temporal event mask
    └── risk_classification  # optional: embedded in config.yaml
```

## Getting started

See `examples/` for annotated templates showing the full config schema.

See `community/` for production-tested configs you can use directly or adapt:
- `community/advanced_codeds.yaml` — 220+ advanced coded phenotypes
  (BraVa UK Biobank phenome-wide study + PMBB-validated additions)
- `community/event_mask.yaml` — obstetric event masking windows (ICD-10)
- `community/risk_classification.yaml` — risk classification rules for
  common quantitative phenotypes (BMI, BP, HDL, LDL)
