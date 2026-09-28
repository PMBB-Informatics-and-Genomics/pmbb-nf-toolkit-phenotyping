#!/usr/bin/env python3
"""
phenotype_manifest.py — List the output columns of one phenotyping result.

Usage:
    phenotype_manifest.py --result BMI.quant.tsv --config BMI.json \
                          --output BMI.manifest.tsv

Writes one row per non-sample_ID column of --result:

    column_name | phenotype | data_type | code

  phenotype : YAML phenotype key (config 'phenotype_name'; chunked coded
              configs keep the parent key)
  data_type : quantitative | categorical | coded | advanced_coded
  code      : coded only — the raw code, i.e. column_name with column_prefix
              stripped; blank for all other data types

main.nf concatenates these into {output_dir}/phenotype_manifest.tsv.
"""

import argparse
import json
import sys

import pandas as pd

MANIFEST_COLUMNS = ['column_name', 'phenotype', 'data_type', 'code']


def run_phenotype_manifest(result_path, config_path, output_path):
    """Write manifest rows for every output column in result_path."""
    with open(config_path) as f:
        cfg = json.load(f)

    phenotype = cfg.get('phenotype_name', '')
    data_type = cfg.get('data_type', '')
    prefix = cfg.get('column_prefix') or ''

    columns = pd.read_csv(result_path, sep='\t', nrows=0).columns
    rows = []
    for col in columns:
        if col == 'sample_ID':
            continue
        code = ''
        if data_type == 'coded':
            code = col[len(prefix):] if prefix and col.startswith(prefix) else col
        rows.append([col, phenotype, data_type, code])

    pd.DataFrame(rows, columns=MANIFEST_COLUMNS).to_csv(output_path, sep='\t', index=False)
    print(f'Wrote: {output_path} ({len(rows)} columns)', file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--result', required=True, help='Per-phenotype result TSV')
    p.add_argument('--config', required=True, help='Per-phenotype JSON config')
    p.add_argument('--output', required=True, help='Output manifest TSV')
    args = p.parse_args()
    run_phenotype_manifest(args.result, args.config, args.output)


if __name__ == '__main__':
    main()
