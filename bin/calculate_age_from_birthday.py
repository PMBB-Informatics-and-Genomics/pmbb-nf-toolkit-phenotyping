"""
calculate_age_from_birthday.py
Calculate age from birthdate column, given a reference date.

# Usage:
    python calculate_age.py \
    --input demographics.tsv \
    --output demographics_with_age.tsv \
    --birth_col birth_datetime \
    --ref_date 2024-10-06 \
    --sep "\t"
"""

import pandas as pd
import argparse

def calculate_age(birthdates, ref_date):
    ref = pd.Timestamp(ref_date)
    b = pd.to_datetime(birthdates)

    return (
        ref.year
        - b.dt.year
        - (
            (ref.month < b.dt.month) |
            ((ref.month == b.dt.month) & (ref.day < b.dt.day))
        )
    )

def main():
    parser = argparse.ArgumentParser(description="Calculate age from birthdate column.")
    parser.add_argument("--input", required=True, help="Input table (csv/tsv)")
    parser.add_argument("--output", required=True, help="Output table")
    parser.add_argument("--birth_col", required=True, help="Birthdate column name")
    parser.add_argument("--ref_date", required=True, help="Reference date (YYYY-MM-DD)")
    parser.add_argument("--sep", default=",", help="Delimiter (default: , ; use '\\t' for TSV)")
    
    args = parser.parse_args()

    sep = "\t" if args.sep == "\\t" else args.sep

    df = pd.read_csv(args.input, sep=sep)

    df["age"] = calculate_age(df[args.birth_col], args.ref_date)

    df.to_csv(args.output, sep=sep, index=False)

if __name__ == "__main__":
    main()