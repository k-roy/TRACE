#!/usr/bin/env python3
"""
Create pooled summary table by combining technical replicates.

Takes a unified comparison table and groups samples by a biological sample ID,
computing mean, SEM, and sum statistics across technical replicates.

Usage:
    python create_pooled_summary.py \
        --input unified_comparison_table.tsv \
        --output pooled_summary.tsv \
        --sample-col sample_id \
        --bio-sample-pattern "plate_(\d+[A-Z])_sample_(\d+)" \
        [--metadata manifest.tsv]

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
"""

import argparse
import re
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional


def extract_bio_sample(sample_id: str, pattern: str) -> Optional[str]:
    """
    Extract biological sample identifier from sample_id using regex pattern.

    Default pattern groups technical replicates by plate+well.
    """
    match = re.search(pattern, sample_id)
    if match:
        return '_'.join(match.groups())
    return None


def compute_pooled_stats(df: pd.DataFrame,
                          bio_sample_col: str,
                          numeric_cols: list,
                          count_col: str = 'total_reads') -> pd.DataFrame:
    """
    Compute pooled statistics for each biological sample.

    For each numeric column:
    - _mean: weighted mean by total reads
    - _sem: standard error of the mean
    - _sum: sum across replicates
    """
    results = []

    for bio_sample, group in df.groupby(bio_sample_col):
        row = {
            bio_sample_col: bio_sample,
            'n_replicates': len(group)
        }

        # Sum of reads
        total_reads = group[count_col].sum()
        row[f'{count_col}_sum'] = total_reads

        for col in numeric_cols:
            if col == count_col:
                continue

            values = group[col].values
            weights = group[count_col].values

            # Check if this is a count column (sum makes sense) or percentage
            is_percentage = 'pct' in col.lower() or col.endswith('%')

            if is_percentage:
                # For percentages, compute weighted mean and SEM
                if total_reads > 0:
                    weighted_mean = np.average(values, weights=weights)
                else:
                    weighted_mean = 0

                # SEM of the weighted values
                if len(values) > 1:
                    sem = np.std(values, ddof=1) / np.sqrt(len(values))
                else:
                    sem = 0

                row[f'{col}_mean'] = round(weighted_mean, 4)
                row[f'{col}_sem'] = round(sem, 4)
            else:
                # For counts, compute sum
                row[f'{col}_sum'] = int(group[col].sum())

        results.append(row)

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description='Create pooled summary table by combining technical replicates'
    )
    parser.add_argument('--input', '-i', type=Path, required=True,
                        help='Input unified comparison table')
    parser.add_argument('--output', '-o', type=Path, required=True,
                        help='Output pooled summary table')
    parser.add_argument('--sample-col', default='sample_id',
                        help='Column name for sample identifiers')
    parser.add_argument('--bio-sample-pattern',
                        default=r'plate_(\d+[A-Z])_sample_(\d+)',
                        help='Regex pattern to extract biological sample ID')
    parser.add_argument('--bio-sample-col', default='bio_sample',
                        help='Column name for biological sample identifier')
    parser.add_argument('--metadata', type=Path,
                        help='Optional metadata TSV to merge in')
    parser.add_argument('--metadata-key', default='sample_id',
                        help='Column in metadata to join on')

    args = parser.parse_args()

    # Load data
    print(f"Loading {args.input}...")
    df = pd.read_csv(args.input, sep='\t')
    print(f"  {len(df)} samples loaded")

    # Extract biological sample ID
    print(f"Extracting biological sample IDs...")
    df[args.bio_sample_col] = df[args.sample_col].apply(
        lambda x: extract_bio_sample(x, args.bio_sample_pattern)
    )

    # Filter out samples without valid bio_sample
    missing = df[args.bio_sample_col].isna().sum()
    if missing > 0:
        print(f"  Warning: {missing} samples did not match pattern, will be excluded")
        df = df[df[args.bio_sample_col].notna()]

    n_bio = df[args.bio_sample_col].nunique()
    print(f"  {n_bio} unique biological samples")

    # Identify numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    numeric_cols = [c for c in numeric_cols if c != args.sample_col]

    # Compute pooled stats
    print("Computing pooled statistics...")
    pooled = compute_pooled_stats(
        df,
        args.bio_sample_col,
        numeric_cols,
        count_col='total_reads'
    )

    # Merge metadata if provided
    if args.metadata:
        print(f"Merging metadata from {args.metadata}...")
        meta = pd.read_csv(args.metadata, sep='\t')

        # Extract bio_sample from metadata too
        if args.metadata_key in meta.columns:
            meta[args.bio_sample_col] = meta[args.metadata_key].apply(
                lambda x: extract_bio_sample(x, args.bio_sample_pattern)
            )

            # Get unique metadata per bio_sample (first occurrence)
            meta_unique = meta.groupby(args.bio_sample_col).first().reset_index()

            # Select non-numeric columns for merge
            meta_cols = [args.bio_sample_col] + [
                c for c in meta_unique.columns
                if c not in numeric_cols and c != args.bio_sample_col
            ]
            meta_unique = meta_unique[meta_cols]

            pooled = pooled.merge(meta_unique, on=args.bio_sample_col, how='left')

    # Reorder columns: bio_sample, metadata, then stats
    cols = pooled.columns.tolist()
    first_cols = [args.bio_sample_col, 'n_replicates']
    other_cols = [c for c in cols if c not in first_cols]
    pooled = pooled[first_cols + other_cols]

    # Save
    print(f"Saving to {args.output}...")
    pooled.to_csv(args.output, sep='\t', index=False)
    print(f"  {len(pooled)} biological samples")

    # Summary
    print("\nSummary:")
    print(f"  Input samples: {len(df)}")
    print(f"  Output bio samples: {len(pooled)}")
    print(f"  Avg replicates: {pooled['n_replicates'].mean():.1f}")


if __name__ == '__main__':
    main()
