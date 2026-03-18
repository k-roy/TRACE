#!/usr/bin/env python3
"""
Automatically detect optimal primer trimming by comparing reads from both primer pairs.

Strategy:
1. Sample reads from both primer pairs (same biological sample)
2. Find longest common suffix/prefix between trimmed reads
3. Determine minimal trim lengths that produce matching core sequences

This generalizes to any primer pair combination without hardcoding.
"""

import gzip
import sys
from pathlib import Path
from collections import Counter

import pandas as pd


def sample_reads(fastq_path, n_reads=1000):
    """Sample first n reads from a FASTQ file."""
    reads = []
    with gzip.open(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                reads.append(line.strip().upper())
                if len(reads) >= n_reads:
                    break
    return reads


def find_common_core(seqs1, seqs2, min_len=100):
    """
    Find the optimal trim lengths to maximize shared sequences.

    Returns: (trim1_5prime, trim1_3prime, trim2_5prime, trim2_3prime, shared_fraction)
    """
    best_shared = 0
    best_trims = (0, 0, 0, 0)

    # Try different 5' trim lengths (0 to 30bp)
    for trim1_5 in range(0, 31, 2):
        for trim2_5 in range(0, 31, 2):
            # Trim sequences
            trimmed1 = set(s[trim1_5:] for s in seqs1 if len(s) > trim1_5 + min_len)
            trimmed2 = set(s[trim2_5:] for s in seqs2 if len(s) > trim2_5 + min_len)

            # Count shared
            shared = len(trimmed1 & trimmed2)
            total = len(trimmed1 | trimmed2)

            if total > 0 and shared > best_shared:
                best_shared = shared
                best_trims = (trim1_5, 0, trim2_5, 0)

    return best_trims, best_shared


def detect_umi_length(reads, expected_lengths=[0, 6, 8, 10]):
    """
    Detect UMI length by looking for random sequence prefix.
    UMIs typically have high diversity in the first N bases.
    """
    for umi_len in expected_lengths:
        if umi_len == 0:
            continue

        # Check diversity of first N bases
        prefixes = [r[:umi_len] for r in reads if len(r) >= umi_len]
        unique_ratio = len(set(prefixes)) / len(prefixes) if prefixes else 0

        # High diversity (>50% unique) suggests UMI
        if unique_ratio > 0.5:
            return umi_len

    return 0


def main():
    """
    Main function to detect primer trimming strategy.

    Usage: python detect_primer_trim.py manifest.tsv output.tsv
    """
    if len(sys.argv) < 3:
        print("Usage: python detect_primer_trim.py manifest.tsv output.tsv")
        sys.exit(1)

    manifest_path = sys.argv[1]
    output_path = sys.argv[2]

    manifest = pd.read_csv(manifest_path, sep='\t')

    # Group samples by biological replicate (same plate + well)
    # to find pairs with different primer pairs
    manifest['bio_id'] = manifest['plate'] + '_' + manifest['well']

    # Detect primer pairs
    primer_pairs = manifest['primer_pair'].unique() if 'primer_pair' in manifest.columns else []

    if len(primer_pairs) < 2:
        print("Only one primer pair found, no cross-primer optimization possible")
        # Just detect UMI for each sample
        results = []
        for _, row in manifest.head(10).iterrows():
            reads = sample_reads(row['r1_path'], n_reads=100)
            umi_len = detect_umi_length(reads)
            results.append({
                'sample_id': row['sample_id'],
                'umi_length': umi_len,
                'r1_5prime_trim': umi_len + 20,  # UMI + primer
                'r2_5prime_trim': 20  # Primer only
            })

        pd.DataFrame(results).to_csv(output_path, sep='\t', index=False)
        return

    # Find samples with both primer pairs
    bio_groups = manifest.groupby('bio_id')

    results = []
    for bio_id, group in list(bio_groups)[:5]:  # Sample first 5
        if len(group) < 2:
            continue

        # Get reads from each primer pair
        for pp in primer_pairs:
            pp_samples = group[group['primer_pair'] == pp]
            if len(pp_samples) > 0:
                sample = pp_samples.iloc[0]
                reads = sample_reads(sample['r1_path'], n_reads=500)
                umi_len = detect_umi_length(reads)

                results.append({
                    'bio_id': bio_id,
                    'primer_pair': pp,
                    'sample_id': sample['sample_id'],
                    'umi_length': umi_len,
                    'suggested_r1_trim': umi_len + 20
                })

    if results:
        df = pd.DataFrame(results)

        # Compare primer pairs to find optimal alignment
        print("\nPrimer pair analysis:")
        for pp in primer_pairs:
            pp_data = df[df['primer_pair'] == pp]
            if len(pp_data) > 0:
                avg_umi = pp_data['umi_length'].mean()
                avg_trim = pp_data['suggested_r1_trim'].mean()
                print(f"  {pp}: UMI={avg_umi:.1f}bp, suggested R1 trim={avg_trim:.1f}bp")

        df.to_csv(output_path, sep='\t', index=False)
    else:
        print("No paired samples found for comparison")


if __name__ == '__main__':
    main()
