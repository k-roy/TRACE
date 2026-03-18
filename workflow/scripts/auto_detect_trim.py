#!/usr/bin/env python3
"""
Automatically detect optimal primer trimming by comparing reads from both primer pairs.

Strategy:
1. Sample reads from representative samples of each primer pair
2. Align reads to find the common core sequence
3. Determine minimal trim lengths that produce matching sequences

This enables automatic 2x cache efficiency by sharing sequences across primer pairs.
"""

import gzip
import json
from pathlib import Path
from collections import Counter

import pandas as pd


def sample_reads(fastq_path, n_reads=500):
    """Sample first n reads from a FASTQ file."""
    reads = []
    try:
        with gzip.open(fastq_path, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:  # Sequence line
                    reads.append(line.strip().upper())
                    if len(reads) >= n_reads:
                        break
    except Exception as e:
        print(f"Warning: Could not read {fastq_path}: {e}")
    return reads


def find_longest_common_substring(s1, s2, min_len=50):
    """Find the longest common substring between two sequences."""
    # Use suffix array approach for efficiency
    best_match = ""

    for i in range(len(s1) - min_len + 1):
        for length in range(min_len, min(len(s1) - i, len(s2)) + 1):
            substr = s1[i:i+length]
            if substr in s2:
                if len(substr) > len(best_match):
                    best_match = substr
                    pos1 = i
                    pos2 = s2.find(substr)

    if best_match:
        return best_match, pos1, pos2
    return None, -1, -1


def detect_optimal_trims(reads1, reads2, primer_pair1, primer_pair2):
    """
    Compare reads from two primer pairs to find optimal trim lengths.

    Returns dict with trim configurations.
    """
    # Get most common reads from each primer pair
    counter1 = Counter(reads1)
    counter2 = Counter(reads2)

    top_reads1 = [r for r, _ in counter1.most_common(10)]
    top_reads2 = [r for r, _ in counter2.most_common(10)]

    # Find common core by comparing top reads
    best_trim1 = 0
    best_trim2 = 0
    best_match_len = 0

    for r1 in top_reads1[:3]:
        for r2 in top_reads2[:3]:
            common, pos1, pos2 = find_longest_common_substring(r1, r2)
            if common and len(common) > best_match_len:
                best_match_len = len(common)
                best_trim1 = pos1
                best_trim2 = pos2

    # Detect UMI by looking at prefix diversity
    def estimate_umi_length(reads, max_umi=10):
        for umi_len in range(max_umi, 0, -1):
            prefixes = [r[:umi_len] for r in reads if len(r) >= umi_len]
            if not prefixes:
                continue
            unique_ratio = len(set(prefixes)) / len(prefixes)
            # High diversity suggests UMI
            if unique_ratio > 0.7:
                return umi_len
        return 0

    umi1 = estimate_umi_length(reads1)
    umi2 = estimate_umi_length(reads2)

    return {
        primer_pair1: {
            'umi_length': umi1,
            'r1_5prime_trim': best_trim1,
            'detected_primer_len': best_trim1 - umi1,
            'common_core_start': best_trim1
        },
        primer_pair2: {
            'umi_length': umi2,
            'r1_5prime_trim': best_trim2,
            'detected_primer_len': best_trim2 - umi2,
            'common_core_start': best_trim2
        },
        'common_core_length': best_match_len,
        'estimated_cache_sharing': f"{min(len(set(reads1)), len(set(reads2)))} sequences could be shared"
    }


def main():
    """
    Snakemake script to auto-detect trimming parameters.
    """
    manifest_path = snakemake.input.manifest
    output_path = snakemake.output.trim_config
    log_path = Path(str(snakemake.log))

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(manifest_path, sep='\t')

    # Group by primer pair
    if 'primer_pair' not in manifest.columns:
        # Try to detect from sample_id
        manifest['primer_pair'] = manifest['sample_id'].apply(
            lambda x: 'KR2478' if 'KR2478' in x else ('KR2476' if 'KR2476' in x else 'unknown')
        )

    primer_pairs = manifest['primer_pair'].unique()

    if len(primer_pairs) < 2:
        # Single primer pair - just detect UMI
        sample = manifest.iloc[0]
        reads = sample_reads(sample['r1_path'])
        umi_len = 0
        for ul in range(10, 0, -1):
            prefixes = [r[:ul] for r in reads if len(r) >= ul]
            if len(set(prefixes)) / len(prefixes) > 0.7:
                umi_len = ul
                break

        config = {
            primer_pairs[0]: {
                'umi_length': umi_len,
                'r1_5prime_trim': umi_len + 20,
                'r2_5prime_trim': 20
            }
        }
    else:
        # Multiple primer pairs - find optimal alignment
        reads_by_pp = {}
        for pp in primer_pairs:
            pp_samples = manifest[manifest['primer_pair'] == pp]
            if len(pp_samples) > 0:
                sample = pp_samples.iloc[0]
                reads_by_pp[pp] = sample_reads(sample['r1_path'])

        # Compare pairs to find optimal trims
        pp_list = list(reads_by_pp.keys())
        config = detect_optimal_trims(
            reads_by_pp[pp_list[0]],
            reads_by_pp[pp_list[1]],
            pp_list[0],
            pp_list[1]
        )

    # Save config
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)

    # Log results
    log_msg = f"Auto-detected trim configuration:\n{json.dumps(config, indent=2)}\n"
    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
