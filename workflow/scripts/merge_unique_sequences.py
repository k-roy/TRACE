#!/usr/bin/env python3
"""
Merge unique sequences from all samples into a global cache.

Creates:
1. global_unique_sequences.tsv.gz - All unique (seq, guide, donor) combinations
2. sequence_to_sample_mapping.tsv.gz - Maps each unique seq to samples that contain it

This is where the major savings happen: if the same sequence appears in 100 samples
with the same guide/donor, we only classify it once.
"""

import gzip
from pathlib import Path
from collections import defaultdict

import pandas as pd


def main():
    collapsed_files = snakemake.input.collapsed
    global_unique_path = snakemake.output.global_unique
    mapping_path = snakemake.output.sample_mapping
    log_path = Path(str(snakemake.log))

    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Track unique sequences globally
    # Key: (merged_seq, guide, donor) -> global count
    global_seqs = defaultdict(int)

    # Track which samples contain each sequence
    # Key: seq_hash -> list of (sample_id, count)
    sample_mapping = defaultdict(list)

    total_sample_seqs = 0

    for collapsed_file in collapsed_files:
        # Extract sample_id from path
        # Path format: .../samples/{sample_id}/collapsed/unique_sequences.tsv.gz
        sample_id = Path(collapsed_file).parent.parent.name

        df = pd.read_csv(collapsed_file, sep='\t', compression='gzip')
        total_sample_seqs += len(df)

        for _, row in df.iterrows():
            key = (row['merged_seq'], row['guide'], row['donor'])
            global_seqs[key] += row['count']
            sample_mapping[row['seq_hash']].append({
                'sample_id': sample_id,
                'count': row['count']
            })

    # Create global unique sequences table
    global_rows = []
    for (seq, guide, donor), total_count in global_seqs.items():
        import hashlib
        seq_hash = hashlib.md5(f"{seq}|{guide}|{donor}".encode()).hexdigest()
        global_rows.append({
            'seq_hash': seq_hash,
            'merged_seq': seq,
            'guide': guide,
            'donor': donor,
            'total_count': total_count
        })

    global_df = pd.DataFrame(global_rows)

    # Create sample mapping table
    mapping_rows = []
    for seq_hash, samples in sample_mapping.items():
        for sample_info in samples:
            mapping_rows.append({
                'seq_hash': seq_hash,
                'sample_id': sample_info['sample_id'],
                'count': sample_info['count']
            })

    mapping_df = pd.DataFrame(mapping_rows)

    # Save outputs
    Path(global_unique_path).parent.mkdir(parents=True, exist_ok=True)
    global_df.to_csv(global_unique_path, sep='\t', index=False, compression='gzip')
    mapping_df.to_csv(mapping_path, sep='\t', index=False, compression='gzip')

    # Calculate savings
    n_samples = len(collapsed_files)
    n_unique_per_sample = total_sample_seqs
    n_global_unique = len(global_seqs)
    compression = n_unique_per_sample / n_global_unique if n_global_unique > 0 else 1

    log_msg = (
        f"Merged unique sequences:\n"
        f"  Number of samples: {n_samples}\n"
        f"  Total unique seqs per sample: {n_unique_per_sample}\n"
        f"  Global unique seqs: {n_global_unique}\n"
        f"  Cross-sample compression: {compression:.1f}x\n"
        f"  \n"
        f"  This means we only need to classify {n_global_unique} sequences\n"
        f"  instead of {n_unique_per_sample} (saved {n_unique_per_sample - n_global_unique} classifications)\n"
    )
    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
