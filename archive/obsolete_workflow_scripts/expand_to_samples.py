#!/usr/bin/env python3
"""
Expand cached classifications to per-sample results.

Takes the global classification cache and creates per-sample output files
by looking up each sample's sequences in the cache.
"""

import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd

# Add TRACE to path
TRACE_ROOT = Path("/oak/stanford/groups/larsms/Users/kevinroy/software/trace")
sys.path.insert(0, str(TRACE_ROOT))


def main():
    classifications_path = snakemake.input.classifications
    collapsed_path = snakemake.input.collapsed
    output_classification = snakemake.output.classification
    output_hdr_detail = snakemake.output.hdr_detail
    log_path = Path(str(snakemake.log))
    sample_id = snakemake.params.sample_id

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_classification).parent.mkdir(parents=True, exist_ok=True)

    # Load cached classifications
    cache_df = pd.read_csv(classifications_path, sep='\t', compression='gzip')
    cache_lookup = {row['seq_hash']: row for _, row in cache_df.iterrows()}

    # Load sample's collapsed sequences
    sample_df = pd.read_csv(collapsed_path, sep='\t', compression='gzip')

    # Expand to per-sample results
    results = []
    per_snv_counts = defaultdict(int)
    total_hdr_reads = 0

    for _, row in sample_df.iterrows():
        seq_hash = row['seq_hash']
        count = row['count']

        if seq_hash not in cache_lookup:
            # Not in cache (shouldn't happen)
            results.append({
                'read_name': f"{seq_hash}_x{count}",
                'outcome': 'UNKNOWN',
                'n_donor_snvs': 0,
                'n_non_donor': 0,
                'donor_fraction': 0.0,
                'donor_snvs_detected': ''
            })
            continue

        cached = cache_lookup[seq_hash]

        # Add entry for this sequence (with count weight)
        results.append({
            'read_name': f"{seq_hash}_x{count}",
            'outcome': cached['outcome'],
            'n_donor_snvs': cached['n_donor_snvs'],
            'n_non_donor': cached['n_non_donor'],
            'donor_fraction': cached['donor_fraction'],
            'donor_snvs_detected': cached['donor_snvs_detected'],
            'count': count  # Store count for weighted aggregation
        })

        # Track per-SNV integration for HDR reads
        if cached['outcome'] in ['HDR_COMPLETE', 'HDR_PARTIAL', 'MIXED']:
            total_hdr_reads += count
            if cached['donor_snvs_detected']:
                for pos in cached['donor_snvs_detected'].split(','):
                    if pos:
                        per_snv_counts[int(pos)] += count

    # Save classification results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_classification, sep='\t', index=False)

    # Save HDR detail (per-SNV frequencies)
    hdr_detail_rows = []
    for pos, count in sorted(per_snv_counts.items()):
        freq = count / total_hdr_reads if total_hdr_reads > 0 else 0
        hdr_detail_rows.append({
            'position': pos,
            'count': count,
            'frequency': freq
        })

    hdr_detail_df = pd.DataFrame(hdr_detail_rows)
    hdr_detail_df.to_csv(output_hdr_detail, sep='\t', index=False)

    # Calculate summary statistics
    total_reads = results_df['count'].sum() if 'count' in results_df.columns else len(results_df)

    outcome_counts = {}
    for outcome in results_df['outcome'].unique():
        mask = results_df['outcome'] == outcome
        if 'count' in results_df.columns:
            outcome_counts[outcome] = results_df.loc[mask, 'count'].sum()
        else:
            outcome_counts[outcome] = mask.sum()

    # Write log
    log_msg = (
        f"Expanded classifications for {sample_id}:\n"
        f"  Unique sequences: {len(sample_df)}\n"
        f"  Total reads: {total_reads}\n"
        f"  \n"
        f"Outcome distribution:\n"
    )

    for outcome, count in sorted(outcome_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / total_reads if total_reads > 0 else 0
        log_msg += f"  {outcome}: {count} ({pct:.1f}%)\n"

    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
