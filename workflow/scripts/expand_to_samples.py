#!/usr/bin/env python3
"""
Expand Level 2 classifications to per-sample results.

Looks up each sample's sequences in the merged classification cache,
filtering by the sample's specific (guide, donor) combination.

NEW: Supports paired-end format with sequencing mode, coverage info,
and insert size metadata from paired alignments.
"""

from pathlib import Path
from collections import defaultdict

import pandas as pd


def main():
    classifications_path = snakemake.input.classifications
    collapsed_path = snakemake.input.collapsed
    output_classification = snakemake.output.classification
    output_hdr_detail = snakemake.output.hdr_detail
    log_path = Path(str(snakemake.log))

    sample_id = snakemake.params.sample_id
    sample_guide = snakemake.params.guide if snakemake.params.guide else ""
    sample_donor = snakemake.params.donor if snakemake.params.donor else ""

    # Clean up nan values
    if sample_guide == "nan":
        sample_guide = ""
    if sample_donor == "nan":
        sample_donor = ""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_classification).parent.mkdir(parents=True, exist_ok=True)

    # Load Level 2 classifications (already filtered to this sample's combo)
    cache_df = pd.read_csv(classifications_path, sep='\t', compression='gzip')

    # Build lookup (no filtering needed - file already has only this combo)
    cache_lookup = {row['seq_hash']: row for _, row in cache_df.iterrows()}

    # Load sample's collapsed sequences
    sample_df = pd.read_csv(collapsed_path, sep='\t', compression='gzip')

    # Expand to per-sample results
    results = []
    per_snv_counts = defaultdict(int)
    total_hdr_reads = 0
    found = 0
    not_found = 0

    for _, row in sample_df.iterrows():
        seq_hash = row['seq_hash']
        count = row['count']

        if seq_hash not in cache_lookup:
            not_found += 1
            results.append({
                'read_name': f"{seq_hash}_x{count}",
                'outcome': 'NOT_IN_CACHE',
                'n_donor_snvs': 0,
                'n_non_donor': 0,
                'donor_fraction': 0.0,
                'donor_snvs_detected': '',
                'sequencing_mode': 'unknown',
                'coverage_mode': 'unknown',
                'uncovered_positions': '',
                'r1_coverage_start': -1,
                'r1_coverage_end': -1,
                'r2_coverage_start': -1,
                'r2_coverage_end': -1,
                'insert_size': 0,
                'is_proper_pair': False,
                'is_donor_capture': False,
                'donor_capture_size': 0,
                'donor_capture_distance': 0,
                'junction_homology_left': 0,
                'junction_homology_right': 0,
                'mmej_class': '',
                'count': count
            })
            continue

        found += 1
        cached = cache_lookup[seq_hash]

        results.append({
            'read_name': f"{seq_hash}_x{count}",
            'outcome': cached['outcome'],
            'n_donor_snvs': cached['n_donor_snvs'],
            'n_non_donor': cached['n_non_donor'],
            'donor_fraction': cached['donor_fraction'],
            'donor_snvs_detected': cached['donor_snvs_detected'] if pd.notna(cached['donor_snvs_detected']) else '',
            'sequencing_mode': cached.get('sequencing_mode', 'unknown'),
            'coverage_mode': cached.get('coverage_mode', 'unknown'),
            'uncovered_positions': cached.get('uncovered_positions', '') if pd.notna(cached.get('uncovered_positions')) else '',
            'r1_coverage_start': cached.get('r1_coverage_start', -1),
            'r1_coverage_end': cached.get('r1_coverage_end', -1),
            'r2_coverage_start': cached.get('r2_coverage_start', -1),
            'r2_coverage_end': cached.get('r2_coverage_end', -1),
            'insert_size': cached.get('insert_size', 0),
            'is_proper_pair': cached.get('is_proper_pair', False),
            'is_donor_capture': cached.get('is_donor_capture', False),
            'donor_capture_size': cached.get('donor_capture_size', 0),
            'donor_capture_distance': cached.get('donor_capture_distance', 0),
            'junction_homology_left': cached.get('junction_homology_left', 0),
            'junction_homology_right': cached.get('junction_homology_right', 0),
            'mmej_class': cached.get('mmej_class', ''),
            'count': count
        })

        # Track per-SNV integration for HDR reads
        if cached['outcome'] in ['HDR_COMPLETE', 'HDR_PARTIAL', 'MIXED']:
            total_hdr_reads += count
            snvs = cached['donor_snvs_detected']
            if pd.notna(snvs) and snvs:
                for pos in str(snvs).split(','):
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
    total_reads = results_df['count'].sum()

    outcome_counts = {}
    for outcome in results_df['outcome'].unique():
        mask = results_df['outcome'] == outcome
        outcome_counts[outcome] = results_df.loc[mask, 'count'].sum()

    # Calculate sequencing mode distribution
    mode_counts = {}
    if 'sequencing_mode' in results_df.columns:
        for mode in results_df['sequencing_mode'].unique():
            mask = results_df['sequencing_mode'] == mode
            mode_counts[mode] = results_df.loc[mask, 'count'].sum()

    # Calculate coverage mode distribution
    coverage_counts = {}
    if 'coverage_mode' in results_df.columns:
        for cov_mode in results_df['coverage_mode'].unique():
            mask = results_df['coverage_mode'] == cov_mode
            coverage_counts[cov_mode] = results_df.loc[mask, 'count'].sum()

    # Write log
    log_msg = (
        f"Expanded classifications for {sample_id}:\n"
        f"  Guide: {sample_guide[:30]}...\n"
        f"  Donor: {sample_donor[:30]}...\n"
        f"  \n"
        f"  Unique sequences: {len(sample_df)}\n"
        f"  Found in cache: {found}\n"
        f"  Not in cache: {not_found}\n"
        f"  Total reads: {total_reads}\n"
        f"  \n"
    )

    # Add sequencing mode info
    if mode_counts:
        log_msg += f"Sequencing modes:\n"
        for mode, count in sorted(mode_counts.items(), key=lambda x: -x[1]):
            pct = 100 * count / total_reads if total_reads > 0 else 0
            log_msg += f"  {mode}: {count} ({pct:.1f}%)\n"
        log_msg += f"  \n"

    # Add coverage info
    if coverage_counts:
        log_msg += f"Coverage modes:\n"
        for cov_mode, count in sorted(coverage_counts.items(), key=lambda x: -x[1]):
            pct = 100 * count / total_reads if total_reads > 0 else 0
            log_msg += f"  {cov_mode}: {count} ({pct:.1f}%)\n"
        log_msg += f"  \n"

    log_msg += f"Outcome distribution:\n"

    for outcome, count in sorted(outcome_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / total_reads if total_reads > 0 else 0
        log_msg += f"  {outcome}: {count} ({pct:.1f}%)\n"

    # Add donor capture statistics
    if 'is_donor_capture' in results_df.columns:
        donor_captures = results_df[results_df['is_donor_capture'] == True]
        if len(donor_captures) > 0:
            capture_count = donor_captures['count'].sum()
            capture_pct = 100 * capture_count / total_reads if total_reads > 0 else 0
            avg_size = donor_captures['donor_capture_size'].mean()
            avg_distance = donor_captures['donor_capture_distance'].mean()
            log_msg += f"  \n"
            log_msg += f"Donor capture detection:\n"
            log_msg += f"  Reads with donor capture: {capture_count} ({capture_pct:.1f}%)\n"
            log_msg += f"  Average insertion size: {avg_size:.0f}bp\n"
            log_msg += f"  Average distance from cut site: {avg_distance:.0f}bp\n"

            # Junction microhomology analysis (NHEJ vs MMEJ classification)
            if 'mmej_class' in donor_captures.columns:
                log_msg += f"  \n"
                log_msg += f"Junction microhomology analysis:\n"
                for mmej_class, group in donor_captures.groupby('mmej_class'):
                    class_count = group['count'].sum()
                    class_pct = 100 * class_count / capture_count
                    avg_homology_l = group['junction_homology_left'].mean()
                    avg_homology_r = group['junction_homology_right'].mean()
                    log_msg += f"    {mmej_class}: {class_count} ({class_pct:.1f}%) | "
                    log_msg += f"Avg homology: {avg_homology_l:.1f}bp (left), {avg_homology_r:.1f}bp (right)\n"

    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
