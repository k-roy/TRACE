#!/usr/bin/env python3
"""
Create unified comparison table combining TRACE and CRISPResso2 results.

Aggregates results by biological sample (plate + well), calculating mean ± SEM
across the two primer pairs (KR2476 and KR2478).

Outputs samples in manifest order.

Author: Kevin R. Roy
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

# QC threshold for flagging low-read samples
MIN_READS = 1000

def parse_crispresso2_results(sample_id, crispresso_dir):
    """
    Parse CRISPResso2 quantification file for a single sample.

    Handles both output formats:
    - CRISPRessoBatch: crispresso_dir/CRISPResso_on_{sample_id}/...
    - Individual runs: crispresso_dir/{sample_id}/CRISPResso_on_{sample_id}/...

    Returns dict with:
        - crispresso_total_reads
        - crispresso_modified_pct
        - crispresso_unmodified_pct
        - crispresso_hdr_pct
        - crispresso_nhej_pct
    """
    # Try CRISPRessoBatch format first (new default)
    quant_file = crispresso_dir / f"CRISPResso_on_{sample_id}" / "CRISPResso_quantification_of_editing_frequency.txt"

    # Fall back to individual run format (legacy)
    if not quant_file.exists():
        quant_file = crispresso_dir / sample_id / f"CRISPResso_on_{sample_id}" / "CRISPResso_quantification_of_editing_frequency.txt"

    if not quant_file.exists():
        return {
            'crispresso_total_reads': np.nan,
            'crispresso_modified_pct': np.nan,
            'crispresso_unmodified_pct': np.nan,
            'crispresso_hdr_pct': np.nan,
            'crispresso_nhej_pct': np.nan
        }

    try:
        df = pd.read_csv(quant_file, sep='\t')

        # Get Reference row (overall editing)
        ref_row = df[df['Amplicon'] == 'Reference'].iloc[0]
        modified_pct = ref_row['Modified%']
        unmodified_pct = ref_row['Unmodified%']
        total_reads = ref_row['Reads_aligned']

        # Get HDR row (if exists)
        hdr_rows = df[df['Amplicon'] == 'HDR']
        if len(hdr_rows) > 0:
            hdr_row = hdr_rows.iloc[0]
            hdr_reads = hdr_row['Reads_aligned']
            hdr_pct = (hdr_reads / total_reads) * 100 if total_reads > 0 else 0.0
        else:
            hdr_pct = 0.0

        # Calculate NHEJ as modified - HDR
        nhej_pct = max(0.0, modified_pct - hdr_pct)

        return {
            'crispresso_total_reads': total_reads,
            'crispresso_modified_pct': modified_pct,
            'crispresso_unmodified_pct': unmodified_pct,
            'crispresso_hdr_pct': hdr_pct,
            'crispresso_nhej_pct': nhej_pct
        }
    except Exception as e:
        print(f"Warning: Error parsing CRISPResso2 results for {sample_id}: {e}", file=sys.stderr)
        return {
            'crispresso_total_reads': np.nan,
            'crispresso_modified_pct': np.nan,
            'crispresso_unmodified_pct': np.nan,
            'crispresso_hdr_pct': np.nan,
            'crispresso_nhej_pct': np.nan
        }


def calculate_mean_sem(values):
    """Calculate mean and SEM, handling NaN values."""
    clean_values = [v for v in values if not pd.isna(v)]
    if len(clean_values) == 0:
        return np.nan, np.nan
    elif len(clean_values) == 1:
        return clean_values[0], 0.0
    else:
        mean = np.mean(clean_values)
        sem = np.std(clean_values, ddof=1) / np.sqrt(len(clean_values))
        return mean, sem


def main():
    # Parse command line arguments
    manifest_file = Path(snakemake.input.manifest)
    trace_results = Path(snakemake.input.trace_results)
    crispresso_dir = Path(snakemake.params.crispresso_dir)
    output_file = Path(snakemake.output.unified_table)

    print(f"Loading manifest: {manifest_file}")
    manifest = pd.read_csv(manifest_file, sep='\t')

    print(f"Loading TRACE results: {trace_results}")
    trace_df = pd.read_csv(trace_results, sep='\t')

    # Create biological sample identifier (plate + well)
    manifest['bio_sample'] = manifest['plate'] + '_' + manifest['well']

    # Get unique biological samples in manifest order
    bio_samples_ordered = manifest.groupby('bio_sample', sort=False).first().reset_index()

    print(f"Found {len(bio_samples_ordered)} unique biological samples")
    print(f"Processing CRISPResso2 results from: {crispresso_dir}")

    # Build results table
    results = []

    for _, bio_sample_row in bio_samples_ordered.iterrows():
        bio_sample = bio_sample_row['bio_sample']
        plate = bio_sample_row['plate']
        well = bio_sample_row['well']

        # Get all technical replicates for this biological sample (KR2476 and KR2478)
        sample_rows = manifest[manifest['bio_sample'] == bio_sample]
        sample_ids = sample_rows['sample_id'].tolist()

        # Collect TRACE metrics for all replicates
        trace_hdr_pct = []
        trace_nhej_pct = []
        trace_wt_pct = []
        trace_donor_capture_pct = []
        trace_total_reads = []

        # Collect CRISPResso2 metrics for all replicates
        crispresso_hdr_pct = []
        crispresso_nhej_pct = []
        crispresso_unmodified_pct = []
        crispresso_total_reads = []

        # Track which replicates have sufficient reads
        replicate_reads = []
        replicate_data = []

        for sample_id in sample_ids:
            # Get TRACE results
            trace_row = trace_df[trace_df['sample_id'] == sample_id]
            if len(trace_row) > 0:
                trace_row = trace_row.iloc[0]
                reads = trace_row['total_reads']
                replicate_reads.append(reads)
                # Get donor_capture_pct if available (may not exist in older results)
                dc_pct = trace_row.get('donor_capture_pct', 0.0) if 'donor_capture_pct' in trace_row.index else 0.0
                replicate_data.append({
                    'sample_id': sample_id,
                    'hdr_pct': trace_row['HDR_pct'],
                    'nhej_pct': trace_row['NHEJ_pct'],
                    'wt_pct': trace_row['WT_pct'],
                    'donor_capture_pct': dc_pct,
                    'total_reads': reads,
                    'low_reads': reads < MIN_READS
                })
                trace_hdr_pct.append(trace_row['HDR_pct'])
                trace_nhej_pct.append(trace_row['NHEJ_pct'])
                trace_wt_pct.append(trace_row['WT_pct'])
                trace_donor_capture_pct.append(dc_pct)
                trace_total_reads.append(reads)
            else:
                trace_hdr_pct.append(np.nan)
                trace_nhej_pct.append(np.nan)
                trace_wt_pct.append(np.nan)
                trace_donor_capture_pct.append(np.nan)
                trace_total_reads.append(np.nan)

            # Get CRISPResso2 results
            crispresso_data = parse_crispresso2_results(sample_id, crispresso_dir)
            crispresso_hdr_pct.append(crispresso_data['crispresso_hdr_pct'])
            crispresso_nhej_pct.append(crispresso_data['crispresso_nhej_pct'])
            crispresso_unmodified_pct.append(crispresso_data['crispresso_unmodified_pct'])
            crispresso_total_reads.append(crispresso_data['crispresso_total_reads'])

        # Determine QC flag based on replicate read counts
        n_replicates_with_data = len(replicate_data)
        n_low_reads = sum(1 for r in replicate_data if r['low_reads'])

        if n_replicates_with_data == 0:
            quality_flag = "no_data"
            n_replicates_used = 0
        elif n_low_reads == 0:
            # All replicates have sufficient reads
            quality_flag = "good"
            n_replicates_used = n_replicates_with_data
        elif n_low_reads == n_replicates_with_data:
            # All replicates have low reads - use them anyway but flag
            quality_flag = "all_low_reads"
            n_replicates_used = n_replicates_with_data
        else:
            # Some replicates have low reads - exclude them from mean calculation
            quality_flag = "one_low_read_removed"
            good_replicates = [r for r in replicate_data if not r['low_reads']]
            n_replicates_used = len(good_replicates)
            # Recalculate means using only good replicates
            trace_hdr_pct = [r['hdr_pct'] for r in good_replicates]
            trace_nhej_pct = [r['nhej_pct'] for r in good_replicates]
            trace_wt_pct = [r['wt_pct'] for r in good_replicates]
            trace_donor_capture_pct = [r['donor_capture_pct'] for r in good_replicates]
            trace_total_reads = [r['total_reads'] for r in good_replicates]

        # Calculate mean ± SEM for each metric
        trace_hdr_mean, trace_hdr_sem = calculate_mean_sem(trace_hdr_pct)
        trace_nhej_mean, trace_nhej_sem = calculate_mean_sem(trace_nhej_pct)
        trace_wt_mean, trace_wt_sem = calculate_mean_sem(trace_wt_pct)
        trace_donor_capture_mean, trace_donor_capture_sem = calculate_mean_sem(trace_donor_capture_pct)
        trace_reads_mean, _ = calculate_mean_sem(trace_total_reads)

        crispresso_hdr_mean, crispresso_hdr_sem = calculate_mean_sem(crispresso_hdr_pct)
        crispresso_nhej_mean, crispresso_nhej_sem = calculate_mean_sem(crispresso_nhej_pct)
        crispresso_unmodified_mean, crispresso_unmodified_sem = calculate_mean_sem(crispresso_unmodified_pct)
        crispresso_reads_mean, _ = calculate_mean_sem(crispresso_total_reads)

        # Get sample description and category from manifest (first replicate)
        sample_description = bio_sample_row.get('sample_description', '')
        category = bio_sample_row.get('category', '')

        results.append({
            'bio_sample': bio_sample,
            'plate': plate,
            'well': well,
            'sample_description': sample_description,
            'category': category,
            'n_replicates': len(sample_ids),
            'n_replicates_used': n_replicates_used,
            'quality_flag': quality_flag,
            'trace_total_reads': trace_reads_mean,
            'trace_hdr_pct_mean': trace_hdr_mean,
            'trace_hdr_pct_sem': trace_hdr_sem,
            'trace_nhej_pct_mean': trace_nhej_mean,
            'trace_nhej_pct_sem': trace_nhej_sem,
            'trace_wt_pct_mean': trace_wt_mean,
            'trace_wt_pct_sem': trace_wt_sem,
            'trace_donor_capture_pct_mean': trace_donor_capture_mean,
            'trace_donor_capture_pct_sem': trace_donor_capture_sem,
            'crispresso_total_reads': crispresso_reads_mean,
            'crispresso_hdr_pct_mean': crispresso_hdr_mean,
            'crispresso_hdr_pct_sem': crispresso_hdr_sem,
            'crispresso_nhej_pct_mean': crispresso_nhej_mean,
            'crispresso_nhej_pct_sem': crispresso_nhej_sem,
            'crispresso_unmodified_pct_mean': crispresso_unmodified_mean,
            'crispresso_unmodified_pct_sem': crispresso_unmodified_sem,
        })

    # Create DataFrame and save
    results_df = pd.DataFrame(results)

    print(f"Writing unified comparison table: {output_file}")
    print(f"  Biological samples: {len(results_df)}")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')

    print("Done!")


if __name__ == '__main__':
    main()
