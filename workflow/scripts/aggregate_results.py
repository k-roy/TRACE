#!/usr/bin/env python3
"""
Aggregate per-sample classification results into summary table.

Called by Snakemake with snakemake object providing:
- snakemake.input: List of classification TSV files
- snakemake.output: Summary TSV file
- snakemake.params: Sample list
"""

import pandas as pd
from pathlib import Path


def main():
    samples = snakemake.params.samples
    samples_dir = Path(snakemake.params.samples_dir)

    results = []
    for sample_id in samples:
        cls_file = samples_dir / sample_id / "classification.tsv"

        if cls_file.exists():
            df = pd.read_csv(cls_file, sep='\t')

            if len(df) == 0:
                # Empty file (control sample or error)
                results.append({
                    'sample_id': sample_id,
                    'total_reads': 0,
                    'HDR_COMPLETE': 0,
                    'HDR_PARTIAL': 0,
                    'HDR_total': 0,
                    'NHEJ_INDEL': 0,
                    'MMEJ_INDEL': 0,
                    'NHEJ_MMEJ_indel': 0,
                    'HDR_PLUS_NHEJ_INDEL': 0,
                    'HDR_PLUS_MMEJ_INDEL': 0,
                    'HDR_PLUS_NHEJ_MMEJ_indel': 0,
                    'HDR_PLUS_OTHER': 0,
                    'NON_NHEJ_INDEL': 0,
                    'NON_DONOR_SNV': 0,
                    'DONOR_CAPTURE': 0,
                    'CHIMERIC': 0,
                    'WT': 0,
                    'UNALIGNED': 0,
                    'HDR_pct': 0.0,
                    'NHEJ_pct': 0.0,
                    'MMEJ_pct': 0.0,
                    'NHEJ_MMEJ_pct': 0.0,
                    'HDR_PLUS_NHEJ_pct': 0.0,
                    'HDR_PLUS_MMEJ_pct': 0.0,
                    'HDR_PLUS_NHEJ_MMEJ_pct': 0.0,
                    'HDR_PLUS_OTHER_pct': 0.0,
                    'WT_pct': 0.0,
                    'edited_pct': 0.0,
                    'off_target_pct': 0.0,
                    'chimeric_pct': 0.0,
                    'donor_capture_pct': 0.0
                })
            else:
                # Count outcomes (sum counts weighted by 'count' column if present)
                if 'count' in df.columns:
                    # Group by outcome and sum counts
                    outcome_counts = df.groupby('outcome')['count'].sum().to_dict()
                    total = df['count'].sum()
                else:
                    outcome_counts = df['outcome'].value_counts().to_dict()
                    total = len(df)

                # HDR outcomes from edit_distance_hdr classifier
                hdr_complete = outcome_counts.get('HDR_COMPLETE', 0)
                hdr_partial = outcome_counts.get('HDR_PARTIAL', 0)
                hdr_total = hdr_complete + hdr_partial

                # NHEJ outcomes (classical pathway: 0-2bp microhomology)
                nhej_indel = outcome_counts.get('NHEJ_INDEL', 0)
                hdr_plus_nhej = outcome_counts.get('HDR_PLUS_NHEJ_INDEL', 0)

                # MMEJ outcomes (alt-NHEJ pathway: >2bp microhomology)
                mmej_indel = outcome_counts.get('MMEJ_INDEL', 0)
                hdr_plus_mmej = outcome_counts.get('HDR_PLUS_MMEJ_INDEL', 0)

                # Mixed/other outcomes
                hdr_plus_other = outcome_counts.get('HDR_PLUS_OTHER', 0)

                # Off-target outcomes
                non_nhej_indel = outcome_counts.get('NON_DONOR_NON_NHEJ_INDEL', 0)
                non_donor_snv = outcome_counts.get('NON_DONOR_SNV', 0)

                # Artifacts
                chimeric = outcome_counts.get('CHIMERIC', 0)

                # Donor capture (NHEJ-mediated donor insertion)
                donor_capture = outcome_counts.get('DONOR_CAPTURE', 0)

                # Control outcomes
                wt = outcome_counts.get('WT', 0)
                unaligned = outcome_counts.get('UNALIGNED', 0)
                no_cigar = outcome_counts.get('NO_CIGAR', 0)

                # Combined NHEJ+MMEJ categories (for comparison with CRISPResso2)
                nhej_mmej_indel = nhej_indel + mmej_indel
                hdr_plus_nhej_mmej = hdr_plus_nhej + hdr_plus_mmej

                # Edited = HDR + NHEJ + MMEJ + mixed outcomes
                edited_total = hdr_total + nhej_mmej_indel + hdr_plus_nhej_mmej + hdr_plus_other

                results.append({
                    'sample_id': sample_id,
                    'total_reads': total,
                    # Core editing outcomes
                    'HDR_COMPLETE': hdr_complete,
                    'HDR_PARTIAL': hdr_partial,
                    'HDR_total': hdr_total,
                    'NHEJ_INDEL': nhej_indel,
                    'MMEJ_INDEL': mmej_indel,
                    'NHEJ_MMEJ_indel': nhej_mmej_indel,  # Combined for CRISPResso2 comparison
                    'HDR_PLUS_NHEJ_INDEL': hdr_plus_nhej,
                    'HDR_PLUS_MMEJ_INDEL': hdr_plus_mmej,
                    'HDR_PLUS_NHEJ_MMEJ_indel': hdr_plus_nhej_mmej,  # Combined for CRISPResso2 comparison
                    'HDR_PLUS_OTHER': hdr_plus_other,
                    # Off-target outcomes
                    'NON_NHEJ_INDEL': non_nhej_indel,
                    'NON_DONOR_SNV': non_donor_snv,
                    # Donor capture
                    'DONOR_CAPTURE': donor_capture,
                    # Artifacts and controls
                    'CHIMERIC': chimeric,
                    'WT': wt,
                    'UNALIGNED': unaligned,
                    # Percentages
                    'HDR_pct': 100.0 * hdr_total / total if total > 0 else 0.0,
                    'NHEJ_pct': 100.0 * nhej_indel / total if total > 0 else 0.0,
                    'MMEJ_pct': 100.0 * mmej_indel / total if total > 0 else 0.0,
                    'NHEJ_MMEJ_pct': 100.0 * nhej_mmej_indel / total if total > 0 else 0.0,  # Combined
                    'HDR_PLUS_NHEJ_pct': 100.0 * hdr_plus_nhej / total if total > 0 else 0.0,
                    'HDR_PLUS_MMEJ_pct': 100.0 * hdr_plus_mmej / total if total > 0 else 0.0,
                    'HDR_PLUS_NHEJ_MMEJ_pct': 100.0 * hdr_plus_nhej_mmej / total if total > 0 else 0.0,  # Combined
                    'HDR_PLUS_OTHER_pct': 100.0 * hdr_plus_other / total if total > 0 else 0.0,
                    'WT_pct': 100.0 * wt / total if total > 0 else 0.0,
                    'edited_pct': 100.0 * edited_total / total if total > 0 else 0.0,
                    'off_target_pct': 100.0 * (non_nhej_indel + non_donor_snv) / total if total > 0 else 0.0,
                    'chimeric_pct': 100.0 * chimeric / total if total > 0 else 0.0,
                    'donor_capture_pct': 100.0 * donor_capture / total if total > 0 else 0.0
                })

    # Write summary
    summary_df = pd.DataFrame(results)
    output_path = Path(snakemake.output.tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
