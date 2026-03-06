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
    output_dir = Path(snakemake.params.output_dir)

    results = []
    for sample_id in samples:
        cls_file = output_dir / sample_id / "classification.tsv"

        if cls_file.exists():
            df = pd.read_csv(cls_file, sep='\t')

            if len(df) == 0:
                # Empty file (control sample or error)
                results.append({
                    'sample_id': sample_id,
                    'total_reads': 0,
                    'HDR_PERFECT': 0,
                    'HDR_IMPERFECT': 0,
                    'HDR_total': 0,
                    'PARTIAL_HDR': 0,
                    'NHEJ_DELETION': 0,
                    'NHEJ_INSERTION': 0,
                    'NHEJ_total': 0,
                    'WILD_TYPE': 0,
                    'LARGE_DELETION': 0,
                    'UNCLASSIFIED': 0,
                    'HDR_pct': 0.0,
                    'NHEJ_pct': 0.0,
                    'WT_pct': 0.0,
                    'PARTIAL_HDR_pct': 0.0
                })
            else:
                # Count outcomes (using correct outcome names from classifier)
                outcome_counts = df['outcome'].value_counts().to_dict()
                total = len(df)

                # HDR = HDR_PERFECT + HDR_IMPERFECT
                hdr_perfect = outcome_counts.get('HDR_PERFECT', 0)
                hdr_imperfect = outcome_counts.get('HDR_IMPERFECT', 0)
                hdr_total = hdr_perfect + hdr_imperfect

                # NHEJ = NHEJ_DELETION + NHEJ_INSERTION
                nhej_del = outcome_counts.get('NHEJ_DELETION', 0)
                nhej_ins = outcome_counts.get('NHEJ_INSERTION', 0)
                nhej_total = nhej_del + nhej_ins

                # WT = WILD_TYPE
                wt = outcome_counts.get('WILD_TYPE', 0)

                # Large deletion
                large_del = outcome_counts.get('LARGE_DELETION', 0)

                # Unclassified
                unclassified = outcome_counts.get('UNCLASSIFIED', 0)

                # Partial HDR: reads with any hdr_match_fraction > 0 but not classified as HDR
                # This captures reads where some donor SNVs are incorporated
                partial_hdr = 0
                if 'hdr_match_fraction' in df.columns:
                    partial_hdr = len(df[(df['hdr_match_fraction'] > 0) &
                                          (~df['outcome'].isin(['HDR_PERFECT', 'HDR_IMPERFECT']))])

                results.append({
                    'sample_id': sample_id,
                    'total_reads': total,
                    'HDR_PERFECT': hdr_perfect,
                    'HDR_IMPERFECT': hdr_imperfect,
                    'HDR_total': hdr_total,
                    'PARTIAL_HDR': partial_hdr,
                    'NHEJ_DELETION': nhej_del,
                    'NHEJ_INSERTION': nhej_ins,
                    'NHEJ_total': nhej_total,
                    'WILD_TYPE': wt,
                    'LARGE_DELETION': large_del,
                    'UNCLASSIFIED': unclassified,
                    'HDR_pct': 100.0 * hdr_total / total if total > 0 else 0.0,
                    'NHEJ_pct': 100.0 * nhej_total / total if total > 0 else 0.0,
                    'WT_pct': 100.0 * wt / total if total > 0 else 0.0,
                    'PARTIAL_HDR_pct': 100.0 * partial_hdr / total if total > 0 else 0.0
                })

    # Write summary
    summary_df = pd.DataFrame(results)
    output_path = Path(snakemake.output.tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
