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
                    'HDR': 0,
                    'NHEJ': 0,
                    'WT': 0,
                    'LARGE_DELETION': 0,
                    'HDR_pct': 0.0
                })
            else:
                # Count outcomes
                outcome_counts = df['outcome'].value_counts().to_dict()
                total = len(df)
                hdr = outcome_counts.get('HDR', 0)

                results.append({
                    'sample_id': sample_id,
                    'total_reads': total,
                    'HDR': hdr,
                    'NHEJ': outcome_counts.get('NHEJ', 0),
                    'WT': outcome_counts.get('WT', 0),
                    'LARGE_DELETION': outcome_counts.get('LARGE_DELETION', 0),
                    'HDR_pct': 100.0 * hdr / total if total > 0 else 0.0
                })

    # Write summary
    summary_df = pd.DataFrame(results)
    output_path = Path(snakemake.output.tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
