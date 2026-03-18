#!/usr/bin/env python3
"""
Aggregate CRISPResso2 results from per-sample runs.

Called by Snakemake with snakemake object providing:
- snakemake.params.samples: List of sample IDs
- snakemake.params.samples_dir: Directory containing per-sample results
- snakemake.output.tsv: Output summary TSV file
"""

import pandas as pd
from pathlib import Path
import json


def extract_crispresso2_results(sample_dir: Path, sample_id: str) -> dict:
    """
    Extract key metrics from CRISPResso2 output.

    Returns dict with sample_id, total_reads, HDR%, NHEJ%, WT%, etc.
    """
    result = {
        'sample_id': sample_id,
        'crispresso2_status': 'NOT_RUN',
        'crispresso2_total_reads': 0,
        'crispresso2_reads_aligned': 0,
        'crispresso2_HDR': 0,
        'crispresso2_NHEJ': 0,
        'crispresso2_unmodified': 0,
        'crispresso2_HDR_pct': 0.0,
        'crispresso2_NHEJ_pct': 0.0,
        'crispresso2_unmodified_pct': 0.0
    }

    crispresso_dir = sample_dir / "crispresso2"
    done_file = crispresso_dir / "CRISPResso_done.txt"

    if not done_file.exists():
        return result

    # Check status
    status = done_file.read_text().strip()
    if status == "SKIPPED: Missing reference or guide sequence":
        result['crispresso2_status'] = 'SKIPPED'
        return result
    elif not status.startswith("SUCCESS"):
        result['crispresso2_status'] = 'FAILED'
        return result

    # Find CRISPResso2 output folder (named CRISPResso_on_<sample_id>)
    output_dirs = list(crispresso_dir.glob(f"CRISPResso_on_{sample_id}*"))
    if not output_dirs:
        result['crispresso2_status'] = 'MISSING_OUTPUT'
        return result

    output_dir = output_dirs[0]

    # Try to read CRISPResso2 quantification results
    # Check for multiple possible output file names
    quant_files = [
        output_dir / "CRISPResso_quantification_of_editing_frequency.txt",
        output_dir / "Quantification_of_editing_frequency.txt"
    ]

    quant_file = None
    for qf in quant_files:
        if qf.exists():
            quant_file = qf
            break

    if quant_file is None:
        # Try reading from mapping stats
        mapping_stats = output_dir / "CRISPResso_mapping_statistics.txt"
        if mapping_stats.exists():
            result['crispresso2_status'] = 'SUCCESS'
            # Parse mapping stats for basic info
            try:
                with open(mapping_stats) as f:
                    for line in f:
                        if 'READS IN INPUTS' in line:
                            parts = line.strip().split(':')
                            if len(parts) >= 2:
                                result['crispresso2_total_reads'] = int(parts[1].strip())
                        elif 'READS AFTER PREPROCESSING' in line:
                            parts = line.strip().split(':')
                            if len(parts) >= 2:
                                result['crispresso2_reads_aligned'] = int(parts[1].strip())
            except Exception:
                pass
        return result

    # Parse quantification file
    # CRISPResso2 output has one row per amplicon (Reference, HDR)
    # Columns: Amplicon, Unmodified%, Modified%, ..., Reads_aligned, Unmodified, Modified, ...
    try:
        df = pd.read_csv(quant_file, sep='\t')
        result['crispresso2_status'] = 'SUCCESS'

        # Get total reads aligned across all amplicons
        if 'Reads_aligned_all_amplicons' in df.columns:
            result['crispresso2_reads_aligned'] = int(df['Reads_aligned_all_amplicons'].iloc[0])
        elif 'Reads_aligned' in df.columns:
            result['crispresso2_reads_aligned'] = int(df['Reads_aligned'].sum())

        # Get Reference row metrics (for NHEJ - modified reads on Reference amplicon)
        ref_row = df[df['Amplicon'] == 'Reference']
        if len(ref_row) > 0:
            ref_row = ref_row.iloc[0]
            # NHEJ = Modified reads on Reference amplicon
            if 'Modified' in df.columns:
                result['crispresso2_NHEJ'] = int(ref_row['Modified'])
            # Unmodified = Unmodified reads on Reference amplicon
            if 'Unmodified' in df.columns:
                result['crispresso2_unmodified'] = int(ref_row['Unmodified'])

        # Get HDR row metrics (reads aligned to HDR amplicon = HDR)
        hdr_row = df[df['Amplicon'] == 'HDR']
        if len(hdr_row) > 0:
            hdr_row = hdr_row.iloc[0]
            # HDR = Total reads aligned to HDR amplicon
            if 'Reads_aligned' in df.columns:
                result['crispresso2_HDR'] = int(hdr_row['Reads_aligned'])

        # Calculate percentages based on total aligned reads
        total = result['crispresso2_reads_aligned']
        if total > 0:
            result['crispresso2_HDR_pct'] = 100.0 * result['crispresso2_HDR'] / total
            result['crispresso2_NHEJ_pct'] = 100.0 * result['crispresso2_NHEJ'] / total
            result['crispresso2_unmodified_pct'] = 100.0 * result['crispresso2_unmodified'] / total

    except Exception as e:
        result['crispresso2_status'] = f'PARSE_ERROR: {str(e)}'

    return result


def main():
    samples = snakemake.params.samples
    samples_dir = Path(snakemake.params.samples_dir)

    results = []
    for sample_id in samples:
        sample_dir = samples_dir / sample_id
        result = extract_crispresso2_results(sample_dir, sample_id)
        results.append(result)

    # Write summary
    summary_df = pd.DataFrame(results)
    output_path = Path(snakemake.output.tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, sep='\t', index=False)

    # Print summary statistics
    success_count = sum(1 for r in results if r['crispresso2_status'] == 'SUCCESS')
    skipped_count = sum(1 for r in results if r['crispresso2_status'] == 'SKIPPED')
    failed_count = sum(1 for r in results if 'FAILED' in r['crispresso2_status'] or 'ERROR' in r['crispresso2_status'])

    print(f"CRISPResso2 aggregation complete:")
    print(f"  SUCCESS: {success_count}/{len(samples)}")
    print(f"  SKIPPED: {skipped_count}/{len(samples)}")
    print(f"  FAILED: {failed_count}/{len(samples)}")


if __name__ == '__main__':
    main()
