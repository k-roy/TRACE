#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate CRISPRessoBatch settings TSV file from TRACE manifest.

Creates a TSV file with per-sample settings for CRISPRessoBatch, enabling:
- Sequence caching across samples (identical reads aligned once)
- Parallel processing
- Per-sample amplicon/guide/HDR sequences

Called by Snakemake with snakemake object providing:
- snakemake.input.manifest: TRACE manifest TSV
- snakemake.input.trimmed_r1: List of trimmed R1 FASTQs
- snakemake.input.trimmed_r2: List of trimmed R2 FASTQs
- snakemake.output.batch_settings: Output batch settings TSV
- snakemake.params.samples: List of sample IDs
- snakemake.params.default_reference: Default reference sequence
- snakemake.params.default_guide: Default guide sequence
- snakemake.params.default_hdr: Default HDR template

Author: Kevin R. Roy
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path


def revcomp(seq):
    """Return reverse complement of sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def extract_amplicon_region(reference, guide, flank_size=140):
    """
    Extract ~300bp amplicon region centered on the guide sequence.
    The guide might be in forward or reverse complement orientation.

    Returns: (amplicon_sequence, guide_in_amplicon, guide_position)
    """
    if not reference or not guide:
        return None, None, -1

    ref_upper = reference.upper()
    guide_upper = guide.upper()
    guide_rc = revcomp(guide_upper)

    # Try forward orientation
    pos = ref_upper.find(guide_upper)
    if pos >= 0:
        start = max(0, pos - flank_size)
        end = min(len(ref_upper), pos + len(guide_upper) + flank_size)
        amplicon = ref_upper[start:end]
        guide_in_amplicon = guide_upper
        guide_pos = pos - start
        return amplicon, guide_in_amplicon, guide_pos

    # Try reverse complement
    pos = ref_upper.find(guide_rc)
    if pos >= 0:
        start = max(0, pos - flank_size)
        end = min(len(ref_upper), pos + len(guide_rc) + flank_size)
        amplicon = ref_upper[start:end]
        guide_in_amplicon = guide_rc
        guide_pos = pos - start
        return amplicon, guide_in_amplicon, guide_pos

    # Guide not found - return full reference
    return ref_upper, guide_upper, -1


def extract_hdr_amplicon(reference_amplicon, hdr_template, guide):
    """
    Extract HDR amplicon by substituting donor edits into the reference amplicon.

    For simple cases, aligns donor to reference and creates hybrid.
    Returns HDR amplicon sequence.
    """
    if not hdr_template or not reference_amplicon:
        return None

    ref_upper = reference_amplicon.upper()
    hdr_upper = hdr_template.upper()
    hdr_rc = revcomp(hdr_upper)

    # Try to find best alignment of HDR template to reference
    best_offset = -1
    best_score = 0
    best_hdr = hdr_upper

    for offset in range(len(ref_upper) - len(hdr_upper) + 1):
        ref_region = ref_upper[offset:offset + len(hdr_upper)]
        matches = sum(a == b for a, b in zip(ref_region, hdr_upper))
        if matches > best_score:
            best_score = matches
            best_offset = offset
            best_hdr = hdr_upper

    # Try reverse complement
    for offset in range(len(ref_upper) - len(hdr_rc) + 1):
        ref_region = ref_upper[offset:offset + len(hdr_rc)]
        matches = sum(a == b for a, b in zip(ref_region, hdr_rc))
        if matches > best_score:
            best_score = matches
            best_offset = offset
            best_hdr = hdr_rc

    if best_offset < 0:
        return None

    # Create HDR amplicon by substituting donor region
    hdr_amplicon = (
        ref_upper[:best_offset] +
        best_hdr +
        ref_upper[best_offset + len(best_hdr):]
    )

    return hdr_amplicon


def main():
    manifest_file = Path(snakemake.input.manifest)
    output_file = Path(snakemake.output.batch_settings)

    samples = snakemake.params.samples
    samples_dir = Path(snakemake.params.samples_dir)
    default_ref = snakemake.params.default_reference
    default_guide = snakemake.params.default_guide
    default_hdr = snakemake.params.default_hdr

    print(f"Generating CRISPRessoBatch settings for {len(samples)} samples")

    # Load manifest for per-sample sequences
    manifest = pd.read_csv(manifest_file, sep='\t')
    manifest_lookup = {row['sample_id']: row for _, row in manifest.iterrows()}

    batch_rows = []
    skipped = 0

    for sample_id in samples:
        row = manifest_lookup.get(sample_id, {})

        # Get trimmed FASTQ paths
        r1_path = samples_dir / sample_id / "trimmed" / "R1_trimmed.fastq.gz"
        r2_path = samples_dir / sample_id / "trimmed" / "R2_trimmed.fastq.gz"

        if not r1_path.exists():
            print(f"SKIP {sample_id}: R1 not found at {r1_path}")
            skipped += 1
            continue

        # Check sequencing mode (default to paired-end for backwards compatibility)
        sequencing_mode = row.get('sequencing_mode', 'paired-end')
        if pd.isna(sequencing_mode) or sequencing_mode == '':
            sequencing_mode = 'paired-end'
        is_paired_end = sequencing_mode.lower() == 'paired-end'

        # Get sequences (per-sample or default)
        reference = row.get('reference', default_ref)
        guide = row.get('guide', default_guide)
        hdr_template = row.get('hdr_template', default_hdr)

        # Handle NaN values
        if pd.isna(reference) or reference == '':
            reference = default_ref
        if pd.isna(guide) or guide == '' or str(guide) == 'nan':
            guide = default_guide
        if pd.isna(hdr_template) or hdr_template == '' or str(hdr_template) == 'nan':
            hdr_template = default_hdr

        # Skip samples without guide (controls)
        if not guide or str(guide) == 'nan':
            print(f"SKIP {sample_id}: No guide sequence")
            skipped += 1
            continue

        # Extract amplicon region centered on guide
        amplicon, guide_in_amplicon, guide_pos = extract_amplicon_region(
            reference, guide, flank_size=140
        )

        if guide_pos < 0:
            print(f"WARN {sample_id}: Guide not found in reference, using full reference")

        # Extract HDR amplicon if template provided
        hdr_amplicon = None
        if hdr_template and str(hdr_template) != 'nan':
            hdr_amplicon = extract_hdr_amplicon(amplicon, hdr_template, guide)

        # Build batch row
        # Only include R2 for paired-end samples (single-end R2 files may be empty)
        batch_row = {
            'name': sample_id,
            'fastq_r1': str(r1_path.absolute()),
            'fastq_r2': str(r2_path.absolute()) if is_paired_end and r2_path.exists() else '',
            'amplicon_seq': amplicon,
            'guide_seq': guide_in_amplicon if guide_in_amplicon else guide.upper(),
        }

        if hdr_amplicon:
            batch_row['expected_hdr_amplicon_seq'] = hdr_amplicon

        batch_rows.append(batch_row)

    # Create DataFrame
    batch_df = pd.DataFrame(batch_rows)

    # Add default columns if missing
    if 'expected_hdr_amplicon_seq' not in batch_df.columns:
        batch_df['expected_hdr_amplicon_seq'] = ''

    print(f"\nGenerated settings for {len(batch_df)} samples (skipped {skipped})")

    # Write output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    batch_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved batch settings to: {output_file}")


if __name__ == '__main__':
    main()
