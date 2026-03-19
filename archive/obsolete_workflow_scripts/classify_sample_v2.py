#!/usr/bin/env python3
"""
Classify reads from a single sample using TRACE with edit distance HDR detection.

This version uses the edit distance approach which:
- Handles any number of SNVs without combinatorial explosion
- Filters to core edits near the cut site (ignores flanking sequences)
- Tracks per-SNV integration for conversion tract analysis
- Handles insertions/deletions via CIGAR parsing

Called by Snakemake with snakemake object providing:
- snakemake.input: Input files (BAM)
- snakemake.output: Output files (classification TSV, HDR detail TSV)
- snakemake.params: Parameters (reference, guide, hdr_template, etc.)
- snakemake.log: Log file path
- snakemake.wildcards: Sample ID
"""

import sys
from pathlib import Path
import pandas as pd
import pysam

# Add TRACE package to path (for script: directive execution)
TRACE_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(TRACE_ROOT))

# Import TRACE modules
try:
    from trace_crispr.config import LocusConfig, NucleaseType, parse_sequence_input
    from trace_crispr.core.classification import classify_read, get_hdr_signature_positions
    from trace_crispr.core.edit_distance_hdr import (
        build_donor_signature, classify_read_edit_distance,
        analyze_donor_integration, align_donor_to_reference
    )
except ModuleNotFoundError as e:
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(f"ERROR: Failed to import trace_crispr: {e}\n")
    sys.exit(1)


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def find_cut_site(ref_seq, guide_seq):
    """Find cut site position from guide location in reference."""
    ref_upper = ref_seq.upper()
    guide_upper = guide_seq.upper()

    # Try forward orientation
    guide_pos = ref_upper.find(guide_upper)
    if guide_pos == -1:
        # Try reverse complement
        guide_rc = reverse_complement(guide_upper)
        guide_pos = ref_upper.find(guide_rc)

    if guide_pos == -1:
        # Guide not found, use middle of reference
        return len(ref_seq) // 2

    # Cut site is 3bp upstream of PAM (end of guide) for Cas9
    return guide_pos + len(guide_seq) - 3


def main():
    # Get parameters
    sample_id = snakemake.wildcards.sample
    reference = snakemake.params.reference
    guide = snakemake.params.guide
    hdr_template = snakemake.params.hdr_template

    # Edit distance specific parameters
    edit_region_size = snakemake.params.get('edit_region_size', 30)  # bp around cut site

    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Check required inputs
    if not guide or pd.isna(guide):
        Path(snakemake.output.tsv).write_text("read_name\toutcome\tn_donor_snvs\tn_non_donor\tdonor_snvs_detected\n")
        if hasattr(snakemake.output, 'hdr_detail'):
            Path(snakemake.output.hdr_detail).write_text("position\tdistance_to_cut\tref_base\tdonor_base\tfrequency\tcount\n")
        log_path.write_text(f"Skipped {sample_id} - no guide sequence\n")
        return

    has_hdr_template = hdr_template and not pd.isna(hdr_template) and str(hdr_template).strip()

    try:
        # Parse sequences
        ref_seq = parse_sequence_input(reference)
        guide_seq = guide.upper()

        # Find cut site
        cut_site = find_cut_site(ref_seq, guide_seq)

        # Build donor signature if HDR template available
        donor_signature = {}
        donor_info = {}

        if has_hdr_template:
            donor_seq = parse_sequence_input(hdr_template)

            # Get donor alignment info
            offset, aligned_donor, is_rc, raw_edits = align_donor_to_reference(ref_seq, donor_seq)

            # Build signature with core edit region only
            edit_region = (cut_site - edit_region_size, cut_site + edit_region_size)
            donor_signature = build_donor_signature(ref_seq, donor_seq, cut_site, edit_region=edit_region)

            donor_info = {
                'offset': offset,
                'is_rc': is_rc,
                'total_edits': len(raw_edits),
                'core_edits': len(donor_signature)
            }

        # Get BAM files - try in order of preference
        aligner_paths = [
            ('bwa', snakemake.input.bwa),
            ('bbmap', snakemake.input.bbmap),
            ('minimap2', snakemake.input.minimap2)
        ]

        results = []
        aligner_used = None

        for aligner_name, bam_path in aligner_paths:
            if not Path(bam_path).exists():
                continue

            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    for read in bam:
                        if read.is_unmapped or read.query_sequence is None:
                            continue

                        # Classify using edit distance approach
                        clf = classify_read_edit_distance(
                            read_seq=read.query_sequence,
                            ref_seq=ref_seq,
                            donor_signature=donor_signature,
                            ref_start=read.reference_start,
                            cigar_ops=read.cigartuples,
                            cut_site=cut_site
                        )

                        # Format detected SNV positions
                        snv_positions = ','.join(str(p) for p in sorted(clf.donor_snvs_detected)) if clf.donor_snvs_detected else ''

                        results.append({
                            'read_name': read.query_name,
                            'outcome': clf.outcome,
                            'n_donor_snvs': clf.n_donor_encoded,
                            'n_non_donor': clf.n_non_donor,
                            'donor_fraction': clf.donor_fraction,
                            'donor_snvs_detected': snv_positions
                        })

                aligner_used = aligner_name
                break

            except Exception as e:
                log_path.write_text(f"Warning: {aligner_name} failed ({e}), trying next...\n")
                continue

        if aligner_used is None:
            raise ValueError("All aligners failed")

        # Deduplicate paired-end reads
        df = pd.DataFrame(results)
        original_count = len(df)

        if len(df) > 0:
            # Keep read with highest donor_fraction per read name
            df = df.sort_values('donor_fraction', ascending=False)
            df = df.drop_duplicates(subset='read_name', keep='first')

        # Save classification results
        df.to_csv(snakemake.output.tsv, sep='\t', index=False)

        # Generate per-SNV frequency analysis if HDR template was used
        if has_hdr_template and hasattr(snakemake.output, 'hdr_detail') and len(donor_signature) > 0:
            # Reconstruct classifications for analysis
            from trace_crispr.core.edit_distance_hdr import EditDistanceClassification

            classifications = []
            for _, row in df.iterrows():
                detected = [int(p) for p in row['donor_snvs_detected'].split(',') if p]
                missing = [p for p in donor_signature.keys() if p not in detected]

                clf = EditDistanceClassification(
                    outcome=row['outcome'],
                    edits=[],
                    n_donor_encoded=row['n_donor_snvs'],
                    n_non_donor=row['n_non_donor'],
                    n_total_edits=row['n_donor_snvs'] + row['n_non_donor'],
                    donor_fraction=row['donor_fraction'],
                    donor_snvs_detected=detected,
                    donor_snvs_missing=missing
                )
                classifications.append(clf)

            # Analyze conversion tract
            analysis = analyze_donor_integration(classifications, donor_signature, cut_site)

            # Save per-position frequency
            hdr_detail_rows = []
            for pos, info in analysis['per_position_frequency'].items():
                hdr_detail_rows.append({
                    'position': pos,
                    'distance_to_cut': info['distance_to_cut'],
                    'ref_base': info['ref_base'],
                    'donor_base': info['donor_base'],
                    'frequency': info['frequency'],
                    'count': info['count']
                })

            hdr_detail_df = pd.DataFrame(hdr_detail_rows)
            hdr_detail_df.to_csv(snakemake.output.hdr_detail, sep='\t', index=False)
        elif hasattr(snakemake.output, 'hdr_detail'):
            # Write empty hdr_detail file for samples without HDR template
            Path(snakemake.output.hdr_detail).write_text("position\tdistance_to_cut\tref_base\tdonor_base\tfrequency\tcount\n")

        # Write log
        outcome_counts = df['outcome'].value_counts().to_dict() if len(df) > 0 else {}
        log_msg = (
            f"Classified {len(df)} reads for {sample_id} using {aligner_used}\n"
            f"  Cut site: {cut_site}\n"
            f"  Donor signature: {len(donor_signature)} SNVs in core region\n"
            f"  Outcomes: {outcome_counts}\n"
            f"  Deduplicated from {original_count} alignments\n"
        )
        if donor_info:
            log_msg += f"  Donor offset: {donor_info['offset']}, RC: {donor_info['is_rc']}\n"
            log_msg += f"  Total donor edits: {donor_info['total_edits']}, Core: {donor_info['core_edits']}\n"

        log_path.write_text(log_msg)

    except Exception as e:
        import traceback
        log_path.write_text(f"ERROR classifying {sample_id}: {e}\n{traceback.format_exc()}\n")
        Path(snakemake.output.tsv).write_text("read_name\toutcome\tn_donor_snvs\tn_non_donor\tdonor_snvs_detected\n")
        if hasattr(snakemake.output, 'hdr_detail'):
            Path(snakemake.output.hdr_detail).write_text("position\tdistance_to_cut\tref_base\tdonor_base\tfrequency\tcount\n")
        raise


if __name__ == '__main__':
    main()