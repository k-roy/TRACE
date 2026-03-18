#!/usr/bin/env python3
"""
Classify unique sequences using edit distance HDR detection.

Each unique (sequence, guide, donor) combination is classified once.
Results are stored in a cache for later expansion to per-sample results.
"""

import sys
from pathlib import Path

# Add TRACE to path
TRACE_ROOT = Path("/oak/stanford/groups/larsms/Users/kevinroy/software/trace")
sys.path.insert(0, str(TRACE_ROOT))

import pysam
import pandas as pd

from trace_crispr.config import parse_sequence_input
from trace_crispr.core.edit_distance_hdr import (
    build_donor_signature, classify_read_edit_distance, align_donor_to_reference
)


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def find_cut_site(ref_seq, guide_seq):
    """Find cut site position from guide location in reference."""
    if not guide_seq:
        return len(ref_seq) // 2

    ref_upper = ref_seq.upper()
    guide_upper = guide_seq.upper()

    guide_pos = ref_upper.find(guide_upper)
    if guide_pos == -1:
        guide_rc = reverse_complement(guide_upper)
        guide_pos = ref_upper.find(guide_rc)

    if guide_pos == -1:
        return len(ref_seq) // 2

    return guide_pos + len(guide_seq) - 3


def main():
    bam_path = snakemake.input.bam
    unique_seqs_path = snakemake.input.unique_seqs
    output_path = snakemake.output.classifications
    log_path = Path(str(snakemake.log))
    edit_region_size = snakemake.params.edit_region_size
    ref_param = snakemake.params.reference

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Load unique sequences
    seqs_df = pd.read_csv(unique_seqs_path, sep='\t', compression='gzip')

    # Create lookup by seq_hash
    seq_info = {}
    for _, row in seqs_df.iterrows():
        seq_info[row['seq_hash']] = {
            'merged_seq': row['merged_seq'],
            'guide': row['guide'] if pd.notna(row['guide']) else "",
            'donor': row['donor'] if pd.notna(row['donor']) else "",
            'total_count': row['total_count']
        }

    # Parse reference
    ref_seq = parse_sequence_input(ref_param) if ref_param else None

    # Cache donor signatures by (guide, donor) to avoid recomputation
    signature_cache = {}

    def get_donor_signature(guide, donor):
        key = (guide, donor)
        if key not in signature_cache:
            if not donor:
                signature_cache[key] = ({}, None, {})
            else:
                cut_site = find_cut_site(ref_seq, guide)
                edit_region = (cut_site - edit_region_size, cut_site + edit_region_size)

                try:
                    donor_sig = build_donor_signature(ref_seq, donor, cut_site, edit_region=edit_region)
                    signature_cache[key] = (donor_sig, cut_site, edit_region)
                except Exception as e:
                    signature_cache[key] = ({}, cut_site, None)

        return signature_cache[key]

    # Process alignments
    results = []
    classified = 0
    unaligned = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            seq_hash = read.query_name

            if seq_hash not in seq_info:
                continue

            info = seq_info[seq_hash]
            guide = info['guide']
            donor = info['donor']

            # Get cached donor signature
            donor_signature, cut_site, edit_region = get_donor_signature(guide, donor)

            if read.is_unmapped or read.query_sequence is None:
                unaligned += 1
                results.append({
                    'seq_hash': seq_hash,
                    'outcome': 'UNALIGNED',
                    'n_donor_snvs': 0,
                    'n_non_donor': 0,
                    'donor_fraction': 0.0,
                    'donor_snvs_detected': '',
                    'total_count': info['total_count']
                })
                continue

            # Classify
            if not cut_site:
                cut_site = len(ref_seq) // 2

            clf = classify_read_edit_distance(
                read_seq=read.query_sequence,
                ref_seq=ref_seq,
                donor_signature=donor_signature,
                ref_start=read.reference_start,
                cigar_ops=read.cigartuples,
                cut_site=cut_site
            )

            classified += 1
            results.append({
                'seq_hash': seq_hash,
                'outcome': clf.outcome,
                'n_donor_snvs': clf.n_donor_encoded,
                'n_non_donor': clf.n_non_donor,
                'donor_fraction': clf.donor_fraction,
                'donor_snvs_detected': ','.join(str(p) for p in sorted(clf.donor_snvs_detected)) if clf.donor_snvs_detected else '',
                'total_count': info['total_count']
            })

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Write log
    log_msg = (
        f"Classified unique sequences:\n"
        f"  Total unique sequences: {len(seqs_df)}\n"
        f"  Classified: {classified}\n"
        f"  Unaligned: {unaligned}\n"
        f"  \n"
        f"  Unique (guide, donor) combinations: {len(signature_cache)}\n"
        f"  \n"
        f"Outcome distribution:\n"
    )

    if len(results_df) > 0:
        outcome_counts = results_df['outcome'].value_counts()
        for outcome, count in outcome_counts.items():
            # Weight by total_count for actual read counts
            weighted = results_df[results_df['outcome'] == outcome]['total_count'].sum()
            log_msg += f"  {outcome}: {count} unique seqs, {weighted} total reads\n"

    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
