#!/usr/bin/env python3
"""
Extract paired-end alignment results from BAM into a lookup table.

NEW APPROACH: Extract BOTH R1 and R2 alignments for each seq_hash.
Stores pairing metadata (insert size, is_proper_pair) for validation.

Supports three modes:
- Paired-end: Both R1 and R2 alignments present
- Single-end: Only R1 alignment present

This allows Level 2 classification to use R1 (primary) with R2 (validation).
"""

from pathlib import Path
from collections import defaultdict

import pysam
import pandas as pd


def cigar_to_string(cigar_tuples):
    """Convert CIGAR tuples to string format."""
    if not cigar_tuples:
        return ""
    op_map = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
    return ''.join(f"{length}{op_map.get(op, '?')}" for op, length in cigar_tuples)


def extract_alignment_info(read):
    """Extract alignment information from a pysam read."""
    if read.is_unmapped or read.query_sequence is None:
        return {
            'is_mapped': False,
            'ref_start': -1,
            'ref_end': -1,
            'cigar_string': '',
            'query_sequence': read.query_sequence or '',
            'is_reverse': False,
            'mapping_quality': 0,
        }
    else:
        return {
            'is_mapped': True,
            'ref_start': read.reference_start,
            'ref_end': read.reference_end,
            'cigar_string': cigar_to_string(read.cigartuples),
            'query_sequence': read.query_sequence,
            'is_reverse': read.is_reverse,
            'mapping_quality': read.mapping_quality,
        }


def main():
    bam_path = snakemake.input.bam
    output_path = snakemake.output.cache
    log_path = Path(str(snakemake.log))

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Group alignments by seq_hash
    # Each seq_hash should have R1 and optionally R2
    pairs = defaultdict(lambda: {'r1': None, 'r2': None, 'pairing': {}})

    total_reads = 0
    r1_count = 0
    r2_count = 0
    single_end_pairs = 0
    paired_end_pairs = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            total_reads += 1
            seq_hash = read.query_name

            # Extract alignment info
            aln_info = extract_alignment_info(read)

            # Determine if R1 or R2
            if read.is_read1:
                r1_count += 1
                pairs[seq_hash]['r1'] = aln_info
            elif read.is_read2:
                r2_count += 1
                pairs[seq_hash]['r2'] = aln_info
            else:
                # Single-end mode: treat as R1
                r1_count += 1
                pairs[seq_hash]['r1'] = aln_info

            # Extract pairing metadata (only if properly paired)
            if read.is_paired and read.is_proper_pair:
                pairs[seq_hash]['pairing'] = {
                    'is_proper_pair': True,
                    'insert_size': read.template_length,
                }

    # Create flattened output
    rows = []
    for seq_hash, pair_data in pairs.items():
        r1 = pair_data['r1']
        r2 = pair_data['r2']
        pairing = pair_data['pairing']

        # Check mode
        if r2 is None:
            single_end_pairs += 1
            mode = 'single_end'
        else:
            paired_end_pairs += 1
            mode = 'paired_end'

        # Flattened row structure
        row = {
            'seq_hash': seq_hash,
            'mode': mode,
        }

        # R1 alignment info
        if r1:
            row.update({
                'r1_is_mapped': r1['is_mapped'],
                'r1_ref_start': r1['ref_start'],
                'r1_ref_end': r1['ref_end'],
                'r1_cigar': r1['cigar_string'],
                'r1_sequence': r1['query_sequence'],
                'r1_is_reverse': r1['is_reverse'],
                'r1_mapping_quality': r1['mapping_quality'],
            })
        else:
            row.update({
                'r1_is_mapped': False,
                'r1_ref_start': -1,
                'r1_ref_end': -1,
                'r1_cigar': '',
                'r1_sequence': '',
                'r1_is_reverse': False,
                'r1_mapping_quality': 0,
            })

        # R2 alignment info
        if r2:
            row.update({
                'r2_is_mapped': r2['is_mapped'],
                'r2_ref_start': r2['ref_start'],
                'r2_ref_end': r2['ref_end'],
                'r2_cigar': r2['cigar_string'],
                'r2_sequence': r2['query_sequence'],
                'r2_is_reverse': r2['is_reverse'],
                'r2_mapping_quality': r2['mapping_quality'],
            })
        else:
            row.update({
                'r2_is_mapped': False,
                'r2_ref_start': -1,
                'r2_ref_end': -1,
                'r2_cigar': '',
                'r2_sequence': '',
                'r2_is_reverse': False,
                'r2_mapping_quality': 0,
            })

        # Pairing metadata
        row.update({
            'is_proper_pair': pairing.get('is_proper_pair', False),
            'insert_size': pairing.get('insert_size', 0),
        })

        rows.append(row)

    # Save cache
    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Calculate mapping statistics
    r1_mapped = sum(1 for row in rows if row['r1_is_mapped'])
    r2_mapped = sum(1 for row in rows if row['r2_is_mapped'])
    proper_pairs = sum(1 for row in rows if row['is_proper_pair'])

    log_msg = (
        f"Paired-end alignment cache extracted:\n"
        f"  Total reads processed: {total_reads}\n"
        f"  Total unique pairs: {len(pairs)}\n"
        f"  \n"
        f"  Mode distribution:\n"
        f"    Single-end pairs: {single_end_pairs} ({100*single_end_pairs/len(pairs):.1f}%)\n"
        f"    Paired-end pairs: {paired_end_pairs} ({100*paired_end_pairs/len(pairs):.1f}%)\n"
        f"  \n"
        f"  R1 alignments:\n"
        f"    Total: {r1_count}\n"
        f"    Mapped: {r1_mapped} ({100*r1_mapped/len(pairs):.1f}%)\n"
        f"  \n"
        f"  R2 alignments:\n"
        f"    Total: {r2_count}\n"
        f"    Mapped: {r2_mapped} ({100*r2_mapped/max(r2_count,1):.1f}% of R2)\n"
        f"  \n"
        f"  Pairing:\n"
        f"    Proper pairs: {proper_pairs} ({100*proper_pairs/len(pairs):.1f}%)\n"
        f"  \n"
        f"  ⚠️  PAIRED-END CACHE: Both R1 and R2 alignments stored\n"
        f"     Classification will use R1 (primary) + R2 (validation)\n"
    )
    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
