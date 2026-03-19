#!/usr/bin/env python3
"""
Collapse paired reads to unique merged sequences with counts.

Key optimizations:
1. UMI handling: KR2478 samples have 6bp UMI at 5' of R1
2. Primer trimming: After UMI removal, trim primer sequences so that
   sequences from both primer pairs become identical for caching

Trimming strategy:
- KR2478 (UMI): 6bp UMI + 20bp primer = 26bp from R1 5', 20bp from R2 5'
- KR2476 (no UMI): 20bp primer from both R1 and R2 5' ends

After trimming, sequences from both primer pairs are identical,
enabling 2x cache efficiency.

Output: TSV with columns: seq_hash, merged_seq, guide, donor, count, ...
"""

import gzip
import hashlib
from pathlib import Path
from collections import defaultdict

import pandas as pd


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def merge_reads(r1_seq, r2_seq, min_overlap=20):
    """
    Merge paired-end reads by finding overlap.
    Returns merged sequence or None if no good overlap found.
    """
    r2_rc = reverse_complement(r2_seq)

    best_overlap = 0
    best_merged = None

    for overlap in range(min_overlap, min(len(r1_seq), len(r2_rc)) + 1):
        r1_end = r1_seq[-overlap:]
        r2_start = r2_rc[:overlap]

        mismatches = sum(1 for a, b in zip(r1_end, r2_start) if a != b)

        if mismatches <= overlap * 0.1:
            if overlap > best_overlap:
                best_overlap = overlap
                best_merged = r1_seq + r2_rc[overlap:]

    return best_merged


def compute_hash(seq):
    """Compute hash for sequence (guide/donor independent for Level 1)."""
    return hashlib.md5(seq.encode()).hexdigest()[:16]


def detect_primer_config(sample_id):
    """
    Detect primer configuration from sample_id.
    Returns: (has_umi, umi_length, r1_trim, r2_trim)

    KR2478: 6bp UMI + 20bp primer on R1, 20bp primer on R2
    KR2476: 20bp primer on both R1 and R2 (no UMI)
    """
    if "KR2478" in sample_id:
        # UMI + primer
        return True, 6, 26, 20  # (umi, umi_len, r1_trim, r2_trim)
    elif "KR2476" in sample_id:
        # Primer only
        return False, 0, 20, 20
    else:
        # Default: assume no UMI, 20bp primers
        return False, 0, 20, 20


def main():
    # Get snakemake parameters
    r1_path = snakemake.input.r1
    r2_path = snakemake.input.r2
    output_path = snakemake.output.collapsed
    log_path = Path(str(snakemake.log))

    guide = str(snakemake.params.guide) if snakemake.params.guide else ""
    donor = str(snakemake.params.hdr_template) if snakemake.params.hdr_template else ""

    # Clean up NaN values
    if guide == "nan":
        guide = ""
    if donor == "nan":
        donor = ""

    # Extract sample_id from path to detect primer config
    sample_id = Path(r1_path).parent.parent.name
    has_umi, umi_length, r1_trim, r2_trim = detect_primer_config(sample_id)

    # Create output directory
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Count unique sequences
    if has_umi:
        # Track (umi, seq) -> count for UMI deduplication
        umi_seq_counts = defaultdict(int)
    else:
        # Track seq -> count
        seq_counts = defaultdict(int)

    total_reads = 0
    merged_reads = 0
    too_short = 0

    with gzip.open(r1_path, 'rt') as f1, gzip.open(r2_path, 'rt') as f2:
        while True:
            # Read R1 record
            r1_header = f1.readline().strip()
            if not r1_header:
                break
            r1_seq = f1.readline().strip()
            f1.readline()  # +
            f1.readline()  # qual

            # Read R2 record
            f2.readline()  # header
            r2_seq = f2.readline().strip()
            f2.readline()  # +
            f2.readline()  # qual

            total_reads += 1

            # Skip if too short after trimming
            if len(r1_seq) <= r1_trim or len(r2_seq) <= r2_trim:
                too_short += 1
                continue

            # Extract UMI and trim primers
            if has_umi:
                umi = r1_seq[:umi_length].upper()
                r1_trimmed = r1_seq[r1_trim:]  # Remove UMI + primer
            else:
                umi = ""
                r1_trimmed = r1_seq[r1_trim:]  # Remove primer only

            r2_trimmed = r2_seq[r2_trim:]  # Remove primer from R2

            # Try to merge trimmed reads
            merged = merge_reads(r1_trimmed, r2_trimmed)
            if merged:
                merged_reads += 1
                seq = merged.upper()
            else:
                # Use concatenated if can't merge
                seq = (r1_trimmed + "NNNNNNNNNN" + reverse_complement(r2_trimmed)).upper()

            # Count by appropriate key
            if has_umi:
                umi_seq_counts[(umi, seq)] += 1
            else:
                seq_counts[seq] += 1

    # Create output dataframe
    rows = []

    if has_umi:
        # For UMI samples: collapse by sequence, count unique UMIs
        seq_to_umis = defaultdict(set)
        seq_total_reads = defaultdict(int)

        for (umi, seq), count in umi_seq_counts.items():
            seq_to_umis[seq].add(umi)
            seq_total_reads[seq] += count

        for seq in seq_to_umis:
            unique_molecules = len(seq_to_umis[seq])
            total_count = seq_total_reads[seq]

            seq_hash = compute_hash(seq)
            rows.append({
                'seq_hash': seq_hash,
                'merged_seq': seq,
                'guide': guide,
                'donor': donor,
                'count': unique_molecules,  # Unique molecules (PCR-deduplicated)
                'total_reads': total_count,
                'has_umi': True,
                'primer_pair': 'KR2478'
            })
    else:
        for seq, count in seq_counts.items():
            seq_hash = compute_hash(seq)
            rows.append({
                'seq_hash': seq_hash,
                'merged_seq': seq,
                'guide': guide,
                'donor': donor,
                'count': count,
                'total_reads': count,
                'has_umi': False,
                'primer_pair': 'KR2476'
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Write log
    if has_umi:
        total_unique_molecules = df['count'].sum() if len(df) > 0 else 0
        total_unique_seqs = len(df)
        pcr_dup_rate = 100 * (1 - total_unique_molecules / total_reads) if total_reads > 0 else 0
        log_msg = (
            f"Collapsed reads (KR2478 with UMI):\n"
            f"  Sample: {sample_id}\n"
            f"  UMI length: {umi_length}bp\n"
            f"  R1 trim: {r1_trim}bp, R2 trim: {r2_trim}bp\n"
            f"  Total reads: {total_reads}\n"
            f"  Too short after trim: {too_short}\n"
            f"  Merged reads: {merged_reads} ({100*merged_reads/(total_reads-too_short):.1f}%)\n"
            f"  Unique sequences: {total_unique_seqs}\n"
            f"  Unique molecules (UMI families): {total_unique_molecules}\n"
            f"  PCR duplication rate: {pcr_dup_rate:.1f}%\n"
            f"  Guide: {guide[:20]}{'...' if len(guide) > 20 else ''}\n"
            f"  Donor: {donor[:20]}{'...' if len(donor) > 20 else ''}\n"
        )
    else:
        total_unique_seqs = len(df)
        compression = total_reads / total_unique_seqs if total_unique_seqs > 0 else 1
        log_msg = (
            f"Collapsed reads (KR2476 no UMI):\n"
            f"  Sample: {sample_id}\n"
            f"  R1 trim: {r1_trim}bp, R2 trim: {r2_trim}bp\n"
            f"  Total reads: {total_reads}\n"
            f"  Too short after trim: {too_short}\n"
            f"  Merged reads: {merged_reads} ({100*merged_reads/(total_reads-too_short):.1f}%)\n"
            f"  Unique sequences: {total_unique_seqs}\n"
            f"  Compression ratio: {compression:.1f}x\n"
            f"  Guide: {guide[:20]}{'...' if len(guide) > 20 else ''}\n"
            f"  Donor: {donor[:20]}{'...' if len(donor) > 20 else ''}\n"
        )

    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
