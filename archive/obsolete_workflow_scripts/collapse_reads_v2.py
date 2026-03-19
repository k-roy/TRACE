#!/usr/bin/env python3
"""
Collapse paired reads to unique merged sequences with counts.

Uses auto-detected trim configuration for optimal cross-primer-pair cache sharing.
Reads trim_config.json to determine:
- UMI length (if any)
- 5' trim length for R1 and R2

This ensures sequences from different primer pairs can share the same cache.
"""

import gzip
import json
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
    """Merge paired-end reads by finding overlap."""
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
    """Compute hash for sequence."""
    return hashlib.md5(seq.encode()).hexdigest()[:16]


def get_trim_config(config_path, sample_id):
    """
    Load trim configuration for this sample.
    Returns (umi_length, r1_trim, r2_trim)
    """
    try:
        with open(config_path) as f:
            config = json.load(f)
    except Exception as e:
        # Fallback to default if config not found
        return 0, 20, 20

    # Detect primer pair from sample_id
    if 'KR2478' in sample_id:
        pp = 'KR2478'
    elif 'KR2476' in sample_id:
        pp = 'KR2476'
    else:
        # Use first available config
        pp = list(config.keys())[0] if config else None

    if pp and pp in config:
        pp_config = config[pp]
        umi_len = pp_config.get('umi_length', 0)
        r1_trim = pp_config.get('r1_5prime_trim', umi_len + 20)
        r2_trim = pp_config.get('r2_5prime_trim', 20)
        return umi_len, r1_trim, r2_trim

    # Default fallback
    return 0, 20, 20


def main():
    # Get snakemake parameters
    r1_path = snakemake.input.r1
    r2_path = snakemake.input.r2
    trim_config_path = snakemake.input.trim_config
    output_path = snakemake.output.collapsed
    log_path = Path(str(snakemake.log))

    sample_id = snakemake.params.sample_id
    guide = str(snakemake.params.guide) if snakemake.params.guide else ""
    donor = str(snakemake.params.hdr_template) if snakemake.params.hdr_template else ""

    # Clean up NaN values
    if guide == "nan":
        guide = ""
    if donor == "nan":
        donor = ""

    # Get trim config
    umi_length, r1_trim, r2_trim = get_trim_config(trim_config_path, sample_id)
    has_umi = umi_length > 0

    # Create output directory
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Count unique sequences
    if has_umi:
        umi_seq_counts = defaultdict(int)  # (umi, seq) -> count
    else:
        seq_counts = defaultdict(int)  # seq -> count

    total_reads = 0
    merged_reads = 0
    too_short = 0

    with gzip.open(r1_path, 'rt') as f1, gzip.open(r2_path, 'rt') as f2:
        while True:
            r1_header = f1.readline().strip()
            if not r1_header:
                break
            r1_seq = f1.readline().strip()
            f1.readline()
            f1.readline()

            f2.readline()
            r2_seq = f2.readline().strip()
            f2.readline()
            f2.readline()

            total_reads += 1

            if len(r1_seq) <= r1_trim or len(r2_seq) <= r2_trim:
                too_short += 1
                continue

            # Extract UMI and trim
            if has_umi:
                umi = r1_seq[:umi_length].upper()
                r1_trimmed = r1_seq[r1_trim:]
            else:
                umi = ""
                r1_trimmed = r1_seq[r1_trim:]

            r2_trimmed = r2_seq[r2_trim:]

            # Merge
            merged = merge_reads(r1_trimmed, r2_trimmed)
            if merged:
                merged_reads += 1
                seq = merged.upper()
            else:
                seq = (r1_trimmed + "NNNNNNNNNN" + reverse_complement(r2_trimmed)).upper()

            if has_umi:
                umi_seq_counts[(umi, seq)] += 1
            else:
                seq_counts[seq] += 1

    # Create output
    rows = []

    if has_umi:
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
                'count': unique_molecules,
                'total_reads': total_count,
                'has_umi': True
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
                'has_umi': False
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Log
    if has_umi:
        total_unique_molecules = df['count'].sum() if len(df) > 0 else 0
        pcr_dup_rate = 100 * (1 - total_unique_molecules / total_reads) if total_reads > 0 else 0
        log_msg = (
            f"Collapsed reads (with UMI):\n"
            f"  Sample: {sample_id}\n"
            f"  Config: UMI={umi_length}bp, R1_trim={r1_trim}bp, R2_trim={r2_trim}bp\n"
            f"  Total reads: {total_reads}\n"
            f"  Too short: {too_short}\n"
            f"  Merged: {merged_reads}\n"
            f"  Unique sequences: {len(df)}\n"
            f"  Unique molecules: {total_unique_molecules}\n"
            f"  PCR dup rate: {pcr_dup_rate:.1f}%\n"
        )
    else:
        compression = total_reads / len(df) if len(df) > 0 else 1
        log_msg = (
            f"Collapsed reads (no UMI):\n"
            f"  Sample: {sample_id}\n"
            f"  Config: R1_trim={r1_trim}bp, R2_trim={r2_trim}bp\n"
            f"  Total reads: {total_reads}\n"
            f"  Too short: {too_short}\n"
            f"  Merged: {merged_reads}\n"
            f"  Unique sequences: {len(df)}\n"
            f"  Compression: {compression:.1f}x\n"
        )

    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
