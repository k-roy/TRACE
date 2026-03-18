#!/usr/bin/env python3
"""
Collapse paired reads to unique R1/R2 pairs with counts.

NEW APPROACH: Store R1 and R2 separately (no merging!) to properly handle:
- Insertions at cut site (donor capture)
- Reads that can't merge due to structural variants
- Paired-end alignment benefits (insert size, proper pairing)

Uses generalized primer configuration for optimal cross-primer-pair cache sharing.
Reads primer_config.json which contains per-primer-pair settings:
- UMI presence and length on BOTH 5' (R1) and 3' (R2) ends
- Total trim lengths for both ends to reach common amplicon core

Supports three modes:
1. Paired asymmetric (e.g., R1=228bp, R2=90bp) - R1 authoritative, R2 validates
2. Paired symmetric (e.g., 2x300bp) - both span edit region
3. Single-end (R1 only) - when R2 absent or low quality

This script:
1. Extracts UMI from R1 5' end AND R2 5' end
2. Trims both reads to remove UMI + primer
3. Stores R1/R2 as separate sequences (NO MERGING)
4. Uses combined UMI for deduplication (if UMIs present)
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


def compute_hash(r1_seq, r2_seq):
    """
    Compute hash for R1/R2 pair.

    Uses combined sequence to create unique identifier.
    Uses full 32-char MD5 to avoid hash collisions with large datasets.
    For single-end mode, r2_seq will be empty string.
    """
    combined = r1_seq + '|' + r2_seq
    return hashlib.md5(combined.encode()).hexdigest()


def detect_primer_pair(sample_id, config):
    """
    Detect which primer pair this sample belongs to.
    Matches against primer pair names in config.
    """
    primer_pairs = config.get('primer_pairs', {})

    for pp_name in primer_pairs.keys():
        if pp_name in sample_id:
            return pp_name

    # Fallback: return first primer pair or None
    if primer_pairs:
        return list(primer_pairs.keys())[0]
    return None


def get_trim_params(config, primer_pair):
    """
    Get trim parameters for a primer pair from config.

    New config format has:
    - umi_5prime: {has_umi, length} - UMI on R1 (5' of merged read)
    - umi_3prime: {has_umi, length} - UMI on R2 (becomes 3' of merged read)
    - total_5prime_trim: UMI + primer offset for R1
    - total_3prime_trim: UMI + primer offset for R2

    Returns: (umi_5prime_len, umi_3prime_len, r1_trim, r2_trim)
    """
    if not primer_pair or 'primer_pairs' not in config:
        return 0, 0, 20, 20  # Default fallback

    pp_config = config['primer_pairs'].get(primer_pair, {})

    # New format: separate UMI info for 5' and 3' ends
    umi_5prime = pp_config.get('umi_5prime', {})
    umi_3prime = pp_config.get('umi_3prime', {})

    umi_5prime_len = umi_5prime.get('length', 0) if umi_5prime.get('has_umi', False) else 0
    umi_3prime_len = umi_3prime.get('length', 0) if umi_3prime.get('has_umi', False) else 0

    # Total trim = UMI + primer offset
    r1_trim = pp_config.get('total_5prime_trim', umi_5prime_len + 20)
    r2_trim = pp_config.get('total_3prime_trim', umi_3prime_len + 20)

    return umi_5prime_len, umi_3prime_len, r1_trim, r2_trim


def main():
    # Get snakemake parameters
    r1_path = snakemake.input.r1
    r2_path = snakemake.input.r2
    config_path = snakemake.input.primer_config
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

    # Check if single-end mode (from snakemake params or file detection)
    single_end_mode = getattr(snakemake.params, 'is_single_end', False)
    if not single_end_mode:
        # Also check if R2 file is empty (created as placeholder for single-end)
        try:
            with gzip.open(r2_path, 'rt') as f:
                first_line = f.readline()
                if not first_line:
                    single_end_mode = True
        except Exception:
            single_end_mode = True

    # Load primer config
    try:
        with open(config_path) as f:
            config = json.load(f)
    except Exception as e:
        config = {}

    # Check library type - Tn5 libraries don't have UMIs
    library_type = config.get('library_type', 'TruSeq')
    is_tn5 = library_type == 'Tn5'

    # Detect primer pair and get trim parameters
    primer_pair = detect_primer_pair(sample_id, config)
    umi_5prime_len, umi_3prime_len, r1_trim, r2_trim = get_trim_params(config, primer_pair)

    # For Tn5: force UMI detection off (Tn5 libraries don't have UMIs)
    if is_tn5:
        umi_5prime_len = 0
        umi_3prime_len = 0
        r1_trim = 0  # No UMI trimming for Tn5
        r2_trim = 0

    has_umi = (umi_5prime_len > 0) or (umi_3prime_len > 0)

    # Create output directory
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Count unique R1/R2 pairs
    if has_umi:
        umi_pair_counts = defaultdict(int)  # (umi, r1_seq, r2_seq) -> count
    else:
        pair_counts = defaultdict(int)  # (r1_seq, r2_seq) -> count

    total_reads = 0
    too_short = 0

    # Read FASTQ files
    if single_end_mode:
        # Single-end mode: R1 only
        with gzip.open(r1_path, 'rt') as f1:
            while True:
                r1_header = f1.readline().strip()
                if not r1_header:
                    break
                r1_seq = f1.readline().strip()
                f1.readline()
                f1.readline()

                total_reads += 1

                if len(r1_seq) <= r1_trim:
                    too_short += 1
                    continue

                # Extract UMI from R1 (if present)
                if has_umi and umi_5prime_len > 0:
                    umi = r1_seq[:umi_5prime_len].upper()
                else:
                    umi = ""

                # Trim R1
                r1_trimmed = r1_seq[r1_trim:].upper()
                r2_trimmed = ""  # Empty for single-end

                if has_umi:
                    umi_pair_counts[(umi, r1_trimmed, r2_trimmed)] += 1
                else:
                    pair_counts[(r1_trimmed, r2_trimmed)] += 1

    else:
        # Paired-end mode: R1 + R2
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

                # Extract UMI from both R1 and R2 (before trimming)
                if has_umi:
                    umi_5 = r1_seq[:umi_5prime_len].upper() if umi_5prime_len > 0 else ""
                    umi_3 = r2_seq[:umi_3prime_len].upper() if umi_3prime_len > 0 else ""
                    # Combine UMIs for deduplication
                    umi = umi_5 + "_" + umi_3 if (umi_5 and umi_3) else (umi_5 or umi_3)
                else:
                    umi = ""

                # Trim to common core (removes UMI + primer from each read)
                r1_trimmed = r1_seq[r1_trim:].upper()
                r2_trimmed = r2_seq[r2_trim:].upper()

                # Store R1/R2 pair (NO MERGING!)
                if has_umi:
                    umi_pair_counts[(umi, r1_trimmed, r2_trimmed)] += 1
                else:
                    pair_counts[(r1_trimmed, r2_trimmed)] += 1

    # Create output
    rows = []

    if has_umi:
        pair_to_umis = defaultdict(set)
        pair_total_reads = defaultdict(int)

        for (umi, r1_seq, r2_seq), count in umi_pair_counts.items():
            pair = (r1_seq, r2_seq)
            pair_to_umis[pair].add(umi)
            pair_total_reads[pair] += count

        for (r1_seq, r2_seq) in pair_to_umis:
            unique_molecules = len(pair_to_umis[(r1_seq, r2_seq)])
            total_count = pair_total_reads[(r1_seq, r2_seq)]

            seq_hash = compute_hash(r1_seq, r2_seq)
            rows.append({
                'seq_hash': seq_hash,
                'r1_sequence': r1_seq,
                'r2_sequence': r2_seq,
                'guide': guide,
                'donor': donor,
                'count': unique_molecules,
                'total_reads': total_count,
                'has_umi': True,
                'primer_pair': primer_pair
            })
    else:
        for (r1_seq, r2_seq), count in pair_counts.items():
            seq_hash = compute_hash(r1_seq, r2_seq)
            rows.append({
                'seq_hash': seq_hash,
                'r1_sequence': r1_seq,
                'r2_sequence': r2_seq,
                'guide': guide,
                'donor': donor,
                'count': count,
                'total_reads': count,
                'has_umi': False,
                'primer_pair': primer_pair
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Log
    cache_sharing = config.get('cache_sharing', {})
    num_primer_pairs = cache_sharing.get('num_primer_pairs', 1)
    mode_str = "single-end (R1 only)" if single_end_mode else "paired-end (R1 + R2)"
    library_str = f"{library_type} library"
    if is_tn5:
        library_str += " (position-based dedup will be applied post-alignment)"

    if has_umi:
        total_unique_molecules = df['count'].sum() if len(df) > 0 else 0
        pcr_dup_rate = 100 * (1 - total_unique_molecules / total_reads) if total_reads > 0 else 0
        umi_str = f"5'={umi_5prime_len}bp, 3'={umi_3prime_len}bp"
        log_msg = (
            f"Collapsed reads (PAIRED-END MODE - NO MERGING):\n"
            f"  Sample: {sample_id}\n"
            f"  Library: {library_str}\n"
            f"  Mode: {mode_str}\n"
            f"  Primer pair: {primer_pair}\n"
            f"  Config: UMI({umi_str}), R1_trim={r1_trim}bp, R2_trim={r2_trim}bp\n"
            f"  Total reads: {total_reads}\n"
            f"  Too short: {too_short}\n"
            f"  Unique R1/R2 pairs: {len(df)}\n"
            f"  Unique molecules (UMI-deduped): {total_unique_molecules}\n"
            f"  PCR dup rate: {pcr_dup_rate:.1f}%\n"
            f"  \n"
            f"  ⚠️  NO READ MERGING - R1 and R2 stored separately\n"
            f"  Cache sharing: {num_primer_pairs} primer pairs detected\n"
            f"  Expected cache efficiency: ~{num_primer_pairs}x\n"
        )
    else:
        compression = total_reads / len(df) if len(df) > 0 else 1
        dedup_note = ""
        if is_tn5:
            dedup_note = "\n  📌 Tn5 library: Position-based deduplication will be applied after alignment"
        log_msg = (
            f"Collapsed reads (PAIRED-END MODE - NO MERGING):\n"
            f"  Sample: {sample_id}\n"
            f"  Library: {library_str}\n"
            f"  Mode: {mode_str}\n"
            f"  Primer pair: {primer_pair}\n"
            f"  Config: R1_trim={r1_trim}bp, R2_trim={r2_trim}bp\n"
            f"  Total reads: {total_reads}\n"
            f"  Too short: {too_short}\n"
            f"  Unique R1/R2 pairs: {len(df)}\n"
            f"  Compression: {compression:.1f}x{dedup_note}\n"
            f"  \n"
            f"  ⚠️  NO READ MERGING - R1 and R2 stored separately\n"
            f"  Cache sharing: {num_primer_pairs} primer pairs detected\n"
            f"  Expected cache efficiency: ~{num_primer_pairs}x\n"
        )

    log_path.write_text(log_msg)


if __name__ == '__main__':
    main()
