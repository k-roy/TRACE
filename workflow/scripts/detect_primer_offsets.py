#!/usr/bin/env python3
"""
Generalized primer offset detection for multi-primer-pair amplicon sequencing.

For merged paired-end reads, the structure is:
  5' end: [R1_UMI][R1_primer][amplicon_core][R2_primer_RC][R2_UMI_RC] :3' end

This script:
1. Detects library type (TruSeq vs Tn5) based on read start position clustering
2. For TruSeq: Detects UMI presence/length on BOTH 5' and 3' ends
3. For Tn5: Skips UMI detection (Tn5 libraries don't have UMIs)
4. Compares primer pairs to find common amplicon core
5. Calculates trim lengths for both ends to produce identical core sequences
6. Outputs configuration for the collapse step

This enables cache sharing across all primer pairs.
"""

import gzip
import json
from pathlib import Path
from collections import Counter
from difflib import SequenceMatcher

import pandas as pd


# Adapter signatures for library type detection
TN5_ADAPTER = "CTGTCTCTTATACACATCT"
TRUSEQ_ADAPTER = "AGATCGGAAGAGC"


def detect_library_type(r1_reads, r2_reads, reference=None, n_reads=1000):
    """
    Auto-detect TruSeq vs Tn5 based on adapter signatures and position clustering.

    Primary detection: Adapter signatures in reads
    - Tn5: CTGTCTCTTATACACATCT (Nextera/Tn5 transposase)
    - TruSeq: AGATCGGAAGAGC

    Secondary detection: Position clustering (if reference provided)
    - TruSeq: Reads cluster at fixed start positions (fixed primers)
    - Tn5: Reads have scattered start positions (random tagmentation)

    Returns: dict with library_type, confidence, and detection_reason
    """
    if not r1_reads:
        return {
            'library_type': 'TruSeq',
            'confidence': 'low',
            'detection_reason': 'no reads found, defaulting to TruSeq'
        }

    # Count adapter signatures
    tn5_count = sum(1 for r in r1_reads[:n_reads] if TN5_ADAPTER in r.upper())
    truseq_count = sum(1 for r in r1_reads[:n_reads] if TRUSEQ_ADAPTER in r.upper())
    reads_checked = min(len(r1_reads), n_reads)

    tn5_fraction = tn5_count / reads_checked if reads_checked > 0 else 0
    truseq_fraction = truseq_count / reads_checked if reads_checked > 0 else 0

    # Primary detection: adapter signatures
    if tn5_fraction > 0.1 and tn5_fraction > truseq_fraction * 2:
        return {
            'library_type': 'Tn5',
            'confidence': 'high' if tn5_fraction > 0.3 else 'medium',
            'detection_reason': f'Nextera adapter in {tn5_fraction:.1%} of reads',
            'tn5_adapter_fraction': tn5_fraction,
            'truseq_adapter_fraction': truseq_fraction
        }
    elif truseq_fraction > 0.1 and truseq_fraction > tn5_fraction * 2:
        return {
            'library_type': 'TruSeq',
            'confidence': 'high' if truseq_fraction > 0.3 else 'medium',
            'detection_reason': f'TruSeq adapter in {truseq_fraction:.1%} of reads',
            'tn5_adapter_fraction': tn5_fraction,
            'truseq_adapter_fraction': truseq_fraction
        }

    # Secondary detection: position clustering (if reference provided)
    if reference:
        position_clustering, n_unique = _calculate_position_clustering(r1_reads, reference)

        if position_clustering >= 0.7:
            return {
                'library_type': 'TruSeq',
                'confidence': 'high' if position_clustering >= 0.85 else 'medium',
                'detection_reason': f'{position_clustering:.0%} of reads cluster at fixed start position',
                'position_clustering': position_clustering,
                'n_unique_positions': n_unique
            }
        elif n_unique >= 20 and position_clustering < 0.5:
            return {
                'library_type': 'Tn5',
                'confidence': 'high' if position_clustering < 0.3 else 'medium',
                'detection_reason': f'scattered start positions ({n_unique} unique positions)',
                'position_clustering': position_clustering,
                'n_unique_positions': n_unique
            }

    # Default to TruSeq if ambiguous
    return {
        'library_type': 'TruSeq',
        'confidence': 'low',
        'detection_reason': 'ambiguous, defaulting to TruSeq'
    }


def _calculate_position_clustering(reads, reference, kmer_size=15):
    """Calculate position clustering for library type detection."""
    ref_upper = reference.upper()
    ref_rc = reverse_complement(ref_upper)

    # Build k-mer index
    ref_kmers = {}
    for i in range(len(ref_upper) - kmer_size + 1):
        kmer = ref_upper[i:i + kmer_size]
        if kmer not in ref_kmers:
            ref_kmers[kmer] = i
    for i in range(len(ref_rc) - kmer_size + 1):
        kmer = ref_rc[i:i + kmer_size]
        if kmer not in ref_kmers:
            ref_kmers[kmer] = len(ref_upper) - i - kmer_size

    # Find start positions
    start_positions = []
    for read in reads[:500]:
        if len(read) < kmer_size + 20:
            continue
        for offset in range(0, min(60, len(read) - kmer_size)):
            read_kmer = read[offset:offset + kmer_size].upper()
            if read_kmer in ref_kmers:
                ref_pos = ref_kmers[read_kmer]
                read_start = ref_pos - offset
                if 0 <= read_start < len(ref_upper):
                    start_positions.append(read_start)
                    break

    if not start_positions:
        return 0.0, 0

    position_counts = Counter(start_positions)
    n_unique = len(position_counts)

    if position_counts:
        dominant_pos, dominant_count = position_counts.most_common(1)[0]
        clustered_count = sum(
            count for pos, count in position_counts.items()
            if abs(pos - dominant_pos) <= 5
        )
        position_clustering = clustered_count / len(start_positions)
        return position_clustering, n_unique

    return 0.0, 0


def sample_reads(fastq_path, n_reads=1000):
    """Sample first n reads from a FASTQ file."""
    reads = []
    try:
        with gzip.open(fastq_path, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    reads.append(line.strip().upper())
                    if len(reads) >= n_reads:
                        break
    except Exception as e:
        print(f"Warning: Could not read {fastq_path}: {e}")
    return reads


def detect_umi_5prime(reads, max_umi_len=12, min_umi_len=4, n_sample=500):
    """
    Detect UMI at 5' end by looking for position where diversity drops.

    Strategy: Find the position where per-base diversity drops sharply.
    - UMI positions: all 4 bases present (diversity ~4)
    - Primer positions: dominated by single base (diversity ~1)

    Returns: (has_umi, umi_length)
    """
    if not reads or len(reads) < 100:
        return False, 0

    # Sample subset for analysis
    import random
    if len(reads) > n_sample:
        reads = random.sample(reads, n_sample)

    # Calculate per-position diversity (number of unique bases / max possible)
    # Look for sharp drop from ~4 unique bases to ~1
    from collections import Counter

    position_diversity = []
    for pos in range(max_umi_len + 5):
        chars = [r[pos] for r in reads if len(r) > pos]
        if len(chars) < 50:
            break
        counts = Counter(chars)
        # Diversity = number of bases seen at least 5% of the time
        significant_bases = sum(1 for c, cnt in counts.items() if cnt >= len(chars) * 0.05)
        position_diversity.append(significant_bases)

    if len(position_diversity) < min_umi_len + 2:
        return False, 0

    # Find UMI boundary: position where diversity drops from >=3 to <=2
    umi_length = 0
    for pos in range(min_umi_len - 1, min(max_umi_len, len(position_diversity) - 1)):
        # Check if this position is diverse and next position is conserved
        current_div = position_diversity[pos]
        next_div = position_diversity[pos + 1] if pos + 1 < len(position_diversity) else 1

        # UMI positions have 3-4 bases; primer positions have 1-2 bases
        if current_div >= 3 and next_div <= 2:
            umi_length = pos + 1  # pos is 0-indexed, length is pos+1
            break

    if umi_length >= min_umi_len:
        # Verify: check that positions before boundary are diverse
        diverse_count = sum(1 for d in position_diversity[:umi_length] if d >= 3)
        if diverse_count >= umi_length * 0.7:  # At least 70% of UMI positions are diverse
            return True, umi_length

    return False, 0


def detect_umi_3prime(reads, max_umi_len=12, min_umi_len=4, n_sample=500):
    """
    Detect UMI at 3' end by looking for position where diversity drops.

    For 3' end, we look from the end of the read inward.
    - UMI positions (at 3' end): all 4 bases present
    - Primer/amplicon positions: dominated by single base

    Returns: (has_umi, umi_length)
    """
    if not reads or len(reads) < 100:
        return False, 0

    # Sample subset for analysis
    import random
    if len(reads) > n_sample:
        reads = random.sample(reads, n_sample)

    # Calculate per-position diversity from 3' end
    from collections import Counter

    position_diversity = []
    for offset in range(max_umi_len + 5):
        # offset=0 is last char, offset=1 is second-to-last, etc.
        chars = [r[-(offset + 1)] for r in reads if len(r) > offset]
        if len(chars) < 50:
            break
        counts = Counter(chars)
        significant_bases = sum(1 for c, cnt in counts.items() if cnt >= len(chars) * 0.05)
        position_diversity.append(significant_bases)

    if len(position_diversity) < min_umi_len + 2:
        return False, 0

    # Find UMI boundary: position where diversity drops from >=3 to <=2
    # (going inward from 3' end)
    umi_length = 0
    for pos in range(min_umi_len - 1, min(max_umi_len, len(position_diversity) - 1)):
        current_div = position_diversity[pos]
        next_div = position_diversity[pos + 1] if pos + 1 < len(position_diversity) else 1

        if current_div >= 3 and next_div <= 2:
            umi_length = pos + 1
            break

    if umi_length >= min_umi_len:
        diverse_count = sum(1 for d in position_diversity[:umi_length] if d >= 3)
        if diverse_count >= umi_length * 0.7:
            return True, umi_length

    return False, 0


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def detect_primer_from_reference(reads, reference, umi_len=0, n_sample=100):
    """
    Detect primer length by finding where reads align to reference.

    Reads may be in reverse complement orientation relative to reference.
    The primer is the portion at the read start that doesn't match reference.

    Returns: primer_length (bp)
    """
    if not reads or not reference:
        return 0

    import random
    if len(reads) > n_sample:
        reads = random.sample(reads, n_sample)

    ref_upper = reference.upper()
    ref_rc = reverse_complement(ref_upper)

    primer_lengths = []
    for read in reads:
        # Strip UMI if present
        read_stripped = read[umi_len:].upper() if umi_len > 0 else read.upper()

        # Try to find a 30bp segment from the read in the reference
        # Start from position 15 onwards (skip potential primer region)
        for start in range(15, min(50, len(read_stripped) - 30)):
            segment = read_stripped[start:start + 30]

            # Check forward orientation
            pos = ref_upper.find(segment)
            if pos >= 0:
                # Read matches at this position - the primer is everything before 'start'
                primer_lengths.append(start)
                break

            # Check reverse complement
            pos_rc = ref_rc.find(segment)
            if pos_rc >= 0:
                primer_lengths.append(start)
                break

    if primer_lengths:
        # Return median primer length
        primer_lengths.sort()
        return primer_lengths[len(primer_lengths) // 2]

    return 0


def find_common_core(read_sets, umi_info_5prime, umi_info_3prime, min_core_len=100, reference=None):
    """
    Find the longest common substring present in reads from all primer pairs.

    First strips UMIs from both ends, then finds common core.
    If reference is provided, uses it to detect primer offsets when primer pairs have identical sequences.
    Returns the core and the offsets (primer lengths) for each primer pair.
    """
    # Strip UMIs from reads
    stripped_reads = {}
    for pp_name, reads in read_sets.items():
        umi_5 = umi_info_5prime[pp_name]['umi_length']
        umi_3 = umi_info_3prime[pp_name]['umi_length']

        stripped = []
        for r in reads:
            if len(r) > umi_5 + umi_3 + min_core_len:
                if umi_3 > 0:
                    stripped.append(r[umi_5:-umi_3])
                else:
                    stripped.append(r[umi_5:])
        stripped_reads[pp_name] = stripped

    if len(stripped_reads) < 2:
        # Single primer pair - use the most common read as reference
        pp_name = list(stripped_reads.keys())[0]
        counter = Counter(stripped_reads[pp_name])
        top_read = counter.most_common(1)[0][0] if counter else ''
        return {
            'common_core': top_read,
            'offsets_5prime': {pp_name: 0},
            'offsets_3prime': {pp_name: 0}
        }

    # Get top reads from each primer pair
    consensus_reads = {}
    for pp_name, reads in stripped_reads.items():
        counter = Counter(reads)
        top_reads = [r for r, _ in counter.most_common(50)]
        consensus_reads[pp_name] = top_reads

    pp_names = list(consensus_reads.keys())
    best_core = ""
    best_offsets_5 = {}
    best_offsets_3 = {}

    # Compare reads from first primer pair against others to find common core
    ref_pp = pp_names[0]
    other_pp = pp_names[1]

    for ref_read in consensus_reads[ref_pp][:10]:
        for other_read in consensus_reads[other_pp][:10]:
            # Find longest common substring
            matcher = SequenceMatcher(None, ref_read, other_read)
            match = matcher.find_longest_match(0, len(ref_read), 0, len(other_read))

            if match.size > len(best_core) and match.size >= min_core_len:
                best_core = ref_read[match.a:match.a + match.size]

                # Calculate offsets for all primer pairs
                offsets_5 = {}
                offsets_3 = {}

                for pp in pp_names:
                    for read in consensus_reads[pp][:20]:
                        pos = read.find(best_core)
                        if pos >= 0:
                            offsets_5[pp] = pos  # 5' offset = primer length
                            offsets_3[pp] = len(read) - pos - len(best_core)  # 3' offset
                            break
                    if pp not in offsets_5:
                        offsets_5[pp] = 0
                        offsets_3[pp] = 0

                best_offsets_5 = offsets_5
                best_offsets_3 = offsets_3

    # If offsets are all 0 and we have a reference, try reference-based detection
    if reference and all(v == 0 for v in best_offsets_5.values()):
        for pp_name, reads in stripped_reads.items():
            umi_len = umi_info_5prime[pp_name]['umi_length']
            primer_len = detect_primer_from_reference(reads, reference, umi_len=0)  # Already UMI-stripped
            if primer_len > 0:
                best_offsets_5[pp_name] = primer_len
                # Assume symmetric primers for 3' end if not detected
                if best_offsets_3.get(pp_name, 0) == 0:
                    best_offsets_3[pp_name] = primer_len

    return {
        'common_core': best_core,
        'offsets_5prime': best_offsets_5,
        'offsets_3prime': best_offsets_3
    }


def analyze_primer_pairs(manifest_df, n_samples_per_pp=5, n_reads=1000):
    """
    Analyze all primer pairs in the manifest.

    Returns configuration dict with:
    - Library type (TruSeq or Tn5)
    - Per primer pair: UMI info (both ends), trim lengths
    - Common core sequence
    - Cache sharing statistics
    """
    # Group samples by primer pair
    if 'primer_pair' not in manifest_df.columns:
        def infer_primer_pair(sample_id):
            patterns = ['KR2478', 'KR2476', 'PP1', 'PP2', 'primer1', 'primer2']
            for p in patterns:
                if p in sample_id:
                    return p
            return 'default'

        manifest_df = manifest_df.copy()
        manifest_df['primer_pair'] = manifest_df['sample_id'].apply(infer_primer_pair)

    primer_pairs = manifest_df['primer_pair'].unique()

    # Get reference sequence from manifest (for primer detection)
    reference = None
    if 'reference' in manifest_df.columns:
        for _, row in manifest_df.iterrows():
            ref = row.get('reference', '')
            if ref and str(ref) != 'nan' and len(str(ref)) > 100:
                reference = str(ref).upper()
                break

    # Sample reads from each primer pair
    read_sets = {}
    umi_info_5prime = {}
    umi_info_3prime = {}
    all_r1_for_libtype = []
    all_r2_for_libtype = []

    for pp in primer_pairs:
        pp_samples = manifest_df[manifest_df['primer_pair'] == pp]

        # Sample R1 reads for 5' UMI detection
        all_r1_reads = []
        all_r2_reads = []
        for _, sample in pp_samples.head(n_samples_per_pp).iterrows():
            r1_reads = sample_reads(sample['r1_path'], n_reads=n_reads)
            all_r1_reads.extend(r1_reads)

            # Also read R2 for 3' UMI detection (R2's 5' = merged read's 3')
            if 'r2_path' in sample and pd.notna(sample['r2_path']):
                r2_reads = sample_reads(sample['r2_path'], n_reads=n_reads)
                all_r2_reads.extend(r2_reads)

        read_sets[pp] = all_r1_reads  # Use R1 for core finding
        all_r1_for_libtype.extend(all_r1_reads[:500])
        all_r2_for_libtype.extend(all_r2_reads[:500])

        # We'll detect UMI after library type is known
        umi_info_5prime[pp] = {'has_umi': False, 'umi_length': 0}
        umi_info_3prime[pp] = {'has_umi': False, 'umi_length': 0}

    # Detect library type FIRST (before UMI detection)
    library_info = detect_library_type(all_r1_for_libtype, all_r2_for_libtype, reference)
    library_type = library_info['library_type']

    # For TruSeq: detect UMIs
    # For Tn5: skip UMI detection (Tn5 libraries don't have UMIs)
    if library_type == 'TruSeq':
        for pp in primer_pairs:
            pp_samples = manifest_df[manifest_df['primer_pair'] == pp]

            all_r1_reads = []
            all_r2_reads = []
            for _, sample in pp_samples.head(n_samples_per_pp).iterrows():
                r1_reads = sample_reads(sample['r1_path'], n_reads=n_reads)
                all_r1_reads.extend(r1_reads)

                if 'r2_path' in sample and pd.notna(sample['r2_path']):
                    r2_reads = sample_reads(sample['r2_path'], n_reads=n_reads)
                    all_r2_reads.extend(r2_reads)

            # Detect UMI on 5' end (from R1)
            has_umi_5, umi_len_5 = detect_umi_5prime(all_r1_reads)
            umi_info_5prime[pp] = {'has_umi': has_umi_5, 'umi_length': umi_len_5}

            # Detect UMI on 3' end (from R2's 5' end)
            if all_r2_reads:
                has_umi_3, umi_len_3 = detect_umi_5prime(all_r2_reads)
            else:
                has_umi_3, umi_len_3 = False, 0
            umi_info_3prime[pp] = {'has_umi': has_umi_3, 'umi_length': umi_len_3}
    # For Tn5: UMI info stays as False/0 (already set above)

    # Find common core across all primer pairs
    core_result = find_common_core(read_sets, umi_info_5prime, umi_info_3prime, reference=reference)

    # Build final configuration
    config = {
        'library_type': library_type,
        'library_detection': library_info,
        'primer_pairs': {},
        'common_core': {
            'sequence': core_result['common_core'][:50] + '...' if len(core_result['common_core']) > 50 else core_result['common_core'],
            'length': len(core_result['common_core'])
        },
        'cache_sharing': {
            'enabled': len(primer_pairs) > 1,
            'num_primer_pairs': len(primer_pairs)
        }
    }

    for pp in primer_pairs:
        umi_len_5 = umi_info_5prime[pp]['umi_length']
        umi_len_3 = umi_info_3prime[pp]['umi_length']
        primer_offset_5 = core_result['offsets_5prime'].get(pp, 0)
        primer_offset_3 = core_result['offsets_3prime'].get(pp, 0)

        # Total trim = UMI + primer offset
        total_5prime_trim = umi_len_5 + primer_offset_5
        total_3prime_trim = umi_len_3 + primer_offset_3

        config['primer_pairs'][pp] = {
            'umi_5prime': {
                'has_umi': umi_info_5prime[pp]['has_umi'],
                'length': umi_len_5
            },
            'umi_3prime': {
                'has_umi': umi_info_3prime[pp]['has_umi'],
                'length': umi_len_3
            },
            'primer_offset_5prime': primer_offset_5,
            'primer_offset_3prime': primer_offset_3,
            'total_5prime_trim': total_5prime_trim,
            'total_3prime_trim': total_3prime_trim,
            'sample_count': len(manifest_df[manifest_df['primer_pair'] == pp])
        }

    return config


def main():
    """Snakemake entry point."""
    manifest_path = snakemake.input.manifest
    output_path = snakemake.output.trim_config
    log_path = Path(str(snakemake.log))

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Load manifest
    manifest = pd.read_csv(manifest_path, sep='\t')

    # Analyze primer pairs
    config = analyze_primer_pairs(manifest)

    # Save config
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)

    # Generate detailed log
    lib_info = config['library_detection']
    log_lines = [
        "=" * 70,
        "Primer Pair Offset Detection Results",
        "=" * 70,
        "",
        f"Library type: {config['library_type']} ({lib_info['detection_reason']})",
        f"  Confidence: {lib_info['confidence']}",
        "",
        f"Number of primer pairs: {config['cache_sharing']['num_primer_pairs']}",
        f"Common core length: {config['common_core']['length']} bp",
        "",
    ]

    if config['library_type'] == 'Tn5':
        log_lines.extend([
            "NOTE: Tn5 library detected - UMI detection skipped",
            "  Tn5 libraries use position-based deduplication (post-alignment)",
            "",
        ])

    log_lines.append("Per-primer-pair configuration:")

    for pp, pp_config in config['primer_pairs'].items():
        umi_5_str = f"Yes ({pp_config['umi_5prime']['length']}bp)" if pp_config['umi_5prime']['has_umi'] else "No"
        umi_3_str = f"Yes ({pp_config['umi_3prime']['length']}bp)" if pp_config['umi_3prime']['has_umi'] else "No"

        log_lines.extend([
            f"  {pp}:",
            f"    Samples: {pp_config['sample_count']}",
            f"    5' UMI: {umi_5_str}",
            f"    3' UMI: {umi_3_str}",
            f"    5' primer offset: {pp_config['primer_offset_5prime']} bp",
            f"    3' primer offset: {pp_config['primer_offset_3prime']} bp",
            f"    Total 5' trim: {pp_config['total_5prime_trim']} bp (UMI + primer)",
            f"    Total 3' trim: {pp_config['total_3prime_trim']} bp (UMI + primer)",
            ""
        ])

    if config['cache_sharing']['enabled']:
        log_lines.extend([
            "Cache sharing: ENABLED",
            f"  All {config['cache_sharing']['num_primer_pairs']} primer pairs will share alignment cache",
            f"  After trimming, all reads will have identical core sequences",
            f"  Expected efficiency gain: ~{config['cache_sharing']['num_primer_pairs']}x",
        ])

    log_path.write_text('\n'.join(log_lines))


if __name__ == '__main__':
    main()
