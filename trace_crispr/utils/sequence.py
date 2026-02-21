"""
Sequence manipulation utilities.

Provides common functions for DNA sequence operations.

Author: Kevin R. Roy
"""

import re
from typing import List, Set, Tuple


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'
    }
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


# Alias for backward compatibility
rev_comp = reverse_complement


def parse_cigar(cigar_str: str) -> List[Tuple[int, str]]:
    """Parse CIGAR string into list of (length, operation) tuples.

    CIGAR operations:
    - M: alignment match (can be match or mismatch)
    - I: insertion to reference
    - D: deletion from reference
    - N: skipped region from reference
    - S: soft clipping (sequence present but not aligned)
    - H: hard clipping (sequence not present)
    - P: padding
    - =: sequence match
    - X: sequence mismatch
    """
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    return [(int(length), op) for length, op in pattern.findall(cigar_str)]


def get_alignment_length(cigar_str: str) -> int:
    """Calculate the reference alignment length from a CIGAR string.

    Counts operations that consume reference: M, D, N, =, X
    """
    parsed = parse_cigar(cigar_str)
    ref_consuming = {'M', 'D', 'N', '=', 'X'}
    return sum(length for length, op in parsed if op in ref_consuming)


def hamming_distance(seq1: str, seq2: str) -> int:
    """Count mismatches between two equal-length sequences.

    Raises ValueError if sequences have different lengths.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be equal length: {len(seq1)} vs {len(seq2)}")
    return sum(a != b for a, b in zip(seq1.upper(), seq2.upper()))


def gc_content(seq: str) -> float:
    """Calculate GC content of a sequence (0.0 to 1.0)."""
    seq = seq.upper()
    gc = sum(1 for base in seq if base in 'GC')
    total = sum(1 for base in seq if base in 'ACGT')
    return gc / total if total > 0 else 0.0


def find_kmer_positions(sequence: str, kmer: str, allow_overlap: bool = True) -> List[int]:
    """Find all positions of a k-mer in a sequence.

    Args:
        sequence: DNA sequence to search
        kmer: K-mer to find
        allow_overlap: If True, overlapping matches are returned

    Returns:
        List of 0-based start positions
    """
    positions = []
    seq_upper = sequence.upper()
    kmer_upper = kmer.upper()
    start = 0

    while True:
        pos = seq_upper.find(kmer_upper, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1 if allow_overlap else pos + len(kmer)

    return positions


def generate_kmers_spanning_position(
    sequence: str,
    position: int,
    kmer_size: int = 12
) -> List[Tuple[str, int]]:
    """Generate all k-mers that span a specific position.

    Args:
        sequence: DNA sequence
        position: Position that must be covered by the k-mer
        kmer_size: Size of k-mers to generate

    Returns:
        List of (kmer, start_position) tuples
    """
    kmers = []
    seq_len = len(sequence)

    # K-mer must contain position, so start can range from
    # position - kmer_size + 1 to position
    start_min = max(0, position - kmer_size + 1)
    start_max = min(position, seq_len - kmer_size)

    for start in range(start_min, start_max + 1):
        kmer = sequence[start:start + kmer_size]
        if len(kmer) == kmer_size:
            kmers.append((kmer.upper(), start))

    return kmers


def generate_edit_kmers(
    reference: str,
    hdr_template: str,
    edit_positions: List[int],
    kmer_size: int = 12
) -> Tuple[Set[str], Set[str]]:
    """Generate WT and HDR k-mers spanning edit positions.

    Args:
        reference: Wild-type reference sequence
        hdr_template: HDR template sequence (same length as reference)
        edit_positions: List of positions where edits occur
        kmer_size: Size of k-mers to generate

    Returns:
        Tuple of (wt_kmers, hdr_kmers) sets including reverse complements
    """
    wt_kmers = set()
    hdr_kmers = set()

    for pos in edit_positions:
        # Generate WT k-mers from reference
        for kmer, _ in generate_kmers_spanning_position(reference, pos, kmer_size):
            wt_kmers.add(kmer)
            wt_kmers.add(reverse_complement(kmer))

        # Generate HDR k-mers from template
        for kmer, _ in generate_kmers_spanning_position(hdr_template, pos, kmer_size):
            hdr_kmers.add(kmer)
            hdr_kmers.add(reverse_complement(kmer))

    return wt_kmers, hdr_kmers


def extract_unique_kmers(
    sequence: str,
    kmer_size: int = 12,
    exclude_sequences: List[str] = None
) -> Set[str]:
    """Extract k-mers unique to a sequence (not in exclude_sequences).

    Useful for generating contaminant-specific k-mers.

    Args:
        sequence: Sequence to extract k-mers from
        kmer_size: Size of k-mers
        exclude_sequences: Sequences whose k-mers should be excluded

    Returns:
        Set of unique k-mers (including reverse complements)
    """
    # Get all k-mers from the sequence
    seq_kmers = set()
    seq_upper = sequence.upper()
    for i in range(len(seq_upper) - kmer_size + 1):
        kmer = seq_upper[i:i + kmer_size]
        seq_kmers.add(kmer)
        seq_kmers.add(reverse_complement(kmer))

    # Get k-mers from exclude sequences
    exclude_kmers = set()
    if exclude_sequences:
        for exclude_seq in exclude_sequences:
            exclude_upper = exclude_seq.upper()
            for i in range(len(exclude_upper) - kmer_size + 1):
                kmer = exclude_upper[i:i + kmer_size]
                exclude_kmers.add(kmer)
                exclude_kmers.add(reverse_complement(kmer))

    # Return k-mers unique to the sequence
    return seq_kmers - exclude_kmers


def translate_codon(codon: str) -> str:
    """Translate a DNA codon to amino acid (single letter).

    Returns '*' for stop codons, 'X' for invalid codons.
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    return codon_table.get(codon.upper(), 'X')


def find_guide_in_sequence(sequence: str, guide: str) -> Tuple[int, str]:
    """Find guide sequence in reference, searching both strands.

    Args:
        sequence: Reference sequence to search
        guide: Guide sequence (without PAM)

    Returns:
        Tuple of (position, strand) where strand is '+' or '-'

    Raises:
        ValueError: If guide not found in sequence
    """
    seq_upper = sequence.upper()
    guide_upper = guide.upper()

    # Search forward strand
    pos = seq_upper.find(guide_upper)
    if pos != -1:
        return pos, '+'

    # Search reverse strand
    guide_rc = reverse_complement(guide_upper)
    pos = seq_upper.find(guide_rc)
    if pos != -1:
        return pos, '-'

    raise ValueError(f"Guide sequence not found in reference: {guide}")
