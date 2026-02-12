"""
Utility modules for CRISPRo.

Author: Kevin R. Roy
"""

from .sequence import (
    reverse_complement,
    rev_comp,
    parse_cigar,
    get_alignment_length,
    hamming_distance,
    gc_content,
    find_kmer_positions,
    generate_kmers_spanning_position,
    generate_edit_kmers,
    extract_unique_kmers,
    translate_codon,
    find_guide_in_sequence,
)

__all__ = [
    'reverse_complement',
    'rev_comp',
    'parse_cigar',
    'get_alignment_length',
    'hamming_distance',
    'gc_content',
    'find_kmer_positions',
    'generate_kmers_spanning_position',
    'generate_edit_kmers',
    'extract_unique_kmers',
    'translate_codon',
    'find_guide_in_sequence',
]
