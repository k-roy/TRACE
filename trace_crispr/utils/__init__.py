"""
Utility modules for CRISPRo.

Author: Kevin R. Roy
"""

from .multi_ref_builder import (
    build_emx1_multi_reference,
    build_multi_reference_fasta,
    calculate_cut_site,
    load_multi_reference_fasta,
)
from .sequence import (
    extract_unique_kmers,
    find_guide_in_sequence,
    find_kmer_positions,
    gc_content,
    generate_edit_kmers,
    generate_kmers_spanning_position,
    get_alignment_length,
    hamming_distance,
    parse_cigar,
    rev_comp,
    reverse_complement,
    translate_codon,
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
    # Multi-reference builder
    'build_multi_reference_fasta',
    'load_multi_reference_fasta',
    'calculate_cut_site',
    'build_emx1_multi_reference',
]
