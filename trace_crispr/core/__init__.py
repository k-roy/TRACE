"""
Core analysis modules for CRISPRo.

Author: Kevin R. Roy
"""

from .cigar import (
    CigarOperation,
    Deletion,
    Insertion,
    parse_cigar_to_operations,
    get_deletions_from_cigar,
    get_insertions_from_cigar,
    count_mismatches_in_region,
)

from .scoring import (
    AlignmentScore,
    score_alignment,
    select_best_alignment,
    select_best_alignment_paired,
    DeduplicationSignature,
    get_dedup_signature,
)

from .classification import (
    EditingOutcome,
    ClassificationResult,
    classify_read,
    summarize_classifications,
    get_hdr_signature_positions,
)

from .kmer import (
    KmerClassification,
    KmerClassifier,
    KmerResults,
    classify_fastq_kmer,
)

from .models import (
    EditingTemplate,
)

from .multi_ref_classifier import (
    Indel,
    AlignmentClassification,
    ClassificationSummary,
    MultiRefClassifier,
    classify_sample,
)

__all__ = [
    # CIGAR
    'CigarOperation',
    'Deletion',
    'Insertion',
    'parse_cigar_to_operations',
    'get_deletions_from_cigar',
    'get_insertions_from_cigar',
    'count_mismatches_in_region',
    # Scoring
    'AlignmentScore',
    'score_alignment',
    'select_best_alignment',
    'select_best_alignment_paired',
    'DeduplicationSignature',
    'get_dedup_signature',
    # Classification
    'EditingOutcome',
    'ClassificationResult',
    'classify_read',
    'summarize_classifications',
    'get_hdr_signature_positions',
    # K-mer
    'KmerClassification',
    'KmerClassifier',
    'KmerResults',
    'classify_fastq_kmer',
    # Models
    'EditingTemplate',
    # Multi-reference alignment classifier
    'Indel',
    'AlignmentClassification',
    'ClassificationSummary',
    'MultiRefClassifier',
    'classify_sample',
]
