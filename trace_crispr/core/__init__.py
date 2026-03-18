"""
Core analysis modules for CRISPRo.

Author: Kevin R. Roy
"""

from .cigar import (
    CigarOperation,
    Deletion,
    Insertion,
    count_mismatches_in_region,
    get_deletions_from_cigar,
    get_insertions_from_cigar,
    parse_cigar_to_operations,
)
from .classification import (
    ClassificationResult,
    EditingOutcome,
    classify_read,
    get_hdr_signature_positions,
    summarize_classifications,
)
from .models import (
    EditingTemplate,
)
from .multi_ref_classifier import (
    AlignmentClassification,
    ClassificationSummary,
    Indel,
    MultiRefClassifier,
    classify_sample,
)
# Obsolete modules (moved to obsolete/ directory):
# - kmer: Superseded by edit_distance_hdr
# - combinatorial_hdr: Not used in current pipeline
from .scoring import (
    AlignmentScore,
    DeduplicationSignature,
    get_dedup_signature,
    score_alignment,
    select_best_alignment,
    select_best_alignment_paired,
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
    # Models
    'EditingTemplate',
    # Multi-reference alignment classifier
    'Indel',
    'AlignmentClassification',
    'ClassificationSummary',
    'MultiRefClassifier',
    'classify_sample',
]
