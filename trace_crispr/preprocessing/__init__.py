"""
Preprocessing modules for read filtering, trimming, and merging.

Author: Kevin R. Roy
"""

from .detection import (
    detect_library_type,
    detect_merge_need,
    detect_crispresso_mode,
    run_auto_detection,
    LibraryDetectionResult,
    MergeDetectionResult,
    CRISPRessoModeResult,
    AutoDetectionResults,
)

from .trimming import (
    trim_adapters,
    run_bbduk_trim,
    TrimmingResult,
    TN5_ADAPTERS,
    TRUSEQ_ADAPTERS,
)

from .contamination import (
    filter_contamination_fastq,
    create_contamination_filter,
    ContaminationFilterResult,
)

__all__ = [
    # Detection
    'detect_library_type',
    'detect_merge_need',
    'detect_crispresso_mode',
    'run_auto_detection',
    'LibraryDetectionResult',
    'MergeDetectionResult',
    'CRISPRessoModeResult',
    'AutoDetectionResults',
    # Trimming
    'trim_adapters',
    'run_bbduk_trim',
    'TrimmingResult',
    'TN5_ADAPTERS',
    'TRUSEQ_ADAPTERS',
    # Contamination
    'filter_contamination_fastq',
    'create_contamination_filter',
    'ContaminationFilterResult',
]
