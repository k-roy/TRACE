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
    PreprocessingMode,
)

from .preprocess import (
    run_preprocessing,
    run_umi_dedup,
    run_read_merge,
    run_collapse,
    PreprocessingResult,
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
    'PreprocessingMode',
    # Preprocessing orchestration
    'run_preprocessing',
    'run_umi_dedup',
    'run_read_merge',
    'run_collapse',
    'PreprocessingResult',
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
