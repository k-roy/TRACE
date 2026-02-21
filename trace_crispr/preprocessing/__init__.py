"""
Preprocessing modules for read filtering, trimming, and merging.

Author: Kevin R. Roy
"""

from .contamination import (
    ContaminationFilterResult,
    create_contamination_filter,
    filter_contamination_fastq,
)
from .detection import (
    AutoDetectionResults,
    CRISPRessoModeResult,
    LibraryDetectionResult,
    MergeDetectionResult,
    PreprocessingMode,
    detect_crispresso_mode,
    detect_library_type,
    detect_merge_need,
    run_auto_detection,
)
from .preprocess import (
    PreprocessingResult,
    run_collapse,
    run_preprocessing,
    run_read_merge,
    run_umi_dedup,
)
from .trimming import (
    TN5_ADAPTERS,
    TRUSEQ_ADAPTERS,
    TrimmingResult,
    run_bbduk_trim,
    trim_adapters,
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
