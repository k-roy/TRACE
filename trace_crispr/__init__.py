"""
TRACE - Triple-aligner Read Analysis for CRISPR Editing.

Author: Kevin R. Roy
"""

__version__ = "0.2.0"
__author__ = "Kevin R. Roy"

from .config import (
    LocusConfig,
    SampleConfig,
    PipelineConfig,
    NucleaseType,
)
from .core.classification import EditingOutcome, ClassificationResult

__all__ = [
    "LocusConfig",
    "SampleConfig",
    "PipelineConfig",
    "NucleaseType",
    "EditingOutcome",
    "ClassificationResult",
    "__version__",
]
