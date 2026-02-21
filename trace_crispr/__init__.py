"""
TRACE - Triple-aligner Read Analysis for CRISPR Editing.

Author: Kevin R. Roy
"""

__version__ = "0.4.0"
__author__ = "Kevin R. Roy"

from .config import (
    LocusConfig,
    NucleaseType,
    PipelineConfig,
    SampleConfig,
)
from .core.classification import ClassificationResult, EditingOutcome

__all__ = [
    "LocusConfig",
    "SampleConfig",
    "PipelineConfig",
    "NucleaseType",
    "EditingOutcome",
    "ClassificationResult",
    "__version__",
]
