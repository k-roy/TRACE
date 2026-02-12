"""
CRISPRo - CRISPR editing outcome analysis with triple-aligner consensus.

Author: Kevin R. Roy
"""

__version__ = "0.1.0"
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
