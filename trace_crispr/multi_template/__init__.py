"""
Multi-template analysis module for TRACE.

Provides support for analyzing samples with multiple possible HDR templates,
such as barcoded CRISPR editing experiments.

Author: Kevin R. Roy
"""

from .batch_processor import BatchMultiTemplateProcessor

__all__ = ['BatchMultiTemplateProcessor']
