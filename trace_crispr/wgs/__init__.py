"""
trace_crispr.wgs
================

WGS editing-outcome backend for TRACE.

Extends TRACE beyond amplicon data to genome-wide MAGESTIC outcome mapping:
GATK-called variants in a donor window are classified as designed / synthetic-
error / spontaneous (:mod:`~trace_crispr.wgs.donor_outcome`), donor-template and
chimeric reads are filtered for SV / plasmid-integration detection, bc1 barcodes
are called from the BAM to bridge each cell to its designed donor, and genome-
wide variants are attributed back to the oligo library.

The WGS algorithms are ported from Shengdi Li's Perl/R pipeline
(``WGS_MAGESTIC_QTL_outcome_mapping`` + ``WGS_MAGESTIC_QTL_bc_calling``;
https://github.com/shli-embl/MAGESTIC-SCORE) by Kevin R. Roy, with bit-for-bit
parity tests against Shengdi's reference outputs. Each ported module carries a
provenance docstring crediting the original.
"""

from .donor_outcome import (
    DESIGNED,
    SPONTANEOUS,
    SYNTHETIC_ERROR,
    DonorSpec,
    SampleOutcome,
    bucket_variant,
    levenshtein,
    precall,
    precall_sample,
)

__all__ = [
    "DESIGNED",
    "SYNTHETIC_ERROR",
    "SPONTANEOUS",
    "DonorSpec",
    "SampleOutcome",
    "bucket_variant",
    "levenshtein",
    "precall",
    "precall_sample",
]
