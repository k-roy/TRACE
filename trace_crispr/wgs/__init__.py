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
from .sv_plasmid import (
    CleanResult,
    SamRecord,
    build_integrated_plasmid_contig,
    classify_donor_reads,
    clean_donor_reads,
    decode_duplicate,
    detect_plasmid_integration,
    integrated_plasmid_junction,
    load_blacklist,
    one_read_inside_large_deletion,
    ref_span_md,
    rev_comp,
    summarize_plasmid_integration,
    within_read_large_deletions,
)

__all__ = [
    # donor_outcome (Stage 1)
    "DESIGNED",
    "SYNTHETIC_ERROR",
    "SPONTANEOUS",
    "DonorSpec",
    "SampleOutcome",
    "bucket_variant",
    "levenshtein",
    "precall",
    "precall_sample",
    # sv_plasmid (Stage 2-i)
    "CleanResult",
    "SamRecord",
    "classify_donor_reads",
    "clean_donor_reads",
    "decode_duplicate",
    "ref_span_md",
    "load_blacklist",
    "detect_plasmid_integration",
    "summarize_plasmid_integration",
    "one_read_inside_large_deletion",
    "within_read_large_deletions",
    "build_integrated_plasmid_contig",
    "integrated_plasmid_junction",
    "rev_comp",
]
