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

from .bc1_consensus import (
    BRIDGE_AMBIGUOUS,
    BRIDGE_EXACT,
    BRIDGE_HD1,
    BRIDGE_NONE,
    BarcodeCall,
    Bc1Call,
    BridgeIndex,
    BridgeResult,
    Extraction,
    TrackCall,
    TrackConfig,
    call_sample,
    caller_confidence,
    confident_barcodes,
    extract_barcode,
    extract_track,
    is_amplicon_read,
    is_read2_hop,
    top_barcode,
    tracks_from_config,
)
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
from .reconcile import (
    CONTAMINATION,
    LOW_CONF,
    MULTI_EDIT,
    NO_EDIT,
    SINGLE_EDIT,
    SITE_FULL,
    SITE_NOCOV,
    SITE_PARTIAL,
    SITE_WT,
    CellReconciliation,
    EditSite,
    SiteCall,
    build_all_edit_sites,
    build_edit_sites,
    classify_site,
    designed_site_vaf,
    load_oligo_windows,
    reconcile_cell,
    reconciliation_to_row,
    window_edit_vaf,
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
    # bc1_consensus (Stage 2-ii)
    "TrackConfig",
    "Extraction",
    "TrackCall",
    "extract_track",
    "extract_barcode",
    "is_amplicon_read",
    "is_read2_hop",
    "BridgeIndex",
    "BridgeResult",
    "BRIDGE_EXACT",
    "BRIDGE_HD1",
    "BRIDGE_AMBIGUOUS",
    "BRIDGE_NONE",
    "Bc1Call",
    "BarcodeCall",
    "top_barcode",
    "caller_confidence",
    "call_sample",
    "confident_barcodes",
    "tracks_from_config",
    # reconcile (Stage 3)
    "NO_EDIT",
    "SINGLE_EDIT",
    "MULTI_EDIT",
    "CONTAMINATION",
    "LOW_CONF",
    "SITE_FULL",
    "SITE_PARTIAL",
    "SITE_WT",
    "SITE_NOCOV",
    "EditSite",
    "SiteCall",
    "CellReconciliation",
    "classify_site",
    "reconcile_cell",
    "designed_site_vaf",
    "window_edit_vaf",
    "build_all_edit_sites",
    "build_edit_sites",
    "load_oligo_windows",
    "reconciliation_to_row",
]
