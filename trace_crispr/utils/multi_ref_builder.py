"""
Multi-Reference FASTA Builder.

Build a FASTA reference containing WT sequence plus all HDR template variants.
Used for alignment-based classification of CRISPR editing outcomes.

Author: Kevin R. Roy
Date: 2026-02-19
"""

from pathlib import Path
from typing import Dict, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def build_multi_reference_fasta(
    wt_reference: str,
    hdr_templates: Dict[str, str],
    output_path: Path,
    guide_seq: Optional[str] = None,
) -> Dict[str, int]:
    """
    Create a multi-sequence FASTA containing WT and all HDR variants.

    Args:
        wt_reference: Wild-type reference sequence (full amplicon)
        hdr_templates: Dict mapping barcode_id -> HDR donor template sequence
            These can be short donors (just the edited region) - they will be
            aligned to WT and used to build full-length edited amplicons.
        output_path: Path to write the FASTA file
        guide_seq: Optional guide sequence for cut site calculation

    Returns:
        Dict mapping sequence_name -> cut_site_position (0-indexed)

    The FASTA will contain:
        >WT
        <wt_reference sequence>
        >HDR_<barcode_id>
        <full-length edited amplicon with barcode>
        ...
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cut_sites = {}
    wt_upper = wt_reference.upper()

    with open(output_path, 'w') as f:
        # Write WT reference
        f.write(">WT\n")
        for i in range(0, len(wt_reference), 80):
            f.write(f"{wt_reference[i:i+80].upper()}\n")

        # Calculate WT cut site if guide provided
        wt_cut_site = None
        if guide_seq:
            wt_cut_site = calculate_cut_site(wt_reference, guide_seq)
            cut_sites['WT'] = wt_cut_site
        else:
            cut_sites['WT'] = None

        # Build full-length HDR templates
        for barcode_id, template_seq in sorted(hdr_templates.items()):
            seq_name = f"HDR_{barcode_id}"
            template_upper = template_seq.upper()

            # Check if template is already full-length
            if len(template_seq) >= len(wt_reference) - 20:
                # Template is already full-length, use as-is
                full_hdr = template_upper
            else:
                # Template is a short donor - build full-length edited amplicon
                full_hdr = build_full_length_hdr(wt_upper, template_upper)

            f.write(f">{seq_name}\n")
            for i in range(0, len(full_hdr), 80):
                f.write(f"{full_hdr[i:i+80]}\n")

            # HDR templates have the guide/PAM region modified - no cut site needed
            cut_sites[seq_name] = None

    n_sequences = 1 + len(hdr_templates)
    logger.info(f"Built multi-reference FASTA with {n_sequences} sequences")
    logger.info(f"  WT length: {len(wt_reference)} bp")
    logger.info(f"  HDR templates: {len(hdr_templates)}")
    logger.info(f"  Output: {output_path}")

    return cut_sites


def calculate_cut_site(sequence: str, guide_seq: str) -> Optional[int]:
    """
    Calculate Cas9 cut site position in a sequence.

    The cut site is 3bp upstream of the PAM (i.e., between positions 17-18
    of the guide for a 20bp guide targeting the sense strand).

    Args:
        sequence: Target sequence
        guide_seq: Guide sequence (20bp)

    Returns:
        Cut site position (0-indexed), or None if guide not found
    """
    sequence_upper = sequence.upper()
    guide_upper = guide_seq.upper()

    # Find guide position
    guide_pos = sequence_upper.find(guide_upper)
    if guide_pos == -1:
        # Try reverse complement
        guide_rc = reverse_complement(guide_upper)
        guide_pos = sequence_upper.find(guide_rc)
        if guide_pos == -1:
            logger.warning(f"Guide not found in sequence")
            return None
        # For antisense guide, cut site is at guide_pos + 3
        return guide_pos + 3

    # For sense guide, cut site is at guide_pos + len(guide) - 3
    # (Cas9 cuts 3bp upstream of PAM, which is 3bp after guide end)
    return guide_pos + len(guide_upper) - 3


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def build_full_length_hdr(wt_reference: str, hdr_donor: str, min_anchor: int = 20) -> str:
    """
    Build a full-length HDR amplicon by inserting the edited region into WT.

    The HDR donor is a short sequence containing:
    - Left homology arm (matches WT)
    - Edited region (barcode + PAM mutation)
    - Right homology arm (matches WT)

    This function finds where the donor aligns in WT and builds a full-length
    sequence by replacing the corresponding WT region with the donor.

    Args:
        wt_reference: Full-length WT amplicon sequence
        hdr_donor: Short HDR donor template
        min_anchor: Minimum anchor length for alignment

    Returns:
        Full-length HDR amplicon sequence
    """
    wt = wt_reference.upper()
    donor = hdr_donor.upper()

    # Find left anchor (start of donor in WT)
    left_anchor = donor[:min_anchor]
    left_pos = wt.find(left_anchor)

    if left_pos == -1:
        # Try shorter anchor
        for anchor_len in range(min_anchor - 1, 9, -1):
            left_anchor = donor[:anchor_len]
            left_pos = wt.find(left_anchor)
            if left_pos != -1:
                break

    if left_pos == -1:
        logger.warning(f"Could not align HDR donor to WT reference (left anchor)")
        return donor  # Return donor as-is

    # Find right anchor (end of donor in WT)
    right_anchor = donor[-min_anchor:]
    right_pos = wt.find(right_anchor)

    if right_pos == -1:
        # Try shorter anchor
        for anchor_len in range(min_anchor - 1, 9, -1):
            right_anchor = donor[-anchor_len:]
            right_pos = wt.find(right_anchor)
            if right_pos != -1:
                break

    if right_pos == -1:
        logger.warning(f"Could not align HDR donor to WT reference (right anchor)")
        return donor  # Return donor as-is

    # Calculate the region to replace
    # The donor spans from left_pos to right_pos + len(right_anchor) in the WT
    # But the donor may have insertions/deletions relative to WT

    # Build full-length HDR:
    # WT[0:left_pos] + donor + WT[right_pos + len(right_anchor):]
    full_hdr = wt[:left_pos] + donor + wt[right_pos + len(right_anchor):]

    return full_hdr


def load_multi_reference_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Load multi-sequence FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dict mapping sequence_name -> sequence
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq).upper()
                # Start new sequence
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq).upper()

    return sequences


def build_emx1_multi_reference(
    reference: str,
    templates_fasta: Path,
    output_path: Path,
    guide_seq: str = "GAGTCCGAGCAGAAGAAGAA",
) -> Dict[str, int]:
    """
    Convenience function for EMX1 project.

    Builds multi-reference FASTA from existing templates FASTA.

    Args:
        reference: EMX1 WT reference sequence
        templates_fasta: Path to HDR templates FASTA (barcode_id -> template)
        output_path: Path to write multi-reference FASTA
        guide_seq: EMX1 guide sequence

    Returns:
        Dict mapping sequence_name -> cut_site_position
    """
    from .barcode_templates import load_templates_fasta

    # Load existing templates
    templates = load_templates_fasta(templates_fasta)

    # Build multi-reference
    return build_multi_reference_fasta(
        wt_reference=reference,
        hdr_templates=templates,
        output_path=output_path,
        guide_seq=guide_seq,
    )
