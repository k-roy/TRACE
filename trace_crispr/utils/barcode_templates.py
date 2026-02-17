"""
Barcode template generation utilities.

Generate HDR templates for barcoded CRISPR editing experiments.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, Optional, List
import pandas as pd
import logging
import re

logger = logging.getLogger(__name__)


def is_valid_barcode(barcode) -> bool:
    """
    Check if barcode is a valid DNA sequence.

    Args:
        barcode: Value to check (can be any type)

    Returns:
        True if barcode is a valid DNA sequence (only ACGT), False otherwise
    """
    if barcode is None or pd.isna(barcode):
        return False
    barcode_str = str(barcode).strip().upper()
    # Must be non-empty and not a placeholder value
    if not barcode_str or barcode_str in ('NA', 'EMPTY', 'NONE', 'NULL', 'CONTROL', 'N/A'):
        return False
    # Check that it's a valid DNA sequence (only ACGT)
    return bool(re.match(r'^[ACGT]+$', barcode_str))


def generate_barcoded_hdr_templates(
    left_homology_arm: str,
    right_homology_arm: str,
    barcodes: Dict[str, str],
    output_fasta: Optional[Path] = None,
) -> Dict[str, str]:
    """
    Generate HDR templates for all barcodes.

    Each template = left_homology_arm + barcode + right_homology_arm

    Args:
        left_homology_arm: Sequence upstream of barcode insertion
        right_homology_arm: Sequence downstream of barcode insertion
        barcodes: Dict mapping barcode_id → barcode sequence
        output_fasta: Optional path to write FASTA file

    Returns:
        Dict mapping barcode_id → full HDR template sequence

    Example:
        >>> templates = generate_barcoded_hdr_templates(
        ...     left_homology_arm="ATCGATCG",
        ...     right_homology_arm="GCTAGCTA",
        ...     barcodes={"bc1": "AAAA", "bc2": "TTTT"}
        ... )
        >>> templates["bc1"]
        'ATCGATCGAAAAGCTAGCTA'
    """
    templates = {}

    for barcode_id, barcode_seq in barcodes.items():
        template = left_homology_arm.upper() + barcode_seq.upper() + right_homology_arm.upper()
        templates[barcode_id] = template

    logger.info(f"Generated {len(templates)} HDR templates")
    logger.info(f"  Template length: {len(left_homology_arm)} + {len(next(iter(barcodes.values())))} + {len(right_homology_arm)} bp")

    if output_fasta:
        write_templates_fasta(templates, output_fasta)

    return templates


def write_templates_fasta(
    templates: Dict[str, str],
    output_path: Path,
) -> None:
    """
    Write HDR templates to FASTA file.

    Args:
        templates: Dict mapping template_id → sequence
        output_path: Path to write FASTA file
    """
    with open(output_path, 'w') as f:
        for template_id, sequence in sorted(templates.items()):
            f.write(f">{template_id}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(sequence), 80):
                f.write(f"{sequence[i:i+80]}\n")

    logger.info(f"Wrote {len(templates)} templates to {output_path}")


def load_templates_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Load HDR templates from FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dict mapping template_id → sequence
    """
    templates = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    templates[current_id] = ''.join(current_seq).upper()

                # Start new sequence
                current_id = line[1:].split()[0]  # Take first word as ID
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None:
            templates[current_id] = ''.join(current_seq).upper()

    logger.info(f"Loaded {len(templates)} templates from {fasta_path}")

    return templates


def load_barcodes_from_tsv(
    tsv_path: Path,
    barcode_column: str = 'barcode',
    id_column: Optional[str] = None,
) -> Dict[str, str]:
    """
    Load barcodes from TSV file.

    Args:
        tsv_path: Path to TSV file
        barcode_column: Column containing barcode sequences
        id_column: Optional column for barcode IDs (defaults to barcode seq)

    Returns:
        Dict mapping barcode_id → barcode sequence
    """
    df = pd.read_csv(tsv_path, sep='\t')

    if barcode_column not in df.columns:
        raise ValueError(f"Column '{barcode_column}' not found in {tsv_path}")

    # Filter valid barcodes (must be valid DNA sequences)
    valid = df[df[barcode_column].apply(is_valid_barcode)].copy()

    if id_column and id_column in df.columns:
        # Use specified ID column
        barcodes = dict(zip(valid[id_column].astype(str), valid[barcode_column].astype(str)))
    else:
        # Use barcode sequence as ID
        unique_barcodes = valid[barcode_column].unique()
        barcodes = {str(bc): str(bc) for bc in unique_barcodes}

    logger.info(f"Loaded {len(barcodes)} barcodes from {tsv_path}")

    return barcodes


def generate_templates_from_keyfiles(
    sample_key_path: Path,
    guide_donor_info_path: Path,
    output_fasta: Path,
    barcode_column: str = 'barcode',
) -> Dict[str, str]:
    """
    Generate HDR templates from project keyfiles.

    Combines sample_key (for barcodes) with guide_donor_info (for homology arms).

    Args:
        sample_key_path: Path to sample_key.tsv
        guide_donor_info_path: Path to guide_donor_and_reference_info.tsv
        output_fasta: Path to write output FASTA
        barcode_column: Column in sample_key containing barcodes

    Returns:
        Dict mapping barcode_id → full HDR template sequence
    """
    # Load homology arms from guide_donor_info
    guide_donor = pd.read_csv(guide_donor_info_path, sep='\t')

    if len(guide_donor) == 0:
        raise ValueError(f"No data in {guide_donor_info_path}")

    # Take first row (assuming single locus)
    row = guide_donor.iloc[0]
    left_arm = str(row['left_homology_arm']).upper()
    right_arm = str(row['right_homology_arm']).upper()

    logger.info(f"Loaded homology arms from {guide_donor_info_path}")
    logger.info(f"  Left arm: {len(left_arm)} bp")
    logger.info(f"  Right arm: {len(right_arm)} bp")

    # Load barcodes from sample_key
    barcodes = load_barcodes_from_tsv(sample_key_path, barcode_column)

    # Generate templates
    templates = generate_barcoded_hdr_templates(
        left_homology_arm=left_arm,
        right_homology_arm=right_arm,
        barcodes=barcodes,
        output_fasta=output_fasta,
    )

    return templates
