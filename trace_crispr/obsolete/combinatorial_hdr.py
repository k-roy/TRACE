"""
Combinatorial HDR detection for partial HDR analysis.

This module generates all possible partial HDR reference sequences
and classifies reads based on which combination of SNVs they contain.

This enables detection of:
- Complete HDR (all SNVs integrated)
- Partial HDR (subset of SNVs integrated)
- Gradient effects based on distance from cut site

Author: Kevin R. Roy
"""

import logging
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class SNVPosition:
    """A single SNV position in the HDR signature."""
    position: int  # 0-based position in reference
    wt_base: str
    hdr_base: str
    distance_to_cut: int = 0  # Signed distance (negative = 5' of cut)

    def __hash__(self):
        return hash((self.position, self.wt_base, self.hdr_base))

    def __eq__(self, other):
        return (self.position == other.position and
                self.wt_base == other.wt_base and
                self.hdr_base == other.hdr_base)


@dataclass
class HDRVariant:
    """An HDR variant with a specific subset of SNVs."""
    name: str
    sequence: str
    snvs_present: Tuple[SNVPosition, ...]  # Which SNVs are in this variant
    n_snvs: int

    @property
    def snv_positions(self) -> Tuple[int, ...]:
        return tuple(snv.position for snv in self.snvs_present)


@dataclass
class HDRClassification:
    """Classification result for a read with detailed SNV information."""
    outcome: str  # 'HDR_COMPLETE', 'HDR_PARTIAL', 'WT', 'NHEJ', etc.
    snvs_detected: Tuple[SNVPosition, ...]  # Which SNVs were found in read
    snvs_missing: Tuple[SNVPosition, ...]  # Which SNVs were NOT found
    n_snvs_detected: int
    n_snvs_total: int
    best_variant: Optional[str] = None  # Name of best-matching variant
    mismatches_to_best: int = 0
    distance_profile: Dict[int, bool] = field(default_factory=dict)  # position -> detected


def get_snv_positions(
    wt_seq: str,
    hdr_seq: str,
    cut_site: int,
    offset: int = 0
) -> List[SNVPosition]:
    """
    Extract SNV positions from WT and HDR sequences.

    Args:
        wt_seq: Wild-type reference sequence
        hdr_seq: HDR template sequence (must be same length as wt_seq)
        cut_site: Position of the cut site in the reference
        offset: Offset to add to positions (for aligned regions)

    Returns:
        List of SNVPosition objects
    """
    if len(wt_seq) != len(hdr_seq):
        raise ValueError(f"WT ({len(wt_seq)}) and HDR ({len(hdr_seq)}) sequences must be same length")

    snvs = []
    for i, (wt_base, hdr_base) in enumerate(zip(wt_seq.upper(), hdr_seq.upper())):
        if wt_base != hdr_base:
            pos = i + offset
            distance = pos - cut_site  # Negative = 5' of cut, positive = 3' of cut
            snvs.append(SNVPosition(
                position=pos,
                wt_base=wt_base,
                hdr_base=hdr_base,
                distance_to_cut=distance
            ))

    return snvs


def generate_hdr_variants(
    wt_seq: str,
    snv_positions: List[SNVPosition],
    max_combinations: int = 1000
) -> Dict[str, HDRVariant]:
    """
    Generate all possible HDR variant sequences.

    For N SNVs, generates 2^N - 1 variants (excluding WT).

    Args:
        wt_seq: Wild-type reference sequence
        snv_positions: List of SNV positions
        max_combinations: Maximum number of combinations to generate

    Returns:
        Dict mapping variant name to HDRVariant object
    """
    n_snvs = len(snv_positions)
    total_combinations = 2**n_snvs - 1

    if total_combinations > max_combinations:
        logger.warning(
            f"Too many SNV combinations ({total_combinations}). "
            f"Consider reducing SNVs or increasing max_combinations."
        )
        # Could implement sampling strategy here

    variants = {}

    # Add WT as reference
    variants['WT'] = HDRVariant(
        name='WT',
        sequence=wt_seq,
        snvs_present=tuple(),
        n_snvs=0
    )

    # Generate all combinations
    for n in range(1, n_snvs + 1):
        for combo in combinations(snv_positions, n):
            # Create sequence with this subset of SNVs
            seq = list(wt_seq.upper())
            for snv in combo:
                seq[snv.position] = snv.hdr_base

            # Create descriptive name
            positions_str = '_'.join(str(snv.position) for snv in combo)
            if n == n_snvs:
                name = 'HDR_COMPLETE'
            else:
                name = f'HDR_{n}of{n_snvs}_pos{positions_str}'

            variants[name] = HDRVariant(
                name=name,
                sequence=''.join(seq),
                snvs_present=combo,
                n_snvs=n
            )

    logger.info(f"Generated {len(variants)} HDR variants from {n_snvs} SNVs")
    return variants


def write_multi_reference_fasta(
    variants: Dict[str, HDRVariant],
    output_path: Path,
    include_metadata: bool = True
) -> None:
    """
    Write all HDR variants to a multi-reference FASTA file.

    Args:
        variants: Dict of variant name to HDRVariant
        output_path: Path to write FASTA file
        include_metadata: If True, include SNV positions in header
    """
    with open(output_path, 'w') as f:
        for name, variant in variants.items():
            if include_metadata and variant.snvs_present:
                positions = ','.join(str(snv.position) for snv in variant.snvs_present)
                header = f">{name} snvs={positions} n={variant.n_snvs}"
            else:
                header = f">{name}"

            f.write(f"{header}\n{variant.sequence}\n")

    logger.info(f"Wrote {len(variants)} reference sequences to {output_path}")


def classify_read_combinatorial(
    read_seq: str,
    variants: Dict[str, HDRVariant],
    snv_positions: List[SNVPosition],
    ref_start: int = 0
) -> HDRClassification:
    """
    Classify a read sequence against all HDR variants.

    Args:
        read_seq: The read sequence to classify
        variants: Dict of all HDR variants
        snv_positions: All SNV positions
        ref_start: Start position of read in reference coordinates

    Returns:
        HDRClassification with detailed SNV information
    """
    read_upper = read_seq.upper()

    # Find best matching variant
    best_variant = None
    best_mismatches = float('inf')

    for name, variant in variants.items():
        # Simple mismatch counting (could use alignment for gapped comparison)
        var_seq = variant.sequence

        # Handle different lengths by comparing overlapping region
        start = max(0, ref_start)
        end = min(len(var_seq), ref_start + len(read_upper))

        if end <= start:
            continue

        var_region = var_seq[start:end]
        read_region = read_upper[:end-start]

        mismatches = sum(1 for a, b in zip(var_region, read_region) if a != b)

        if mismatches < best_mismatches:
            best_mismatches = mismatches
            best_variant = name

    if best_variant is None:
        return HDRClassification(
            outcome='UNMAPPED',
            snvs_detected=tuple(),
            snvs_missing=tuple(snv_positions),
            n_snvs_detected=0,
            n_snvs_total=len(snv_positions)
        )

    # Determine which SNVs are present in the read
    snvs_detected = []
    snvs_missing = []
    distance_profile = {}

    for snv in snv_positions:
        read_pos = snv.position - ref_start

        if 0 <= read_pos < len(read_upper):
            read_base = read_upper[read_pos]
            if read_base == snv.hdr_base:
                snvs_detected.append(snv)
                distance_profile[snv.distance_to_cut] = True
            else:
                snvs_missing.append(snv)
                distance_profile[snv.distance_to_cut] = False
        else:
            # SNV position not covered by read
            snvs_missing.append(snv)

    # Determine outcome
    n_detected = len(snvs_detected)
    n_total = len(snv_positions)

    if n_detected == 0:
        outcome = 'WT'
    elif n_detected == n_total:
        outcome = 'HDR_COMPLETE'
    else:
        outcome = 'HDR_PARTIAL'

    return HDRClassification(
        outcome=outcome,
        snvs_detected=tuple(snvs_detected),
        snvs_missing=tuple(snvs_missing),
        n_snvs_detected=n_detected,
        n_snvs_total=n_total,
        best_variant=best_variant,
        mismatches_to_best=best_mismatches,
        distance_profile=distance_profile
    )


def analyze_conversion_tract(
    classifications: List[HDRClassification],
    snv_positions: List[SNVPosition]
) -> Dict:
    """
    Analyze conversion tract patterns across all reads.

    Returns statistics on:
    - Per-SNV integration frequency
    - Gradient from cut site
    - Correlation between adjacent SNVs

    Args:
        classifications: List of HDRClassification results
        snv_positions: All SNV positions

    Returns:
        Dict with analysis results
    """
    # Count detections per SNV
    snv_counts = {snv.position: 0 for snv in snv_positions}
    total_hdr = 0

    for clf in classifications:
        if clf.outcome in ('HDR_COMPLETE', 'HDR_PARTIAL'):
            total_hdr += 1
            for snv in clf.snvs_detected:
                snv_counts[snv.position] += 1

    # Calculate frequencies
    snv_frequencies = {}
    for snv in snv_positions:
        freq = snv_counts[snv.position] / total_hdr if total_hdr > 0 else 0
        snv_frequencies[snv.position] = {
            'frequency': freq,
            'count': snv_counts[snv.position],
            'distance_to_cut': snv.distance_to_cut,
            'wt_base': snv.wt_base,
            'hdr_base': snv.hdr_base
        }

    # Analyze co-occurrence patterns
    co_occurrence = {}
    for i, snv1 in enumerate(snv_positions):
        for snv2 in snv_positions[i+1:]:
            key = (snv1.position, snv2.position)
            both_present = 0
            either_present = 0

            for clf in classifications:
                pos1_present = snv1 in clf.snvs_detected
                pos2_present = snv2 in clf.snvs_detected

                if pos1_present and pos2_present:
                    both_present += 1
                if pos1_present or pos2_present:
                    either_present += 1

            co_occurrence[key] = {
                'both': both_present,
                'either': either_present,
                'jaccard': both_present / either_present if either_present > 0 else 0
            }

    # Partial HDR breakdown
    partial_breakdown = {}
    for clf in classifications:
        if clf.outcome == 'HDR_PARTIAL':
            n = clf.n_snvs_detected
            partial_breakdown[n] = partial_breakdown.get(n, 0) + 1

    return {
        'total_hdr_reads': total_hdr,
        'snv_frequencies': snv_frequencies,
        'co_occurrence': co_occurrence,
        'partial_breakdown': partial_breakdown,
        'complete_hdr_count': sum(1 for c in classifications if c.outcome == 'HDR_COMPLETE'),
        'partial_hdr_count': sum(1 for c in classifications if c.outcome == 'HDR_PARTIAL'),
    }
