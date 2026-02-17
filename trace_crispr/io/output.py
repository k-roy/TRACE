"""
Output generation for TRACE results.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import pandas as pd
import logging

logger = logging.getLogger(__name__)


@dataclass
class SampleResult:
    """Results for a single sample."""
    sample_id: str
    total_reads: int = 0
    aligned_reads: int = 0
    classifiable_reads: int = 0
    duplicate_rate: float = 0.0

    # Deduplicated counts
    dedup_wt_count: int = 0
    dedup_hdr_count: int = 0
    dedup_nhej_count: int = 0
    dedup_large_del_count: int = 0

    # Non-deduplicated counts (for comparison)
    nondedup_hdr_count: int = 0
    nondedup_nhej_count: int = 0

    # K-mer method results
    kmer_wt_count: int = 0
    kmer_hdr_count: int = 0
    kmer_contamination_count: int = 0

    # CRISPResso results
    crispresso_hdr_rate: Optional[float] = None
    crispresso_nhej_rate: Optional[float] = None

    # Metadata
    metadata: Dict = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

    @property
    def dedup_total(self) -> int:
        return self.dedup_wt_count + self.dedup_hdr_count + self.dedup_nhej_count + self.dedup_large_del_count

    @property
    def dedup_wt_pct(self) -> float:
        total = self.dedup_total
        return (self.dedup_wt_count / total * 100) if total > 0 else 0

    @property
    def dedup_hdr_pct(self) -> float:
        total = self.dedup_total
        return (self.dedup_hdr_count / total * 100) if total > 0 else 0

    @property
    def dedup_nhej_pct(self) -> float:
        total = self.dedup_total
        return (self.dedup_nhej_count / total * 100) if total > 0 else 0

    @property
    def dedup_large_del_pct(self) -> float:
        total = self.dedup_total
        return (self.dedup_large_del_count / total * 100) if total > 0 else 0

    @property
    def nondedup_total(self) -> int:
        return self.aligned_reads

    @property
    def nondedup_hdr_pct(self) -> float:
        total = self.nondedup_total
        return (self.nondedup_hdr_count / total * 100) if total > 0 else 0

    @property
    def nondedup_nhej_pct(self) -> float:
        total = self.nondedup_total
        return (self.nondedup_nhej_count / total * 100) if total > 0 else 0

    @property
    def kmer_hdr_rate(self) -> float:
        total = self.kmer_wt_count + self.kmer_hdr_count
        return (self.kmer_hdr_count / total) if total > 0 else 0


@dataclass
class GranularResult:
    """
    Per-template result for a single sample.

    Used for multi-template analysis where we track outcomes for each
    possible barcode/template separately.
    """
    sample_id: str
    template_id: str  # Template/barcode ID, or 'NA' for non-template outcomes
    outcome: str  # 'HDR_PERFECT', 'HDR_IMPERFECT', 'NHEJ', 'LARGE_DELETION', 'WT'
    count: int
    rate: float  # Percentage of total classifiable reads
    is_expected: bool  # True if this template is the expected one for this sample
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for output."""
        return {
            'sample_id': self.sample_id,
            'template_id': self.template_id,
            'outcome': self.outcome,
            'count': self.count,
            'rate': f"{self.rate:.4f}",
            'is_expected': str(self.is_expected).upper() if self.template_id != 'NA' else 'NA',
            **self.metadata,
        }


@dataclass
class MultiTemplateSampleResult:
    """
    Complete results for a sample with multiple possible templates.

    Tracks HDR outcomes per-template, allowing detection of expected
    barcodes and contaminating barcodes.
    """
    sample_id: str
    total_reads: int
    aligned_reads: int

    # Per-template HDR counts
    hdr_counts_by_template: Dict[str, int] = field(default_factory=dict)

    # HDR subtypes (perfect vs imperfect)
    hdr_perfect_by_template: Dict[str, int] = field(default_factory=dict)
    hdr_imperfect_by_template: Dict[str, int] = field(default_factory=dict)

    # Non-template outcomes
    wt_count: int = 0
    nhej_count: int = 0
    large_del_count: int = 0
    ambiguous_count: int = 0

    # Expected template (if specified in sample key)
    expected_template: Optional[str] = None

    # Metadata from sample key
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def classifiable_reads(self) -> int:
        """Total classifiable reads (WT + all HDR + NHEJ + LgDel)."""
        return (
            self.wt_count +
            sum(self.hdr_counts_by_template.values()) +
            self.nhej_count +
            self.large_del_count
        )

    @property
    def total_hdr_count(self) -> int:
        """Sum of HDR counts across all templates."""
        return sum(self.hdr_counts_by_template.values())

    @property
    def total_hdr_rate(self) -> float:
        """HDR rate as percentage of classifiable reads."""
        total = self.classifiable_reads
        return (self.total_hdr_count / total * 100) if total > 0 else 0

    @property
    def wt_rate(self) -> float:
        """WT rate as percentage of classifiable reads."""
        total = self.classifiable_reads
        return (self.wt_count / total * 100) if total > 0 else 0

    @property
    def nhej_rate(self) -> float:
        """NHEJ rate as percentage of classifiable reads."""
        total = self.classifiable_reads
        return (self.nhej_count / total * 100) if total > 0 else 0

    @property
    def large_del_rate(self) -> float:
        """Large deletion rate as percentage of classifiable reads."""
        total = self.classifiable_reads
        return (self.large_del_count / total * 100) if total > 0 else 0

    def hdr_rate_for_template(self, template_id: str) -> float:
        """HDR rate for a specific template as percentage."""
        total = self.classifiable_reads
        count = self.hdr_counts_by_template.get(template_id, 0)
        return (count / total * 100) if total > 0 else 0

    @property
    def expected_template_rate(self) -> Optional[float]:
        """HDR rate for the expected template (if specified)."""
        if self.expected_template is None:
            return None
        return self.hdr_rate_for_template(self.expected_template)

    def get_granular_results(self) -> List[GranularResult]:
        """
        Convert to list of GranularResult objects.

        Returns one GranularResult per (template, outcome) combination.
        """
        results = []
        total = self.classifiable_reads

        # Add WT result
        if self.wt_count > 0:
            results.append(GranularResult(
                sample_id=self.sample_id,
                template_id='NA',
                outcome='WT',
                count=self.wt_count,
                rate=(self.wt_count / total * 100) if total > 0 else 0,
                is_expected=False,
                metadata=self.metadata,
            ))

        # Add per-template HDR results
        for template_id in sorted(self.hdr_counts_by_template.keys()):
            perfect_count = self.hdr_perfect_by_template.get(template_id, 0)
            imperfect_count = self.hdr_imperfect_by_template.get(template_id, 0)
            total_hdr = self.hdr_counts_by_template.get(template_id, 0)
            is_expected = (template_id == self.expected_template)

            # If we have perfect/imperfect breakdown, output separately
            if perfect_count > 0 or imperfect_count > 0:
                if perfect_count > 0:
                    results.append(GranularResult(
                        sample_id=self.sample_id,
                        template_id=template_id,
                        outcome='HDR_PERFECT',
                        count=perfect_count,
                        rate=(perfect_count / total * 100) if total > 0 else 0,
                        is_expected=is_expected,
                        metadata=self.metadata,
                    ))
                if imperfect_count > 0:
                    results.append(GranularResult(
                        sample_id=self.sample_id,
                        template_id=template_id,
                        outcome='HDR_IMPERFECT',
                        count=imperfect_count,
                        rate=(imperfect_count / total * 100) if total > 0 else 0,
                        is_expected=is_expected,
                        metadata=self.metadata,
                    ))
            elif total_hdr > 0:
                # No perfect/imperfect breakdown, just output total HDR
                results.append(GranularResult(
                    sample_id=self.sample_id,
                    template_id=template_id,
                    outcome='HDR',
                    count=total_hdr,
                    rate=(total_hdr / total * 100) if total > 0 else 0,
                    is_expected=is_expected,
                    metadata=self.metadata,
                ))

        # Add NHEJ result
        if self.nhej_count > 0:
            results.append(GranularResult(
                sample_id=self.sample_id,
                template_id='NA',
                outcome='NHEJ',
                count=self.nhej_count,
                rate=(self.nhej_count / total * 100) if total > 0 else 0,
                is_expected=False,
                metadata=self.metadata,
            ))

        # Add large deletion result
        if self.large_del_count > 0:
            results.append(GranularResult(
                sample_id=self.sample_id,
                template_id='NA',
                outcome='LARGE_DELETION',
                count=self.large_del_count,
                rate=(self.large_del_count / total * 100) if total > 0 else 0,
                is_expected=False,
                metadata=self.metadata,
            ))

        return results


def write_results_tsv(
    results: List[SampleResult],
    output_path: Path,
    include_metadata: bool = True,
) -> Path:
    """
    Write results to TSV file.

    Args:
        results: List of SampleResult objects
        output_path: Path for output TSV
        include_metadata: Include metadata columns

    Returns:
        Path to written file
    """
    rows = []

    for r in results:
        row = {
            'sample': r.sample_id,
            'total_reads': r.total_reads,
            'aligned_reads': r.aligned_reads,
            'classifiable_reads': r.classifiable_reads,
            'duplicate_rate': f"{r.duplicate_rate:.4f}",
            'Dedup_WT_%': f"{r.dedup_wt_pct:.2f}",
            'Dedup_HDR_%': f"{r.dedup_hdr_pct:.2f}",
            'Dedup_NHEJ_%': f"{r.dedup_nhej_pct:.2f}",
            'Dedup_LgDel_%': f"{r.dedup_large_del_pct:.3f}",
            'NonDedup_HDR_%': f"{r.nondedup_hdr_pct:.2f}",
            'NonDedup_NHEJ_%': f"{r.nondedup_nhej_pct:.2f}",
            'kmer_hdr_rate': f"{r.kmer_hdr_rate:.4f}",
        }

        # Add CRISPResso results if available
        if r.crispresso_hdr_rate is not None:
            row['crispresso_hdr_rate'] = f"{r.crispresso_hdr_rate:.4f}"
        if r.crispresso_nhej_rate is not None:
            row['crispresso_nhej_rate'] = f"{r.crispresso_nhej_rate:.4f}"

        # Add metadata
        if include_metadata and r.metadata:
            for k, v in r.metadata.items():
                row[k] = v

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)

    logger.info(f"Wrote results to {output_path}")

    return output_path


def write_per_read_classifications(
    classifications: List[Dict],
    output_path: Path,
) -> Path:
    """
    Write per-read classification details to TSV.

    Args:
        classifications: List of classification dicts
        output_path: Path for output TSV

    Returns:
        Path to written file
    """
    df = pd.DataFrame(classifications)
    df.to_csv(output_path, sep='\t', index=False)

    logger.info(f"Wrote {len(classifications)} read classifications to {output_path}")

    return output_path


def generate_summary_report(
    results: List[SampleResult],
    output_path: Path,
) -> Path:
    """
    Generate a summary report in markdown format.

    Args:
        results: List of SampleResult objects
        output_path: Path for output markdown file

    Returns:
        Path to written file
    """
    with open(output_path, 'w') as f:
        f.write("# CRISPRo Analysis Summary\n\n")

        # Overall statistics
        total_samples = len(results)
        total_reads = sum(r.total_reads for r in results)
        total_aligned = sum(r.aligned_reads for r in results)

        f.write("## Overview\n\n")
        f.write(f"- **Samples analyzed:** {total_samples}\n")
        f.write(f"- **Total reads:** {total_reads:,}\n")
        f.write(f"- **Aligned reads:** {total_aligned:,} ({total_aligned/total_reads*100:.1f}%)\n\n")

        # HDR summary
        hdr_rates = [r.dedup_hdr_pct for r in results]
        if hdr_rates:
            f.write("## HDR Rates (Deduplicated)\n\n")
            f.write(f"- **Mean:** {sum(hdr_rates)/len(hdr_rates):.2f}%\n")
            f.write(f"- **Max:** {max(hdr_rates):.2f}%\n")
            f.write(f"- **Min:** {min(hdr_rates):.2f}%\n\n")

        # Top samples by HDR
        f.write("## Top 10 Samples by HDR Rate\n\n")
        f.write("| Sample | HDR % | NHEJ % | WT % |\n")
        f.write("|--------|-------|--------|------|\n")

        sorted_results = sorted(results, key=lambda r: r.dedup_hdr_pct, reverse=True)
        for r in sorted_results[:10]:
            f.write(f"| {r.sample_id} | {r.dedup_hdr_pct:.2f} | {r.dedup_nhej_pct:.2f} | {r.dedup_wt_pct:.2f} |\n")

        f.write("\n")

    logger.info(f"Wrote summary report to {output_path}")

    return output_path


def write_multi_template_summary_tsv(
    results: List[MultiTemplateSampleResult],
    output_path: Path,
) -> Path:
    """
    Write multi-template summary results to TSV file.

    One row per sample with overall editing rates and expected template rate.

    Args:
        results: List of MultiTemplateSampleResult objects
        output_path: Path for output TSV

    Returns:
        Path to written file
    """
    rows = []

    for r in results:
        row = {
            'sample_id': r.sample_id,
            'total_reads': r.total_reads,
            'aligned_reads': r.aligned_reads,
            'classifiable_reads': r.classifiable_reads,
            'WT_%': f"{r.wt_rate:.2f}",
            'HDR_%': f"{r.total_hdr_rate:.2f}",
            'NHEJ_%': f"{r.nhej_rate:.2f}",
            'LgDel_%': f"{r.large_del_rate:.3f}",
            'expected_template': r.expected_template or 'NA',
            'expected_template_%': f"{r.expected_template_rate:.2f}" if r.expected_template_rate is not None else 'NA',
            'ambiguous_count': r.ambiguous_count,
        }

        # Add metadata
        for k, v in r.metadata.items():
            row[k] = v

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)

    logger.info(f"Wrote multi-template summary to {output_path}")

    return output_path


def write_granular_results_tsv(
    results: List[MultiTemplateSampleResult],
    output_path: Path,
) -> Path:
    """
    Write granular per-template results to TSV file.

    One row per (sample, template, outcome) combination.

    Args:
        results: List of MultiTemplateSampleResult objects
        output_path: Path for output TSV

    Returns:
        Path to written file
    """
    rows = []

    for r in results:
        granular = r.get_granular_results()
        for g in granular:
            rows.append(g.to_dict())

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)

    logger.info(f"Wrote {len(rows)} granular results to {output_path}")

    return output_path


def write_multi_template_results(
    results: List[MultiTemplateSampleResult],
    output_dir: Path,
    prefix: str = "",
) -> Dict[str, Path]:
    """
    Write all multi-template result files.

    Args:
        results: List of MultiTemplateSampleResult objects
        output_dir: Directory for output files
        prefix: Optional prefix for output filenames

    Returns:
        Dict mapping file type to path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}

    # Summary table
    summary_path = output_dir / f"{prefix}per_sample_editing_outcomes_all_methods.tsv"
    paths['summary'] = write_multi_template_summary_tsv(results, summary_path)

    # Granular table
    granular_path = output_dir / f"{prefix}per_sample_per_template_outcomes.tsv"
    paths['granular'] = write_granular_results_tsv(results, granular_path)

    return paths
