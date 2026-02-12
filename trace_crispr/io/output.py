"""
Output generation for CRISPRo results.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
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
