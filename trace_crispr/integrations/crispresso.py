"""
CRISPResso2 integration wrapper.

Author: Kevin R. Roy
"""

import json
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


@dataclass
class CRISPRessoResult:
    """Result from CRISPResso2 analysis."""
    success: bool
    output_dir: Path
    hdr_rate: Optional[float] = None
    nhej_rate: Optional[float] = None
    wt_rate: Optional[float] = None
    total_reads: int = 0
    aligned_reads: int = 0
    error_message: Optional[str] = None
    raw_results: Optional[Dict] = None


class CRISPRessoRunner:
    """Handle CRISPResso2 execution via conda, module, or container."""

    def __init__(self, container_path: Optional[Path] = None):
        self.container_path = container_path
        self.execution_mode = self._detect_mode()

    def _detect_mode(self) -> str:
        """Detect how to run CRISPResso2."""
        # 1. Check for container
        if self.container_path and self.container_path.exists():
            return 'singularity'

        # 2. Check for direct installation
        try:
            subprocess.run(
                ['CRISPResso', '--version'],
                capture_output=True, timeout=10
            )
            return 'direct'
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        return 'unavailable'

    def is_available(self) -> bool:
        """Check if CRISPResso2 is available."""
        return self.execution_mode != 'unavailable'

    def run(
        self,
        r1_fastq: Path,
        r2_fastq: Optional[Path],
        amplicon_seq: str,
        guide_seq: str,
        output_dir: Path,
        hdr_seq: Optional[str] = None,
        name: str = "CRISPResso",
        quantification_window: int = 10,
        min_identity: int = 70,
        **kwargs
    ) -> CRISPRessoResult:
        """
        Run CRISPResso2 analysis.

        Args:
            r1_fastq: Path to R1 FASTQ
            r2_fastq: Optional path to R2 FASTQ
            amplicon_seq: Amplicon sequence
            guide_seq: Guide sequence (without PAM)
            output_dir: Output directory
            hdr_seq: Optional HDR expected sequence
            name: Analysis name
            quantification_window: Window size for quantification
            min_identity: Minimum alignment identity score
            **kwargs: Additional CRISPResso parameters

        Returns:
            CRISPRessoResult with analysis results
        """
        if not self.is_available():
            return CRISPRessoResult(
                success=False,
                output_dir=output_dir,
                error_message="CRISPResso2 is not available"
            )

        # Build command
        cmd = self._build_command(
            r1_fastq=r1_fastq,
            r2_fastq=r2_fastq,
            amplicon_seq=amplicon_seq,
            guide_seq=guide_seq,
            output_dir=output_dir,
            hdr_seq=hdr_seq,
            name=name,
            quantification_window=quantification_window,
            min_identity=min_identity,
            **kwargs
        )

        try:
            logger.info(f"Running CRISPResso2: {' '.join(cmd[:5])}...")
            subprocess.run(
                cmd,
                capture_output=True,
                check=True,
                cwd=str(output_dir.parent)
            )

            # Parse results
            return self._parse_results(output_dir, name)

        except subprocess.CalledProcessError as e:
            return CRISPRessoResult(
                success=False,
                output_dir=output_dir,
                error_message=e.stderr.decode() if e.stderr else str(e)
            )
        except Exception as e:
            return CRISPRessoResult(
                success=False,
                output_dir=output_dir,
                error_message=str(e)
            )

    def _build_command(
        self,
        r1_fastq: Path,
        r2_fastq: Optional[Path],
        amplicon_seq: str,
        guide_seq: str,
        output_dir: Path,
        hdr_seq: Optional[str],
        name: str,
        quantification_window: int,
        min_identity: int,
        **kwargs
    ) -> list:
        """Build CRISPResso command."""
        if self.execution_mode == 'singularity':
            cmd = ['singularity', 'exec', str(self.container_path), 'CRISPResso']
        else:
            cmd = ['CRISPResso']

        # Add required arguments
        cmd.extend([
            '--fastq_r1', str(r1_fastq),
            '--amplicon_seq', amplicon_seq,
            '--guide_seq', guide_seq,
            '--output_folder', str(output_dir),
            '--name', name,
            '--quantification_window_size', str(quantification_window),
            '--min_identity_score', str(min_identity),
        ])

        # Add R2 if provided
        if r2_fastq:
            cmd.extend(['--fastq_r2', str(r2_fastq)])

        # Add HDR sequence if provided
        if hdr_seq:
            cmd.extend(['--expected_hdr_amplicon_seq', hdr_seq])

        # Add any additional kwargs
        for key, value in kwargs.items():
            if value is not None:
                cmd.extend([f'--{key}', str(value)])

        return cmd

    def _parse_results(self, output_dir: Path, name: str) -> CRISPRessoResult:
        """Parse CRISPResso2 output files."""
        # Find the output folder
        crispresso_dir = output_dir / f"CRISPResso_on_{name}"

        if not crispresso_dir.exists():
            # Try alternate naming
            possible_dirs = list(output_dir.glob("CRISPResso_on_*"))
            if possible_dirs:
                crispresso_dir = possible_dirs[0]
            else:
                return CRISPRessoResult(
                    success=False,
                    output_dir=output_dir,
                    error_message="CRISPResso output directory not found"
                )

        # Try to read mapping statistics
        mapping_stats = crispresso_dir / "Mapping_statistics.txt"

        raw_results = {}
        total_reads = 0
        aligned_reads = 0
        hdr_rate = None
        nhej_rate = None
        wt_rate = None

        # Parse mapping statistics
        if mapping_stats.exists():
            with open(mapping_stats) as f:
                for line in f:
                    if ':' in line:
                        key, value = line.strip().split(':', 1)
                        raw_results[key.strip()] = value.strip()
                        if 'READS IN INPUTS' in key.upper():
                            try:
                                total_reads = int(value.strip())
                            except ValueError:
                                pass
                        if 'READS ALIGNED' in key.upper():
                            try:
                                aligned_reads = int(value.strip())
                            except ValueError:
                                pass

        # Try to parse JSON summary if available
        json_file = crispresso_dir / "CRISPResso2_info.json"
        if json_file.exists():
            try:
                with open(json_file) as f:
                    json_data = json.load(f)
                    raw_results['json'] = json_data

                    # Extract rates from JSON
                    results = json_data.get('results', {}).get('alignment_stats', {})
                    if 'hdr' in results:
                        hdr_rate = results['hdr'].get('pct', 0) / 100
                    if 'nhej' in results:
                        nhej_rate = results['nhej'].get('pct', 0) / 100
                    if 'unmodified' in results:
                        wt_rate = results['unmodified'].get('pct', 0) / 100

            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Could not parse CRISPResso JSON: {e}")

        # Parse alleles frequency table as fallback
        alleles_file = crispresso_dir / "Alleles_frequency_table.zip"
        if alleles_file.exists() and hdr_rate is None:
            # Would need to unzip and parse - skip for now
            pass

        return CRISPRessoResult(
            success=True,
            output_dir=crispresso_dir,
            hdr_rate=hdr_rate,
            nhej_rate=nhej_rate,
            wt_rate=wt_rate,
            total_reads=total_reads,
            aligned_reads=aligned_reads,
            raw_results=raw_results
        )


def run_crispresso_batch(
    samples: list,
    amplicon_seq: str,
    guide_seq: str,
    output_dir: Path,
    hdr_seq: Optional[str] = None,
    threads: int = 4,
    container_path: Optional[Path] = None,
) -> Dict[str, CRISPRessoResult]:
    """
    Run CRISPResso2 on multiple samples.

    Args:
        samples: List of dicts with 'sample_id', 'r1_path', 'r2_path' keys
        amplicon_seq: Amplicon sequence
        guide_seq: Guide sequence
        output_dir: Base output directory
        hdr_seq: Optional HDR sequence
        threads: Number of threads
        container_path: Optional path to Singularity container

    Returns:
        Dict mapping sample_id to CRISPRessoResult
    """
    runner = CRISPRessoRunner(container_path)

    if not runner.is_available():
        logger.warning("CRISPResso2 is not available, skipping")
        return {}

    results = {}

    for sample in samples:
        sample_id = sample['sample_id']
        sample_output = output_dir / sample_id

        logger.info(f"Running CRISPResso2 for {sample_id}")

        result = runner.run(
            r1_fastq=Path(sample['r1_path']),
            r2_fastq=Path(sample['r2_path']) if sample.get('r2_path') else None,
            amplicon_seq=amplicon_seq,
            guide_seq=guide_seq,
            output_dir=sample_output,
            hdr_seq=hdr_seq,
            name=sample_id,
        )

        results[sample_id] = result

    return results
