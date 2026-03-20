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
    output_dir: Path,
    default_amplicon_seq: Optional[str] = None,
    default_guide_seq: Optional[str] = None,
    default_hdr_seq: Optional[str] = None,
    threads: int = 4,
    skip_finished: bool = True,
    container_path: Optional[Path] = None,
    quantification_window: int = 10,
    min_identity: int = 70,
) -> Dict[str, CRISPRessoResult]:
    """
    Run CRISPRessoBatch on multiple samples efficiently.

    Uses CRISPRessoBatch for:
    - Sequence caching (identical reads aligned once across all samples)
    - Parallel processing (-p flag)
    - Skip already-completed samples (--skip_finished)

    Args:
        samples: List of dicts with keys:
            - 'sample_id': Sample name (required)
            - 'r1_path': Path to R1 FASTQ (required)
            - 'r2_path': Path to R2 FASTQ (optional)
            - 'amplicon_seq': Per-sample amplicon (optional, uses default)
            - 'guide_seq': Per-sample guide (optional, uses default)
            - 'hdr_seq': Per-sample HDR template (optional, uses default)
        output_dir: Base output directory
        default_amplicon_seq: Default amplicon sequence
        default_guide_seq: Default guide sequence
        default_hdr_seq: Default HDR sequence
        threads: Number of parallel processes
        skip_finished: Skip already-completed samples
        container_path: Optional path to Singularity container

    Returns:
        Dict mapping sample_id to CRISPRessoResult
    """
    runner = CRISPRessoRunner(container_path)

    if not runner.is_available():
        logger.warning("CRISPResso2 is not available, skipping")
        return {}

    if not samples:
        logger.warning("No samples provided for CRISPRessoBatch")
        return {}

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate batch settings file
    batch_file = output_dir / "crispresso_batch_settings.tsv"
    _write_batch_settings(
        batch_file=batch_file,
        samples=samples,
        default_amplicon_seq=default_amplicon_seq,
        default_guide_seq=default_guide_seq,
        default_hdr_seq=default_hdr_seq,
        quantification_window=quantification_window,
        min_identity=min_identity,
    )

    # Build CRISPRessoBatch command
    cmd = _build_batch_command(
        runner=runner,
        batch_file=batch_file,
        output_dir=output_dir,
        threads=threads,
        skip_finished=skip_finished,
    )

    logger.info(f"Running CRISPRessoBatch on {len(samples)} samples with {threads} threads...")

    try:
        subprocess.run(
            cmd,
            capture_output=True,
            check=True,
            cwd=str(output_dir),
        )
        logger.info("CRISPRessoBatch completed successfully")

    except subprocess.CalledProcessError as e:
        logger.error(f"CRISPRessoBatch failed: {e.stderr.decode() if e.stderr else str(e)}")
        # Continue to parse whatever results exist

    # Parse results for each sample
    results = {}
    for sample in samples:
        sample_id = sample['sample_id']
        sample_result = runner._parse_results(output_dir, sample_id)
        results[sample_id] = sample_result

    return results


def _write_batch_settings(
    batch_file: Path,
    samples: list,
    default_amplicon_seq: Optional[str],
    default_guide_seq: Optional[str],
    default_hdr_seq: Optional[str],
    quantification_window: int,
    min_identity: int,
) -> None:
    """
    Write CRISPRessoBatch settings TSV file.

    Each row can have per-sample parameters that override defaults.
    """
    # Determine which columns we need
    has_r2 = any(sample.get('r2_path') for sample in samples)
    _has_per_sample_amplicon = any(sample.get('amplicon_seq') for sample in samples)  # noqa: F841
    _has_per_sample_guide = any(sample.get('guide_seq') for sample in samples)  # noqa: F841
    has_per_sample_hdr = any(sample.get('hdr_seq') for sample in samples)

    # Build header
    columns = ['name', 'fastq_r1']
    if has_r2:
        columns.append('fastq_r2')
    columns.append('amplicon_seq')
    columns.append('guide_seq')
    if default_hdr_seq or has_per_sample_hdr:
        columns.append('expected_hdr_amplicon_seq')
    columns.extend(['quantification_window_size', 'min_identity_score'])

    with open(batch_file, 'w') as f:
        # Header
        f.write('\t'.join(columns) + '\n')

        # Data rows
        for sample in samples:
            row = [sample['sample_id']]

            # R1 path (required)
            row.append(str(sample['r1_path']))

            # R2 path (optional)
            if has_r2:
                row.append(str(sample.get('r2_path', '')))

            # Amplicon sequence (per-sample or default)
            amplicon = sample.get('amplicon_seq', default_amplicon_seq) or ''
            row.append(amplicon)

            # Guide sequence (per-sample or default)
            guide = sample.get('guide_seq', default_guide_seq) or ''
            row.append(guide)

            # HDR sequence (per-sample or default)
            if default_hdr_seq or has_per_sample_hdr:
                hdr = sample.get('hdr_seq', default_hdr_seq) or ''
                row.append(hdr)

            # Quantification window and min identity
            row.append(str(quantification_window))
            row.append(str(min_identity))

            f.write('\t'.join(row) + '\n')

    logger.info(f"Wrote batch settings to {batch_file} ({len(samples)} samples)")


def _build_batch_command(
    runner: CRISPRessoRunner,
    batch_file: Path,
    output_dir: Path,
    threads: int,
    skip_finished: bool,
) -> list:
    """Build CRISPRessoBatch command."""
    if runner.execution_mode == 'singularity':
        cmd = ['singularity', 'exec', str(runner.container_path), 'CRISPRessoBatch']
    else:
        cmd = ['CRISPRessoBatch']

    cmd.extend([
        '--batch_settings', str(batch_file),
        '--output_folder', str(output_dir),
        '-p', str(threads),
    ])

    if skip_finished:
        cmd.append('--skip_finished')

    return cmd


# Legacy function for backwards compatibility
def run_crispresso_sequential(
    samples: list,
    amplicon_seq: str,
    guide_seq: str,
    output_dir: Path,
    hdr_seq: Optional[str] = None,
    threads: int = 4,
    container_path: Optional[Path] = None,
) -> Dict[str, CRISPRessoResult]:
    """
    Run CRISPResso2 on multiple samples sequentially (legacy method).

    Use run_crispresso_batch() instead for better performance.
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
