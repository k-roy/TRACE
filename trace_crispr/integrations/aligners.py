"""
Aligner wrappers for BWA-MEM, BBMap, and minimap2.

Author: Kevin R. Roy
"""

import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class AlignerResult:
    """Result from running an aligner."""
    aligner: str
    bam_path: Path
    success: bool
    error_message: Optional[str] = None
    reads_aligned: int = 0
    reads_total: int = 0


class AlignerManager:
    """Manage external aligner availability and execution."""

    REQUIRED_ALIGNERS = ['bwa', 'bbmap', 'minimap2']
    OPTIONAL_TOOLS = ['samtools']

    def __init__(self):
        self.available = {}
        self._detect_aligners()

    def _detect_aligners(self):
        """Check which aligners are available."""
        # Check BWA
        try:
            result = subprocess.run(['bwa'], capture_output=True, timeout=5)
            self.available['bwa'] = True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.available['bwa'] = False

        # Check BBMap
        try:
            result = subprocess.run(['bbmap.sh', '--version'], capture_output=True, timeout=5)
            self.available['bbmap'] = True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.available['bbmap'] = False

        # Check minimap2
        try:
            result = subprocess.run(['minimap2', '--version'], capture_output=True, timeout=5)
            self.available['minimap2'] = True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.available['minimap2'] = False

        # Check samtools
        try:
            result = subprocess.run(['samtools', '--version'], capture_output=True, timeout=5)
            self.available['samtools'] = True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.available['samtools'] = False

    def check_requirements(self) -> List[str]:
        """Return list of missing required aligners."""
        missing = [a for a in self.REQUIRED_ALIGNERS if not self.available.get(a)]
        return missing

    def get_installation_instructions(self, missing: List[str]) -> str:
        """Return installation instructions for missing tools."""
        return f"""
Missing aligners: {', '.join(missing)}

Install via conda:
    conda install -c bioconda bwa bbmap minimap2 samtools

Or via mamba (faster):
    mamba install -c bioconda bwa bbmap minimap2 samtools
"""


def run_bwa_mem(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    reference: Path,
    output_bam: Path,
    threads: int = 4,
) -> AlignerResult:
    """
    Run BWA-MEM alignment.

    Args:
        r1_fastq: Path to R1 FASTQ
        r2_fastq: Optional path to R2 FASTQ
        reference: Path to reference FASTA
        output_bam: Path for output BAM file
        threads: Number of threads

    Returns:
        AlignerResult with alignment info
    """
    try:
        # Check if index exists, create if not
        index_files = [reference.with_suffix(reference.suffix + ext)
                      for ext in ['.bwt', '.pac', '.ann', '.amb', '.sa']]
        if not all(f.exists() for f in index_files):
            logger.info(f"Creating BWA index for {reference}")
            subprocess.run(
                ['bwa', 'index', str(reference)],
                check=True, capture_output=True
            )

        # Build alignment command
        cmd = ['bwa', 'mem', '-t', str(threads), str(reference), str(r1_fastq)]
        if r2_fastq:
            cmd.append(str(r2_fastq))

        # Run alignment and convert to BAM
        with tempfile.NamedTemporaryFile(suffix='.sam', delete=False) as sam_file:
            sam_path = Path(sam_file.name)

        # Run BWA
        with open(sam_path, 'w') as sam_out:
            result = subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, check=True)

        # Convert to BAM and sort
        subprocess.run(
            ['samtools', 'view', '-bS', str(sam_path), '-o', str(output_bam)],
            check=True, capture_output=True
        )

        # Clean up
        sam_path.unlink()

        return AlignerResult(
            aligner='bwa',
            bam_path=output_bam,
            success=True,
        )

    except subprocess.CalledProcessError as e:
        return AlignerResult(
            aligner='bwa',
            bam_path=output_bam,
            success=False,
            error_message=e.stderr.decode() if e.stderr else str(e),
        )
    except Exception as e:
        return AlignerResult(
            aligner='bwa',
            bam_path=output_bam,
            success=False,
            error_message=str(e),
        )


def run_bbmap(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    reference: Path,
    output_bam: Path,
    threads: int = 4,
    max_indel: int = 600,
) -> AlignerResult:
    """
    Run BBMap alignment.

    Args:
        r1_fastq: Path to R1 FASTQ
        r2_fastq: Optional path to R2 FASTQ
        reference: Path to reference FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        max_indel: Maximum indel size to allow

    Returns:
        AlignerResult with alignment info
    """
    try:
        # Build command
        cmd = [
            'bbmap.sh',
            f'ref={reference}',
            f'in={r1_fastq}',
            f'out={output_bam}',
            f'threads={threads}',
            f'maxindel={max_indel}',
            'sam=1.3',
            'nodisk=t',
            'local=t',
        ]

        if r2_fastq:
            cmd.append(f'in2={r2_fastq}')

        # Run BBMap
        result = subprocess.run(cmd, capture_output=True, check=True)

        return AlignerResult(
            aligner='bbmap',
            bam_path=output_bam,
            success=True,
        )

    except subprocess.CalledProcessError as e:
        return AlignerResult(
            aligner='bbmap',
            bam_path=output_bam,
            success=False,
            error_message=e.stderr.decode() if e.stderr else str(e),
        )
    except Exception as e:
        return AlignerResult(
            aligner='bbmap',
            bam_path=output_bam,
            success=False,
            error_message=str(e),
        )


def run_minimap2(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    reference: Path,
    output_bam: Path,
    threads: int = 4,
    end_bonus: int = 50,
) -> AlignerResult:
    """
    Run minimap2 alignment.

    Args:
        r1_fastq: Path to R1 FASTQ
        r2_fastq: Optional path to R2 FASTQ
        reference: Path to reference FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        end_bonus: Bonus for extending alignments to ends

    Returns:
        AlignerResult with alignment info
    """
    try:
        # Build command
        cmd = [
            'minimap2',
            '-ax', 'sr',  # Short reads preset
            '-t', str(threads),
            f'--end-bonus={end_bonus}',
            str(reference),
            str(r1_fastq),
        ]

        if r2_fastq:
            cmd.append(str(r2_fastq))

        # Run minimap2 and convert to BAM
        with tempfile.NamedTemporaryFile(suffix='.sam', delete=False) as sam_file:
            sam_path = Path(sam_file.name)

        with open(sam_path, 'w') as sam_out:
            result = subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, check=True)

        # Convert to BAM
        subprocess.run(
            ['samtools', 'view', '-bS', str(sam_path), '-o', str(output_bam)],
            check=True, capture_output=True
        )

        # Clean up
        sam_path.unlink()

        return AlignerResult(
            aligner='minimap2',
            bam_path=output_bam,
            success=True,
        )

    except subprocess.CalledProcessError as e:
        return AlignerResult(
            aligner='minimap2',
            bam_path=output_bam,
            success=False,
            error_message=e.stderr.decode() if e.stderr else str(e),
        )
    except Exception as e:
        return AlignerResult(
            aligner='minimap2',
            bam_path=output_bam,
            success=False,
            error_message=str(e),
        )


def run_triple_alignment(
    r1_fastq: Path,
    r2_fastq: Optional[Path],
    reference: Path,
    output_dir: Path,
    threads: int = 4,
    aligners: List[str] = None,
) -> Dict[str, AlignerResult]:
    """
    Run all three aligners and return results.

    Args:
        r1_fastq: Path to R1 FASTQ
        r2_fastq: Optional path to R2 FASTQ
        reference: Path to reference FASTA
        output_dir: Directory for output files
        threads: Number of threads per aligner
        aligners: List of aligners to run (default: all three)

    Returns:
        Dict mapping aligner name to AlignerResult
    """
    if aligners is None:
        aligners = ['bwa', 'bbmap', 'minimap2']

    output_dir.mkdir(parents=True, exist_ok=True)
    results = {}

    aligner_funcs = {
        'bwa': run_bwa_mem,
        'bbmap': run_bbmap,
        'minimap2': run_minimap2,
    }

    for aligner in aligners:
        if aligner not in aligner_funcs:
            logger.warning(f"Unknown aligner: {aligner}")
            continue

        output_bam = output_dir / f"{aligner}.bam"
        logger.info(f"Running {aligner}...")

        result = aligner_funcs[aligner](
            r1_fastq=r1_fastq,
            r2_fastq=r2_fastq,
            reference=reference,
            output_bam=output_bam,
            threads=threads,
        )

        results[aligner] = result

        if result.success:
            logger.info(f"{aligner} completed successfully")
        else:
            logger.warning(f"{aligner} failed: {result.error_message}")

    return results


def create_reference_fasta(sequence: str, output_path: Path, name: str = "reference") -> Path:
    """
    Create a FASTA file from a sequence string.

    Args:
        sequence: DNA sequence
        output_path: Path for output FASTA
        name: Sequence name for FASTA header

    Returns:
        Path to created FASTA file
    """
    with open(output_path, 'w') as f:
        f.write(f">{name}\n")
        # Write sequence in 80-character lines
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")

    return output_path
