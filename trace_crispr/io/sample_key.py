"""
Sample key parsing and validation.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass
import pandas as pd
import logging

logger = logging.getLogger(__name__)


@dataclass
class Sample:
    """Represents a single sample.

    Attributes:
        sample_id: Unique sample identifier
        r1_path: Path to R1 FASTQ file
        r2_path: Optional path to R2 FASTQ file
        metadata: Additional metadata columns from sample key

        Per-sample locus overrides (optional):
        reference: Reference amplicon (DNA sequence or FASTA path)
        hdr_template: HDR template (DNA sequence or FASTA path)
        guide: Guide sequence (DNA sequence)

        If per-sample locus fields are not provided, the CLI defaults are used.
    """
    sample_id: str
    r1_path: Path
    r2_path: Optional[Path] = None
    metadata: Dict = None

    # Per-sample locus overrides (optional)
    reference: Optional[str] = None
    hdr_template: Optional[str] = None
    guide: Optional[str] = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

    def has_custom_locus(self) -> bool:
        """Check if sample has any per-sample locus overrides."""
        return any([self.reference, self.hdr_template, self.guide])

    def validate(self) -> List[str]:
        """Validate sample configuration. Returns list of errors."""
        errors = []

        if not self.r1_path.exists():
            errors.append(f"R1 file not found: {self.r1_path}")

        if self.r2_path and not self.r2_path.exists():
            errors.append(f"R2 file not found: {self.r2_path}")

        # Validate that if any locus field is set, all required fields are set
        locus_fields = [self.reference, self.hdr_template, self.guide]
        has_any = any(locus_fields)
        has_all = all(locus_fields)

        if has_any and not has_all:
            missing = []
            if not self.reference:
                missing.append('reference')
            if not self.hdr_template:
                missing.append('hdr_template')
            if not self.guide:
                missing.append('guide')
            errors.append(f"Partial locus override - missing: {', '.join(missing)}")

        return errors


def load_sample_key(
    path: Path,
    validate: bool = True
) -> List[Sample]:
    """
    Load samples from a sample key TSV file.

    Required columns:
    - sample_id: Unique sample identifier
    - r1_path: Path to R1 FASTQ file

    Optional columns:
    - r2_path: Path to R2 FASTQ file

    Per-sample locus override columns (optional):
    - reference: Reference amplicon (DNA sequence or FASTA path)
    - hdr_template: HDR template (DNA sequence or FASTA path)
    - guide: Guide sequence

    If a sample has any locus override column filled in, it must have all three.
    Empty values use the CLI-provided defaults.

    Additional columns are stored as metadata.

    Args:
        path: Path to sample key TSV
        validate: If True, validate that files exist

    Returns:
        List of Sample objects
    """
    df = pd.read_csv(path, sep='\t')

    # Validate required columns
    if 'sample_id' not in df.columns:
        raise ValueError("Sample key must have 'sample_id' column")

    if 'r1_path' not in df.columns:
        raise ValueError("Sample key must have 'r1_path' column")

    # Standard columns and locus override columns
    standard_cols = {'sample_id', 'r1_path', 'r2_path'}
    locus_cols = {'reference', 'hdr_template', 'guide'}

    samples = []
    errors = []

    for _, row in df.iterrows():
        sample_id = str(row['sample_id'])
        r1_path = Path(row['r1_path'])
        r2_path = Path(row['r2_path']) if pd.notna(row.get('r2_path', None)) else None

        # Parse per-sample locus overrides (if columns exist)
        reference = str(row['reference']).strip() if 'reference' in row and pd.notna(row.get('reference')) else None
        hdr_template = str(row['hdr_template']).strip() if 'hdr_template' in row and pd.notna(row.get('hdr_template')) else None
        guide = str(row['guide']).strip() if 'guide' in row and pd.notna(row.get('guide')) else None

        # Collect metadata from other columns
        metadata = {
            k: v for k, v in row.items()
            if k not in (standard_cols | locus_cols) and pd.notna(v)
        }

        sample = Sample(
            sample_id=sample_id,
            r1_path=r1_path,
            r2_path=r2_path,
            metadata=metadata,
            reference=reference,
            hdr_template=hdr_template,
            guide=guide,
        )

        if validate:
            sample_errors = sample.validate()
            for err in sample_errors:
                errors.append(f"{sample_id}: {err}")

        samples.append(sample)

    if errors:
        logger.warning(f"Sample key validation found {len(errors)} errors:")
        for err in errors[:10]:
            logger.warning(f"  {err}")
        if len(errors) > 10:
            logger.warning(f"  ... and {len(errors) - 10} more")

    # Log if per-sample loci detected
    samples_with_custom_locus = sum(1 for s in samples if s.has_custom_locus())
    if samples_with_custom_locus > 0:
        logger.info(f"Detected {samples_with_custom_locus} samples with per-sample locus overrides")

    return samples


def create_sample_key_template(output_path: Path, include_locus_columns: bool = False):
    """Create a template sample key file.

    Args:
        output_path: Path to write template file
        include_locus_columns: If True, include reference/hdr_template/guide columns
    """
    if include_locus_columns:
        template = """sample_id\tr1_path\tr2_path\tcondition\treplicate\treference\thdr_template\tguide
sample_1\t/path/to/sample_1_R1.fastq.gz\t/path/to/sample_1_R2.fastq.gz\ttreatment\t1\t\t\t
sample_2\t/path/to/sample_2_R1.fastq.gz\t/path/to/sample_2_R2.fastq.gz\tcontrol\t1\t\t\t
sample_3\t/path/to/sample_3_R1.fastq.gz\t/path/to/sample_3_R2.fastq.gz\ttreatment\t2\tlocus2.fasta\tlocus2_hdr.fasta\tACGTACGTACGTACGTACGT
sample_4\t/path/to/sample_4_R1.fastq.gz\t/path/to/sample_4_R2.fastq.gz\tcontrol\t2\tlocus2.fasta\tlocus2_hdr.fasta\tACGTACGTACGTACGTACGT
"""
    else:
        template = """sample_id\tr1_path\tr2_path\tcondition\treplicate
sample_1\t/path/to/sample_1_R1.fastq.gz\t/path/to/sample_1_R2.fastq.gz\ttreatment\t1
sample_2\t/path/to/sample_2_R1.fastq.gz\t/path/to/sample_2_R2.fastq.gz\tcontrol\t1
sample_3\t/path/to/sample_3_R1.fastq.gz\t/path/to/sample_3_R2.fastq.gz\ttreatment\t2
sample_4\t/path/to/sample_4_R1.fastq.gz\t/path/to/sample_4_R2.fastq.gz\tcontrol\t2
"""
    with open(output_path, 'w') as f:
        f.write(template)

    logger.info(f"Created sample key template: {output_path}")
