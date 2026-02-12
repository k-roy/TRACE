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
    """Represents a single sample."""
    sample_id: str
    r1_path: Path
    r2_path: Optional[Path] = None
    metadata: Dict = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

    def validate(self) -> List[str]:
        """Validate sample configuration. Returns list of errors."""
        errors = []

        if not self.r1_path.exists():
            errors.append(f"R1 file not found: {self.r1_path}")

        if self.r2_path and not self.r2_path.exists():
            errors.append(f"R2 file not found: {self.r2_path}")

        return errors


def load_sample_key(
    path: Path,
    validate: bool = True
) -> List[Sample]:
    """
    Load samples from a sample key TSV file.

    Expected columns:
    - sample_id: Unique sample identifier
    - r1_path: Path to R1 FASTQ file
    - r2_path: (Optional) Path to R2 FASTQ file

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

    samples = []
    errors = []

    for _, row in df.iterrows():
        sample_id = str(row['sample_id'])
        r1_path = Path(row['r1_path'])
        r2_path = Path(row['r2_path']) if pd.notna(row.get('r2_path', None)) else None

        # Collect metadata from other columns
        metadata = {
            k: v for k, v in row.items()
            if k not in ('sample_id', 'r1_path', 'r2_path') and pd.notna(v)
        }

        sample = Sample(
            sample_id=sample_id,
            r1_path=r1_path,
            r2_path=r2_path,
            metadata=metadata,
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

    return samples


def create_sample_key_template(output_path: Path):
    """Create a template sample key file."""
    template = """sample_id\tr1_path\tr2_path\tcondition\treplicate
sample_1\t/path/to/sample_1_R1.fastq.gz\t/path/to/sample_1_R2.fastq.gz\ttreatment\t1
sample_2\t/path/to/sample_2_R1.fastq.gz\t/path/to/sample_2_R2.fastq.gz\tcontrol\t1
sample_3\t/path/to/sample_3_R1.fastq.gz\t/path/to/sample_3_R2.fastq.gz\ttreatment\t2
sample_4\t/path/to/sample_4_R1.fastq.gz\t/path/to/sample_4_R2.fastq.gz\tcontrol\t2
"""
    with open(output_path, 'w') as f:
        f.write(template)

    logger.info(f"Created sample key template: {output_path}")
