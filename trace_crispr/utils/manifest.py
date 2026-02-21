"""
Sample manifest generation utilities.

Generate TRACE-compatible sample manifests from project keyfiles.

Author: Kevin R. Roy
"""

import logging
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

logger = logging.getLogger(__name__)


def generate_trace_manifest(
    sample_key_path: Path,
    plate_key_path: Path,
    raw_data_dir: Path,
    output_path: Optional[Path] = None,
    locus_filter: Optional[str] = None,
    replicate_filter: Optional[str] = None,
) -> pd.DataFrame:
    """
    Generate TRACE-compatible sample manifest from project keyfiles.

    Joins sample_key with plate_key to resolve FASTQ paths.

    Args:
        sample_key_path: Path to sample_key.tsv
        plate_key_path: Path to plate_key.tsv
        raw_data_dir: Path to raw_data directory containing FASTQ files
        output_path: Optional path to write output TSV
        locus_filter: Optional locus to filter to (e.g., 'EMX1_exon_3_SCL819_SCL820')
        replicate_filter: Optional replicate to filter to (e.g., 'rep_1')

    Returns:
        DataFrame with TRACE-compatible sample manifest

    FASTQ naming convention (expected):
        {sequencing_date}_{outer_primers}_{sample_number}_R1.fastq.gz
        {sequencing_date}_{outer_primers}_{sample_number}_R2.fastq.gz
    """
    # Load keyfiles
    sample_key = pd.read_csv(sample_key_path, sep='\t')
    plate_key = pd.read_csv(plate_key_path, sep='\t')

    # Apply filters if specified
    if locus_filter:
        plate_key = plate_key[plate_key['locus'] == locus_filter]
        logger.info(f"Filtered to locus: {locus_filter} ({len(plate_key)} plate entries)")

    if replicate_filter:
        plate_key = plate_key[plate_key['replicate'] == replicate_filter]
        logger.info(f"Filtered to replicate: {replicate_filter} ({len(plate_key)} plate entries)")

    # Join sample_key with plate_key on gDNA_plate_num
    # This creates one row per (sample, plate/sequencing run) combination
    merged = sample_key.merge(
        plate_key,
        on='gDNA_plate_num',
        how='inner'
    )

    # Generate sample_id that includes sequencing context
    merged['sample_id'] = (
        merged['sequencing_date'].astype(str) + '_' +
        merged['outer_primers'] + '_' +
        merged['sample_number']
    )

    # Generate FASTQ paths
    fastq_dir = Path(raw_data_dir) / 'fastq'

    merged['r1_path'] = merged.apply(
        lambda row: str(fastq_dir / f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_R1.fastq.gz"),
        axis=1
    )
    merged['r2_path'] = merged.apply(
        lambda row: str(fastq_dir / f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_R2.fastq.gz"),
        axis=1
    )

    # Rename barcode column for TRACE compatibility
    merged = merged.rename(columns={'barcode': 'expected_barcode'})

    # Select and order columns for TRACE
    trace_columns = [
        'sample_id',
        'r1_path',
        'r2_path',
        'expected_barcode',
        'sample_number',
        'gDNA_plate_num',
        'Well ID',
        'sample_ID',  # Original detailed sample description
        'date',
        'locus',
        'sequencing_date',
        'outer_primers',
        'replicate',
    ]

    # Only include columns that exist
    available_columns = [c for c in trace_columns if c in merged.columns]
    result = merged[available_columns].copy()

    # Validate FASTQ files exist
    missing_files = []
    for _, row in result.iterrows():
        r1 = Path(row['r1_path'])
        r2 = Path(row['r2_path'])
        if not r1.exists():
            missing_files.append(str(r1))
        if not r2.exists():
            missing_files.append(str(r2))

    if missing_files:
        logger.warning(f"Missing {len(missing_files)} FASTQ files:")
        for f in missing_files[:10]:
            logger.warning(f"  {f}")
        if len(missing_files) > 10:
            logger.warning(f"  ... and {len(missing_files) - 10} more")

    # Write output if path specified
    if output_path:
        result.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Wrote {len(result)} samples to {output_path}")

    return result


def extract_barcodes_from_sample_key(
    sample_key_path: Path,
    barcode_column: str = 'barcode',
) -> Dict[str, str]:
    """
    Extract unique barcodes from sample_key.

    Args:
        sample_key_path: Path to sample_key.tsv
        barcode_column: Column name containing barcodes

    Returns:
        Dict mapping barcode_id to barcode sequence
        (barcode_id is the barcode sequence itself for uniqueness)
    """
    sample_key = pd.read_csv(sample_key_path, sep='\t')

    # Filter out NA/empty barcodes
    valid_barcodes = sample_key[
        sample_key[barcode_column].notna() &
        (sample_key[barcode_column] != 'NA') &
        (sample_key[barcode_column] != '')
    ][barcode_column].unique()

    # Create dict mapping barcode sequence to itself
    # (could also derive node-based names from sample_ID column)
    barcodes = {str(bc): str(bc) for bc in valid_barcodes}

    logger.info(f"Extracted {len(barcodes)} unique barcodes")

    return barcodes


def get_sample_expected_barcode(
    manifest: pd.DataFrame,
    sample_id: str,
    expected_barcode_column: str = 'expected_barcode',
) -> Optional[str]:
    """
    Get the expected barcode for a sample from the manifest.

    Args:
        manifest: TRACE manifest DataFrame
        sample_id: Sample ID to look up
        expected_barcode_column: Column containing expected barcodes

    Returns:
        Expected barcode sequence, or None if not found/NA
    """
    matches = manifest[manifest['sample_id'] == sample_id]
    if len(matches) == 0:
        return None

    barcode = matches[expected_barcode_column].iloc[0]
    if pd.isna(barcode) or barcode == 'NA' or barcode == '':
        return None

    return str(barcode)
