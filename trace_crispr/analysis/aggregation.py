"""
Data aggregation functions for TRACE analysis.

Aggregates SampleResult data by condition, batch, or other grouping variables
for statistical comparison and visualization.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Callable, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from ..io.output import SampleResult
from ..io.sample_key import Sample
from .types import ConditionStats


def results_to_dataframe(
    results: List[SampleResult],
    samples: Optional[List[Sample]] = None,
) -> pd.DataFrame:
    """
    Convert list of SampleResult objects to a pandas DataFrame.

    Merges SampleResult metrics with Sample metadata if provided.
    Includes all new comprehensive classification columns.

    Args:
        results: List of SampleResult objects from TRACE pipeline
        samples: Optional list of Sample objects for metadata merging

    Returns:
        DataFrame with one row per sample, columns for all metrics and metadata

    Example:
        >>> results = pipeline.run_all(samples)
        >>> df = results_to_dataframe(results, samples)
        >>> df.groupby('condition')['dedup_hdr_pct'].mean()
    """
    rows = []
    for r in results:
        row = {
            "sample_id": r.sample_id,
            "total_reads": r.total_reads,
            "aligned_reads": r.aligned_reads,
            "classifiable_reads": r.classifiable_reads,
            "duplicate_rate": r.duplicate_rate,

            # === New comprehensive classification counts ===
            # HDR categories
            "dedup_hdr_complete_count": r.dedup_hdr_complete_count,
            "dedup_hdr_partial_count": r.dedup_hdr_partial_count,
            "dedup_hdr_plus_nhej_count": r.dedup_hdr_plus_nhej_count,
            "dedup_hdr_plus_mmej_count": r.dedup_hdr_plus_mmej_count,
            "dedup_hdr_plus_other_count": r.dedup_hdr_plus_other_count,
            "dedup_donor_capture_count": r.dedup_donor_capture_count,

            # NHEJ/MMEJ categories
            "dedup_nhej_indel_count": r.dedup_nhej_indel_count,
            "dedup_mmej_indel_count": r.dedup_mmej_indel_count,

            # Other categories
            "dedup_wt_count": r.dedup_wt_count,
            "dedup_non_donor_snv_count": r.dedup_non_donor_snv_count,
            "dedup_unclassified_count": r.dedup_unclassified_count,
            "dedup_unmapped_count": r.dedup_unmapped_count,

            # === Aggregated counts ===
            "dedup_hdr_total": r.dedup_hdr_total,
            "dedup_nhej_mmej_total": r.dedup_nhej_mmej_total,
            "dedup_edited_total": r.dedup_edited_total,
            "dedup_total": r.dedup_total,

            # === New comprehensive classification percentages ===
            # HDR categories
            "dedup_hdr_complete_pct": r.dedup_hdr_complete_pct,
            "dedup_hdr_partial_pct": r.dedup_hdr_partial_pct,
            "dedup_hdr_plus_nhej_pct": r.dedup_hdr_plus_nhej_pct,
            "dedup_hdr_plus_mmej_pct": r.dedup_hdr_plus_mmej_pct,
            "dedup_hdr_plus_other_pct": r.dedup_hdr_plus_other_pct,
            "dedup_donor_capture_pct": r.dedup_donor_capture_pct,

            # NHEJ/MMEJ categories
            "dedup_nhej_indel_pct": r.dedup_nhej_indel_pct,
            "dedup_mmej_indel_pct": r.dedup_mmej_indel_pct,

            # Other categories
            "dedup_wt_pct": r.dedup_wt_pct,
            "dedup_non_donor_snv_pct": r.dedup_non_donor_snv_pct,
            "dedup_unclassified_pct": r.dedup_unclassified_pct,

            # Aggregated percentages
            "dedup_hdr_total_pct": r.dedup_hdr_total_pct,
            "dedup_nhej_mmej_total_pct": r.dedup_nhej_mmej_total_pct,
            "dedup_edited_total_pct": r.dedup_edited_total_pct,

            # === Legacy fields (backwards compatibility) ===
            "dedup_hdr_count": r.dedup_hdr_count,
            "dedup_nhej_count": r.dedup_nhej_count,
            "dedup_large_del_count": r.dedup_large_del_count,
            "dedup_hdr_pct": r.dedup_hdr_pct,
            "dedup_nhej_pct": r.dedup_nhej_pct,
            "dedup_large_del_pct": r.dedup_large_del_pct,

            # K-mer method
            "kmer_hdr_rate": r.kmer_hdr_rate,

            # CRISPResso comparison
            "crispresso_hdr_rate": r.crispresso_hdr_rate,
            "crispresso_nhej_rate": r.crispresso_nhej_rate,
        }

        # Add SampleResult metadata
        if r.metadata:
            for k, v in r.metadata.items():
                row[k] = v

        rows.append(row)

    df = pd.DataFrame(rows)

    # Merge sample metadata if provided
    if samples:
        sample_meta = {}
        for s in samples:
            sample_meta[s.sample_id] = s.metadata

        # Add sample metadata columns that don't already exist
        for idx, row in df.iterrows():
            sample_id = row["sample_id"]
            if sample_id in sample_meta:
                for key, value in sample_meta[sample_id].items():
                    if key not in df.columns:
                        df[key] = None
                    df.at[idx, key] = value

    return df


def get_condition_stats(
    df: pd.DataFrame,
    condition_col: str,
    metric: str,
    condition_value: Optional[str] = None,
    filter_func: Optional[Callable[[pd.DataFrame], pd.Series]] = None,
) -> Union[ConditionStats, Dict[str, ConditionStats]]:
    """
    Calculate statistics for a metric grouped by condition.

    Args:
        df: DataFrame with sample data
        condition_col: Column name containing condition labels
        metric: Column name for the metric to analyze
        condition_value: If specified, return stats for only this condition
        filter_func: Optional function to filter rows before aggregation

    Returns:
        Single ConditionStats if condition_value specified,
        otherwise Dict mapping condition names to ConditionStats

    Example:
        >>> stats = get_condition_stats(df, 'condition', 'dedup_hdr_pct')
        >>> print(f"Control HDR: {stats['control'].mean:.1f}%")

        >>> treatment_stats = get_condition_stats(
        ...     df, 'condition', 'dedup_hdr_pct',
        ...     condition_value='treatment'
        ... )
    """
    # Validate columns exist
    if condition_col not in df.columns:
        raise ValueError(
            f"Condition column '{condition_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    if metric not in df.columns:
        raise ValueError(
            f"Metric '{metric}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    # Apply filter if provided
    if filter_func is not None:
        mask = filter_func(df)
        df = df[mask].copy()

    if condition_value is not None:
        subset = df[df[condition_col] == condition_value]
        return ConditionStats.from_series(
            condition=condition_value,
            metric=metric,
            data=subset[metric],
            sample_ids=(
                subset["sample_id"].tolist() if "sample_id" in subset.columns else []
            ),
        )

    # Return stats for all conditions
    result = {}
    for condition in df[condition_col].unique():
        subset = df[df[condition_col] == condition]
        result[condition] = ConditionStats.from_series(
            condition=str(condition),
            metric=metric,
            data=subset[metric],
            sample_ids=(
                subset["sample_id"].tolist() if "sample_id" in subset.columns else []
            ),
        )

    return result


def aggregate_by_group(
    df: pd.DataFrame,
    group_cols: Union[str, List[str]],
    metrics: Union[str, List[str]],
    agg_funcs: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Aggregate metrics by one or more grouping columns.

    Args:
        df: DataFrame with sample data
        group_cols: Column(s) to group by
        metrics: Metric column(s) to aggregate
        agg_funcs: Aggregation functions (default: ['mean', 'std', 'sem', 'count'])

    Returns:
        DataFrame with aggregated statistics

    Example:
        >>> summary = aggregate_by_group(
        ...     df,
        ...     group_cols=['condition', 'date'],
        ...     metrics=['dedup_hdr_pct', 'dedup_nhej_pct']
        ... )
    """
    if isinstance(group_cols, str):
        group_cols = [group_cols]
    if isinstance(metrics, str):
        metrics = [metrics]
    if agg_funcs is None:
        agg_funcs = ["mean", "std", "sem", "count"]

    def sem(x):
        return x.std() / np.sqrt(len(x)) if len(x) > 1 else 0

    agg_dict = {}
    for metric in metrics:
        for func in agg_funcs:
            if func == "sem":
                agg_dict[f"{metric}_{func}"] = (metric, sem)
            else:
                agg_dict[f"{metric}_{func}"] = (metric, func)

    grouped = df.groupby(group_cols, as_index=False).agg(**agg_dict)

    return grouped


# ============================================================================
# Unified Comparison Table Generation
# ============================================================================

# QC thresholds
MIN_READS_THRESHOLD = 1000


def _assign_quality_flag(
    replicate_reads: List[int],
    min_reads: int = MIN_READS_THRESHOLD,
) -> tuple:
    """
    Assign quality flag based on replicate read counts.

    Args:
        replicate_reads: List of read counts for each replicate
        min_reads: Minimum reads threshold for a replicate to be "good"

    Returns:
        Tuple of (quality_flag, n_replicates_used, valid_replicate_indices)

    Quality flags:
        - good: All replicates have >= min_reads
        - all_low_reads: All replicates have < min_reads (still used in mean)
        - one_low_read_removed: Low-read replicates excluded from mean
        - no_data: No data available
    """
    if not replicate_reads:
        return "no_data", 0, []

    valid_indices = [i for i, reads in enumerate(replicate_reads) if reads >= min_reads]
    low_count = sum(1 for reads in replicate_reads if reads < min_reads)
    total = len(replicate_reads)

    if low_count == 0:
        # All replicates good
        return "good", total, list(range(total))
    elif low_count == total:
        # All replicates low - still use them but flag
        return "all_low_reads", total, list(range(total))
    else:
        # Some replicates low - exclude them
        return "one_low_read_removed", len(valid_indices), valid_indices


def generate_unified_comparison_table(
    df: pd.DataFrame,
    bio_sample_col: str = "bio_sample",
    min_reads: int = MIN_READS_THRESHOLD,
    include_crispresso: bool = True,
    metadata_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Generate unified comparison table with replicate aggregation and QC flags.

    Groups samples by bio_sample identifier, calculates mean/SEM for all metrics,
    applies QC flags based on read counts, and optionally includes CRISPResso comparison.

    Column Naming Convention:
        - TRACE output columns have no prefix (e.g., hdr_complete_pct_mean)
        - CRISPResso columns are prefixed with "crispresso_" (e.g., crispresso_hdr_rate_mean)

    Aggregated Metric Calculations:
        hdr_total_pct = (
            hdr_complete_pct +
            hdr_partial_pct +
            hdr_plus_nhej_pct +
            hdr_plus_mmej_pct +
            hdr_plus_other_pct
        )

        nhej_mmej_total_pct = nhej_indel_pct + mmej_indel_pct

        edited_total_pct = hdr_total_pct + nhej_mmej_total_pct + donor_capture_pct

    Note on DONOR_CAPTURE:
        DONOR_CAPTURE is NOT included in hdr_total_pct because it represents
        over-integration of donor material, which is a distinct biological outcome
        from successful HDR. DONOR_CAPTURE IS included in edited_total_pct.

    Note on CRISPResso Comparison:
        CRISPResso does not distinguish DONOR_CAPTURE from HDR. CRISPResso's HDR
        rate includes reads that TRACE classifies as DONOR_CAPTURE. For direct
        comparison, use: trace_hdr_total + trace_donor_capture ≈ crispresso_hdr.

    Args:
        df: DataFrame with per-sample results (from results_to_dataframe)
        bio_sample_col: Column name containing biological sample identifier
        min_reads: Minimum reads threshold for QC
        include_crispresso: Include CRISPResso comparison columns
        metadata_cols: Optional list of metadata columns to include

    Returns:
        DataFrame with one row per bio_sample, columns for:
        - Sample identifiers and metadata
        - QC flags (quality_flag, n_replicates_used)
        - Mean/SEM for all classification metrics
        - CRISPResso comparison (if available)

    Example:
        >>> df = results_to_dataframe(results, samples)
        >>> unified = generate_unified_comparison_table(df, bio_sample_col='bio_sample')
        >>> unified.to_csv('unified_comparison_table.tsv', sep='\\t', index=False)
    """
    if bio_sample_col not in df.columns:
        raise ValueError(
            f"bio_sample column '{bio_sample_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    # Define metrics to aggregate
    # Maps internal field name -> output column name (without "dedup_" prefix)
    pct_metrics = {
        # New comprehensive classification percentages
        "dedup_hdr_complete_pct": "hdr_complete_pct",
        "dedup_hdr_partial_pct": "hdr_partial_pct",
        "dedup_hdr_plus_nhej_pct": "hdr_plus_nhej_pct",
        "dedup_hdr_plus_mmej_pct": "hdr_plus_mmej_pct",
        "dedup_hdr_plus_other_pct": "hdr_plus_other_pct",
        "dedup_donor_capture_pct": "donor_capture_pct",
        "dedup_nhej_indel_pct": "nhej_indel_pct",
        "dedup_mmej_indel_pct": "mmej_indel_pct",
        "dedup_wt_pct": "wt_pct",
        "dedup_non_donor_snv_pct": "non_donor_snv_pct",
        "dedup_unclassified_pct": "unclassified_pct",
        # Aggregated percentages
        "dedup_hdr_total_pct": "hdr_total_pct",
        "dedup_nhej_mmej_total_pct": "nhej_mmej_total_pct",
        "dedup_edited_total_pct": "edited_total_pct",
    }

    count_metrics = {
        "total_reads": "total_reads",
        "aligned_reads": "aligned_reads",
        "classifiable_reads": "classifiable_reads",
        "dedup_total": "trace_total_reads",
    }

    # CRISPResso metrics (if available)
    crispresso_metrics = ["crispresso_hdr_rate", "crispresso_nhej_rate"]

    rows = []
    grouped = df.groupby(bio_sample_col)

    for bio_sample, group in grouped:
        row = {bio_sample_col: bio_sample}

        # Get replicate read counts for QC
        if "dedup_total" in group.columns:
            replicate_reads = group["dedup_total"].tolist()
        elif "total_reads" in group.columns:
            replicate_reads = group["total_reads"].tolist()
        else:
            replicate_reads = []

        # Assign QC flag
        quality_flag, n_replicates_used, valid_indices = _assign_quality_flag(
            replicate_reads, min_reads
        )

        row["quality_flag"] = quality_flag
        row["n_replicates"] = len(group)
        row["n_replicates_used"] = n_replicates_used

        # Filter to valid replicates for mean calculation
        if valid_indices and quality_flag != "all_low_reads":
            valid_group = group.iloc[valid_indices]
        else:
            # Use all replicates if all_low_reads or no filtering needed
            valid_group = group

        # Aggregate count metrics (sum)
        for internal_name, output_name in count_metrics.items():
            if internal_name in valid_group.columns:
                row[f"{output_name}_total"] = valid_group[internal_name].sum()

        # Aggregate percentage metrics (mean and SEM)
        for internal_name, output_name in pct_metrics.items():
            if internal_name in valid_group.columns:
                values = valid_group[internal_name].dropna()
                if len(values) > 0:
                    row[f"{output_name}_mean"] = values.mean()
                    row[f"{output_name}_sem"] = values.std() / np.sqrt(len(values)) if len(values) > 1 else 0
                else:
                    row[f"{output_name}_mean"] = np.nan
                    row[f"{output_name}_sem"] = np.nan

        # Aggregate CRISPResso metrics (if available and requested)
        if include_crispresso:
            for metric in crispresso_metrics:
                if metric in valid_group.columns:
                    values = valid_group[metric].dropna()
                    if len(values) > 0:
                        row[f"{metric}_mean"] = values.mean()
                        row[f"{metric}_sem"] = values.std() / np.sqrt(len(values)) if len(values) > 1 else 0
                    else:
                        row[f"{metric}_mean"] = np.nan
                        row[f"{metric}_sem"] = np.nan

        # Add metadata columns (take first value from group)
        if metadata_cols:
            for col in metadata_cols:
                if col in group.columns:
                    row[col] = group[col].iloc[0]

        # Auto-detect and add common metadata columns
        common_metadata = ["plate", "well", "experiment_date", "category", "sample_description"]
        for col in common_metadata:
            if col in group.columns and col not in row:
                row[col] = group[col].iloc[0]

        rows.append(row)

    result_df = pd.DataFrame(rows)

    # Reorder columns for readability
    column_order = _get_unified_table_column_order(result_df, bio_sample_col, include_crispresso)
    existing_cols = [c for c in column_order if c in result_df.columns]
    remaining_cols = [c for c in result_df.columns if c not in existing_cols]
    result_df = result_df[existing_cols + remaining_cols]

    return result_df


def _get_unified_table_column_order(
    df: pd.DataFrame,
    bio_sample_col: str,
    include_crispresso: bool,
) -> List[str]:
    """
    Get preferred column order for unified comparison table.

    Column naming convention:
    - TRACE output columns have no prefix (e.g., hdr_complete_pct_mean)
    - CRISPResso columns are prefixed with "crispresso_" (e.g., crispresso_hdr_rate_mean)

    HDR Total calculation:
        hdr_total = hdr_complete + hdr_partial + hdr_plus_nhej + hdr_plus_mmej + hdr_plus_other

    Note: DONOR_CAPTURE is NOT included in hdr_total because it represents
    over-integration of donor material (a distinct biological outcome).
    """
    order = [
        # Identifiers
        bio_sample_col,
        "sample_description",
        "plate",
        "well",
        "experiment_date",
        "category",

        # QC
        "quality_flag",
        "n_replicates",
        "n_replicates_used",

        # Read counts
        "total_reads_total",
        "aligned_reads_total",
        "classifiable_reads_total",
        "trace_total_reads_total",

        # HDR categories (mean and SEM)
        "hdr_complete_pct_mean",
        "hdr_complete_pct_sem",
        "hdr_partial_pct_mean",
        "hdr_partial_pct_sem",
        "hdr_plus_nhej_pct_mean",
        "hdr_plus_nhej_pct_sem",
        "hdr_plus_mmej_pct_mean",
        "hdr_plus_mmej_pct_sem",
        "hdr_plus_other_pct_mean",
        "hdr_plus_other_pct_sem",

        # HDR total (= hdr_complete + hdr_partial + hdr_plus_nhej + hdr_plus_mmej + hdr_plus_other)
        "hdr_total_pct_mean",
        "hdr_total_pct_sem",

        # Donor capture (separate from HDR total)
        "donor_capture_pct_mean",
        "donor_capture_pct_sem",

        # NHEJ/MMEJ categories
        "nhej_indel_pct_mean",
        "nhej_indel_pct_sem",
        "mmej_indel_pct_mean",
        "mmej_indel_pct_sem",
        "nhej_mmej_total_pct_mean",
        "nhej_mmej_total_pct_sem",

        # Other categories
        "wt_pct_mean",
        "wt_pct_sem",
        "non_donor_snv_pct_mean",
        "non_donor_snv_pct_sem",
        "unclassified_pct_mean",
        "unclassified_pct_sem",

        # Edited total (= hdr_total + nhej_mmej_total + donor_capture)
        "edited_total_pct_mean",
        "edited_total_pct_sem",
    ]

    # Add CRISPResso columns if requested
    if include_crispresso:
        order.extend([
            "crispresso_hdr_rate_mean",
            "crispresso_hdr_rate_sem",
            "crispresso_nhej_rate_mean",
            "crispresso_nhej_rate_sem",
        ])

    return order


def write_unified_comparison_table(
    df: pd.DataFrame,
    output_path: Path,
    float_format: str = "%.4f",
) -> Path:
    """
    Write unified comparison table to TSV file.

    Args:
        df: DataFrame from generate_unified_comparison_table()
        output_path: Path for output TSV file
        float_format: Format string for floating point numbers

    Returns:
        Path to written file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_path, sep='\t', index=False, float_format=float_format)

    return output_path
