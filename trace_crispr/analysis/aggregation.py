"""
Data aggregation functions for TRACE analysis.

Aggregates SampleResult data by condition, batch, or other grouping variables
for statistical comparison and visualization.

Author: Kevin R. Roy
"""

from typing import List, Dict, Optional, Union, Callable
import pandas as pd
import numpy as np

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
            "dedup_wt_count": r.dedup_wt_count,
            "dedup_hdr_count": r.dedup_hdr_count,
            "dedup_nhej_count": r.dedup_nhej_count,
            "dedup_large_del_count": r.dedup_large_del_count,
            "dedup_total": r.dedup_total,
            "dedup_wt_pct": r.dedup_wt_pct,
            "dedup_hdr_pct": r.dedup_hdr_pct,
            "dedup_nhej_pct": r.dedup_nhej_pct,
            "dedup_large_del_pct": r.dedup_large_del_pct,
            "kmer_hdr_rate": r.kmer_hdr_rate,
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
