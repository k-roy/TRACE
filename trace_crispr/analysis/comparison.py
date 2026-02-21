"""
Condition comparison workflows for TRACE.

Provides high-level functions for comparing editing outcomes across conditions.

Author: Kevin R. Roy
"""

from typing import List, Optional

import pandas as pd

from ..io.output import SampleResult
from ..io.sample_key import Sample
from .aggregation import get_condition_stats, results_to_dataframe
from .statistics import compare_conditions
from .types import ComparisonSet


def compare_metric_by_condition(
    results: List[SampleResult],
    samples: Optional[List[Sample]] = None,
    condition_col: str = "condition",
    metric: str = "dedup_hdr_pct",
    base_condition: Optional[str] = None,
    fdr_correction: bool = True,
    min_replicates: int = 2,
) -> ComparisonSet:
    """
    Compare a metric across conditions with statistical testing.

    This is the main entry point for condition comparison analysis.

    Args:
        results: List of SampleResult objects from TRACE pipeline
        samples: Optional list of Sample objects for metadata
        condition_col: Column name containing condition labels
        metric: Metric to compare (e.g., 'dedup_hdr_pct', 'dedup_nhej_pct')
        base_condition: Reference condition for comparison. If None, uses
                        the first condition alphabetically.
        fdr_correction: Apply Benjamini-Hochberg FDR correction
        min_replicates: Minimum replicates required for statistical test

    Returns:
        ComparisonSet with statistical results for all comparisons

    Example:
        >>> from trace_crispr.analysis import compare_metric_by_condition
        >>>
        >>> results = pipeline.run_all(samples)
        >>> comparisons = compare_metric_by_condition(
        ...     results, samples,
        ...     condition_col='treatment',
        ...     metric='dedup_hdr_pct',
        ...     base_condition='untreated'
        ... )
        >>>
        >>> print(comparisons.to_dataframe())
        >>> print(f"Significant conditions: {comparisons.significant_conditions}")
    """
    # Convert to DataFrame
    df = results_to_dataframe(results, samples)

    if condition_col not in df.columns:
        raise ValueError(
            f"Condition column '{condition_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    if metric not in df.columns:
        raise ValueError(
            f"Metric '{metric}' not found. "
            f"Available metrics: {[c for c in df.columns if 'pct' in c or 'rate' in c]}"
        )

    # Get stats for each condition
    stats_dict = get_condition_stats(df, condition_col, metric)

    # Filter conditions with insufficient replicates
    stats_dict = {k: v for k, v in stats_dict.items() if v.n >= min_replicates}

    if not stats_dict:
        raise ValueError(
            f"No conditions have at least {min_replicates} replicates"
        )

    # Determine base condition
    if base_condition is None:
        base_condition = sorted(stats_dict.keys())[0]

    if base_condition not in stats_dict:
        raise ValueError(
            f"Base condition '{base_condition}' not found or has < {min_replicates} replicates"
        )

    # Perform comparisons
    return compare_conditions(
        stats_dict=stats_dict,
        base_condition=base_condition,
        metric=metric,
        fdr_correction=fdr_correction,
    )


def get_condition_summary(
    results: List[SampleResult],
    samples: Optional[List[Sample]] = None,
    condition_col: str = "condition",
    metrics: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Generate summary statistics for each condition.

    Args:
        results: List of SampleResult objects
        samples: Optional list of Sample objects for metadata
        condition_col: Column name containing condition labels
        metrics: List of metrics to summarize. Default includes all rate metrics.

    Returns:
        DataFrame with mean, std, sem, n for each metric by condition

    Example:
        >>> summary = get_condition_summary(results, samples)
        >>> print(summary[['condition', 'dedup_hdr_pct_mean', 'dedup_hdr_pct_sem']])
    """
    if metrics is None:
        metrics = [
            "dedup_hdr_pct",
            "dedup_nhej_pct",
            "dedup_wt_pct",
            "dedup_large_del_pct",
            "kmer_hdr_rate",
        ]

    df = results_to_dataframe(results, samples)

    if condition_col not in df.columns:
        raise ValueError(
            f"Condition column '{condition_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    # Calculate summary stats
    summary_rows = []
    for condition in df[condition_col].unique():
        subset = df[df[condition_col] == condition]
        row = {"condition": condition, "n": len(subset)}

        for metric in metrics:
            if metric in subset.columns:
                values = subset[metric].dropna()
                row[f"{metric}_mean"] = values.mean() if len(values) > 0 else None
                row[f"{metric}_std"] = values.std() if len(values) > 1 else None
                row[f"{metric}_sem"] = (
                    values.std() / (len(values) ** 0.5) if len(values) > 1 else None
                )

        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def compare_dataframe_by_condition(
    df: pd.DataFrame,
    condition_col: str = "condition",
    metric: str = "dedup_hdr_pct",
    base_condition: Optional[str] = None,
    fdr_correction: bool = True,
    min_replicates: int = 2,
) -> ComparisonSet:
    """
    Compare a metric across conditions from an existing DataFrame.

    Alternative to compare_metric_by_condition that takes a DataFrame directly.

    Args:
        df: DataFrame with sample data (must have condition_col and metric columns)
        condition_col: Column name containing condition labels
        metric: Metric to compare
        base_condition: Reference condition for comparison
        fdr_correction: Apply Benjamini-Hochberg FDR correction
        min_replicates: Minimum replicates required for statistical test

    Returns:
        ComparisonSet with statistical results
    """
    if condition_col not in df.columns:
        raise ValueError(
            f"Condition column '{condition_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    if metric not in df.columns:
        raise ValueError(f"Metric '{metric}' not found in DataFrame columns")

    # Get stats for each condition
    stats_dict = get_condition_stats(df, condition_col, metric)

    # Filter conditions with insufficient replicates
    stats_dict = {k: v for k, v in stats_dict.items() if v.n >= min_replicates}

    if not stats_dict:
        raise ValueError(
            f"No conditions have at least {min_replicates} replicates"
        )

    # Determine base condition
    if base_condition is None:
        base_condition = sorted(stats_dict.keys())[0]

    if base_condition not in stats_dict:
        raise ValueError(
            f"Base condition '{base_condition}' not found or has < {min_replicates} replicates"
        )

    # Perform comparisons
    return compare_conditions(
        stats_dict=stats_dict,
        base_condition=base_condition,
        metric=metric,
        fdr_correction=fdr_correction,
    )
