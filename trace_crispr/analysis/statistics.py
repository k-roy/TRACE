"""
Statistical analysis functions for TRACE.

Provides t-tests, ANOVA, and multiple testing correction (FDR).

Author: Kevin R. Roy
"""

from typing import List, Dict, Optional, Tuple
import numpy as np
import pandas as pd

from .types import ConditionStats, ComparisonResult, ComparisonSet


def ttest_vs_base(
    condition_stats: ConditionStats,
    base_stats: ConditionStats,
    test_type: str = "independent",
) -> Tuple[float, float]:
    """
    Perform t-test comparing condition to base.

    Args:
        condition_stats: Stats for the test condition
        base_stats: Stats for the reference/control condition
        test_type: 'independent' for two-sample t-test (Welch's),
                   'onesample' for one-sample t-test vs base mean

    Returns:
        Tuple of (t_statistic, p_value)

    Raises:
        ImportError: If scipy is not installed
    """
    try:
        from scipy import stats
    except ImportError:
        raise ImportError(
            "scipy is required for statistical tests. "
            "Install with: pip install trace-crispr[visualization]"
        )

    if test_type == "independent":
        if len(condition_stats.values) < 2 or len(base_stats.values) < 2:
            return np.nan, np.nan
        t_stat, p_value = stats.ttest_ind(
            condition_stats.values, base_stats.values, equal_var=False  # Welch's t-test
        )
    elif test_type == "onesample":
        if len(condition_stats.values) < 2:
            return np.nan, np.nan
        t_stat, p_value = stats.ttest_1samp(condition_stats.values, base_stats.mean)
    else:
        raise ValueError(f"Unknown test_type: {test_type}")

    return float(t_stat), float(p_value)


def p_value_to_stars(
    p_value: float, thresholds: Optional[Dict[str, float]] = None
) -> str:
    """
    Convert p-value to significance stars.

    Args:
        p_value: The p-value to convert
        thresholds: Optional custom thresholds dict
                   Default: {'***': 0.001, '**': 0.01, '*': 0.05}

    Returns:
        Significance string ('***', '**', '*', or 'ns')

    Example:
        >>> p_value_to_stars(0.003)
        '**'
        >>> p_value_to_stars(0.08)
        'ns'
    """
    if thresholds is None:
        thresholds = {"***": 0.001, "**": 0.01, "*": 0.05}

    if pd.isna(p_value):
        return ""

    for stars, threshold in sorted(thresholds.items(), key=lambda x: x[1]):
        if p_value < threshold:
            return stars

    return "ns"


def benjamini_hochberg(p_values: List[float], alpha: float = 0.05) -> List[float]:
    """
    Apply Benjamini-Hochberg FDR correction.

    Args:
        p_values: List of p-values
        alpha: FDR threshold (default: 0.05)

    Returns:
        List of adjusted p-values

    Example:
        >>> p_adj = benjamini_hochberg([0.01, 0.04, 0.03, 0.5])
        >>> [f'{p:.3f}' for p in p_adj]
        ['0.040', '0.053', '0.053', '0.500']
    """
    p_array = np.array(p_values)
    n = len(p_array)

    # Handle NaN values
    valid_mask = ~np.isnan(p_array)
    valid_p = p_array[valid_mask]
    n_valid = len(valid_p)

    if n_valid == 0:
        return list(p_values)

    # Sort p-values and get ranking
    sorted_indices = np.argsort(valid_p)
    sorted_p = valid_p[sorted_indices]

    # Calculate BH adjusted p-values
    # p_adj[i] = min(p[i] * n / rank[i], p_adj[i+1])
    adjusted = np.zeros(n_valid)
    adjusted[-1] = sorted_p[-1]

    for i in range(n_valid - 2, -1, -1):
        rank = i + 1
        adjusted[i] = min(sorted_p[i] * n_valid / rank, adjusted[i + 1])

    # Cap at 1.0
    adjusted = np.minimum(adjusted, 1.0)

    # Unsort
    unsorted_adjusted = np.zeros(n_valid)
    unsorted_adjusted[sorted_indices] = adjusted

    # Put back into original array with NaN positions
    result = np.full(n, np.nan)
    result[valid_mask] = unsorted_adjusted

    return result.tolist()


def compare_conditions(
    stats_dict: Dict[str, ConditionStats],
    base_condition: str,
    metric: str,
    test_type: str = "independent",
    fdr_correction: bool = True,
    alpha: float = 0.05,
) -> ComparisonSet:
    """
    Compare all conditions to a base condition with statistical tests.

    Args:
        stats_dict: Dict mapping condition names to ConditionStats
        base_condition: Name of the reference/control condition
        metric: Name of the metric being compared
        test_type: Type of t-test ('independent' or 'onesample')
        fdr_correction: Whether to apply Benjamini-Hochberg FDR correction
        alpha: Significance threshold

    Returns:
        ComparisonSet containing all comparison results

    Example:
        >>> stats = get_condition_stats(df, 'condition', 'dedup_hdr_pct')
        >>> comparisons = compare_conditions(stats, base_condition='control')
        >>> print(f"Significant: {comparisons.significant_conditions}")
    """
    if base_condition not in stats_dict:
        raise ValueError(f"Base condition '{base_condition}' not found in stats_dict")

    base_stats = stats_dict[base_condition]
    comparisons = []

    # Perform comparisons for all non-base conditions
    for condition, condition_stats in stats_dict.items():
        if condition == base_condition:
            continue

        t_stat, p_value = ttest_vs_base(condition_stats, base_stats, test_type)

        # Calculate fold change
        if base_stats.mean > 0:
            fold_change = condition_stats.mean / base_stats.mean
        else:
            fold_change = None

        comparison = ComparisonResult(
            condition=condition,
            base_condition=base_condition,
            metric=metric,
            condition_stats=condition_stats,
            base_stats=base_stats,
            difference=condition_stats.mean - base_stats.mean,
            fold_change=fold_change,
            t_statistic=t_stat,
            p_value=p_value,
        )
        comparisons.append(comparison)

    # Apply FDR correction
    if fdr_correction and comparisons:
        p_values = [c.p_value for c in comparisons]
        p_adjusted = benjamini_hochberg(p_values, alpha)

        for c, p_adj in zip(comparisons, p_adjusted):
            c.p_adjusted = p_adj
            c.significance = p_value_to_stars(p_adj)
    else:
        for c in comparisons:
            c.significance = p_value_to_stars(c.p_value)

    return ComparisonSet(
        comparisons=comparisons,
        metric=metric,
        base_condition=base_condition,
        fdr_method="benjamini-hochberg" if fdr_correction else None,
        alpha=alpha,
    )
