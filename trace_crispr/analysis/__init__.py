"""
Analysis and visualization module for TRACE.

Provides functions for:
- Aggregating sample results by condition
- Statistical comparison between conditions
- Publication-quality visualizations

Author: Kevin R. Roy

Example usage:

    from trace_crispr.analysis import (
        results_to_dataframe,
        compare_metric_by_condition,
        plot_condition_comparison,
        get_condition_stats,
    )

    # Run TRACE pipeline
    results = pipeline.run_all(samples)

    # Compare HDR rates by condition
    comparisons = compare_metric_by_condition(
        results, samples,
        condition_col='treatment',
        metric='dedup_hdr_pct',
        base_condition='untreated'
    )

    # View results
    print(comparisons.to_dataframe())

    # Create visualization
    df = results_to_dataframe(results, samples)
    stats = get_condition_stats(df, 'treatment', 'dedup_hdr_pct')
    fig = plot_condition_comparison(
        stats, comparisons,
        title='HDR Rate by Treatment'
    )
    fig.savefig('comparison.png')
"""

from .aggregation import (
    aggregate_by_group,
    get_condition_stats,
    results_to_dataframe,
)
from .comparison import (
    compare_dataframe_by_condition,
    compare_metric_by_condition,
    get_condition_summary,
)
from .statistics import (
    benjamini_hochberg,
    compare_conditions,
    p_value_to_stars,
    ttest_vs_base,
)
from .types import (
    ComparisonResult,
    ComparisonSet,
    ConditionStats,
)

# Plotting functions - import only if visualization deps available
try:
    from .plotting import (
        plot_comparison_summary,
        plot_condition_comparison,
        plot_multi_metric_comparison,
        plot_replicate_correlation,
    )

    _PLOTTING_AVAILABLE = True
except ImportError:
    _PLOTTING_AVAILABLE = False

    def plot_condition_comparison(*args, **kwargs):
        raise ImportError(
            "Plotting requires matplotlib and seaborn. "
            "Install with: pip install trace-crispr[visualization]"
        )

    def plot_comparison_summary(*args, **kwargs):
        raise ImportError(
            "Plotting requires matplotlib and seaborn. "
            "Install with: pip install trace-crispr[visualization]"
        )

    def plot_replicate_correlation(*args, **kwargs):
        raise ImportError(
            "Plotting requires matplotlib and seaborn. "
            "Install with: pip install trace-crispr[visualization]"
        )

    def plot_multi_metric_comparison(*args, **kwargs):
        raise ImportError(
            "Plotting requires matplotlib and seaborn. "
            "Install with: pip install trace-crispr[visualization]"
        )


__all__ = [
    # Types
    "ConditionStats",
    "ComparisonResult",
    "ComparisonSet",
    # Aggregation
    "results_to_dataframe",
    "get_condition_stats",
    "aggregate_by_group",
    # Statistics
    "ttest_vs_base",
    "p_value_to_stars",
    "benjamini_hochberg",
    "compare_conditions",
    # Comparison workflows
    "compare_metric_by_condition",
    "compare_dataframe_by_condition",
    "get_condition_summary",
    # Plotting
    "plot_condition_comparison",
    "plot_comparison_summary",
    "plot_replicate_correlation",
    "plot_multi_metric_comparison",
]
