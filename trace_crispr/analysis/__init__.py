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

from .types import (
    ConditionStats,
    ComparisonResult,
    ComparisonSet,
)

from .aggregation import (
    results_to_dataframe,
    get_condition_stats,
    aggregate_by_group,
)

from .statistics import (
    ttest_vs_base,
    p_value_to_stars,
    benjamini_hochberg,
    compare_conditions,
)

from .comparison import (
    compare_metric_by_condition,
    compare_dataframe_by_condition,
    get_condition_summary,
)

# Plotting functions - import only if visualization deps available
try:
    from .plotting import (
        plot_condition_comparison,
        plot_comparison_summary,
        plot_replicate_correlation,
        plot_multi_metric_comparison,
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
