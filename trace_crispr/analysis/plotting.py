"""
Plotting functions for TRACE analysis.

Provides publication-quality visualizations for condition comparisons.
Requires optional dependencies: matplotlib, seaborn

Author: Kevin R. Roy
"""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .types import ComparisonSet, ConditionStats


def _check_plotting_deps():
    """Check that plotting dependencies are available."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        return plt, sns
    except ImportError:
        raise ImportError(
            "Plotting requires matplotlib and seaborn. "
            "Install with: pip install trace-crispr[visualization]"
        )


def plot_condition_comparison(
    stats_dict: Dict[str, ConditionStats],
    comparisons: Optional[ComparisonSet] = None,
    base_condition: Optional[str] = None,
    title: Optional[str] = None,
    ylabel: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
    bar_color: str = "#22c55e",
    base_color: str = "#888888",
    point_color: str = "black",
    point_alpha: float = 0.6,
    jitter: float = 0.15,
    show_significance: bool = True,
    condition_order: Optional[List[str]] = None,
    ax: Optional[Any] = None,
) -> Any:
    """
    Create bar plot with individual data points overlay.

    Generates a bar chart showing mean values with error bars (SEM),
    overlaid with individual replicate points with horizontal jitter.
    If comparisons provided, shows significance stars above bars.

    Args:
        stats_dict: Dict mapping condition names to ConditionStats
        comparisons: Optional ComparisonSet with statistical results
        base_condition: Highlight the base condition differently
        title: Plot title
        ylabel: Y-axis label (default: metric name from first ConditionStats)
        figsize: Figure size as (width, height)
        bar_color: Color for treatment bars (hex or named color)
        base_color: Color for base condition bar
        point_color: Color for individual data points
        point_alpha: Transparency for individual points
        jitter: Amount of horizontal jitter for points (0-0.5)
        show_significance: Show significance stars if comparisons provided
        condition_order: Optional list specifying condition order
        ax: Optional matplotlib Axes to plot on

    Returns:
        matplotlib Figure object

    Example:
        >>> stats = get_condition_stats(df, 'condition', 'dedup_hdr_pct')
        >>> comparisons = compare_conditions(stats, base_condition='control')
        >>> fig = plot_condition_comparison(
        ...     stats, comparisons,
        ...     title='HDR Rate by Treatment',
        ...     ylabel='HDR Rate (%)'
        ... )
        >>> fig.savefig('hdr_comparison.png', dpi=150)
    """
    plt, sns = _check_plotting_deps()

    # Determine condition order
    if condition_order is None:
        condition_order = sorted(stats_dict.keys())

    # Filter to conditions in stats_dict
    condition_order = [c for c in condition_order if c in stats_dict]

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    x_positions = np.arange(len(condition_order))
    means = [stats_dict[c].mean for c in condition_order]
    sems = [stats_dict[c].sem for c in condition_order]

    # Determine bar colors
    bar_colors = []
    for cond in condition_order:
        if cond == base_condition:
            bar_colors.append(base_color)
        else:
            bar_colors.append(bar_color)

    # Plot bars
    ax.bar(
        x_positions,
        means,
        yerr=sems,
        capsize=3,
        color=bar_colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=1,
    )

    # Overlay individual points with jitter
    for i, cond in enumerate(condition_order):
        stats = stats_dict[cond]
        if stats.values:
            x_jitter = np.random.uniform(-jitter, jitter, len(stats.values))
            ax.scatter(
                i + x_jitter,
                stats.values,
                color=point_color,
                alpha=point_alpha,
                s=30,
                zorder=5,
                edgecolors="white",
                linewidth=0.5,
            )

    # Add significance stars
    if show_significance and comparisons is not None:
        max_y = max(
            stats_dict[c].mean + stats_dict[c].sem for c in condition_order
        )
        y_offset = max_y * 0.05

        for comparison in comparisons.comparisons:
            if comparison.significance and comparison.significance != "ns":
                try:
                    idx = condition_order.index(comparison.condition)
                    y_pos = (
                        stats_dict[comparison.condition].mean
                        + stats_dict[comparison.condition].sem
                        + y_offset
                    )
                    ax.text(
                        idx,
                        y_pos,
                        comparison.significance,
                        ha="center",
                        va="bottom",
                        fontsize=12,
                        fontweight="bold",
                    )
                except ValueError:
                    pass  # Condition not in order

    # Labels and formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(condition_order, rotation=45, ha="right")

    if ylabel is None and stats_dict:
        ylabel = list(stats_dict.values())[0].metric
    ax.set_ylabel(ylabel, fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)

    ax.grid(True, alpha=0.3, axis="y")
    ax.set_axisbelow(True)

    plt.tight_layout()

    return fig


def plot_comparison_summary(
    comparison_set: ComparisonSet,
    plot_type: str = "bar",
    show_error_bars: bool = True,
    figsize: Tuple[float, float] = (12, 6),
    ax: Optional[Any] = None,
) -> Any:
    """
    Create summary plot from ComparisonSet results.

    Args:
        comparison_set: ComparisonSet from compare_conditions()
        plot_type: 'bar' for bar chart, 'forest' for forest plot
        show_error_bars: Include error bars
        figsize: Figure size
        ax: Optional matplotlib Axes

    Returns:
        matplotlib Figure object
    """
    plt, sns = _check_plotting_deps()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    df = comparison_set.to_dataframe()
    df = df.sort_values("fold_change", ascending=True)

    if plot_type == "bar":
        # Horizontal bar plot showing fold change
        colors = ["#22c55e" if fc > 1 else "#ef4444" for fc in df["fold_change"]]

        y_pos = np.arange(len(df))
        ax.barh(y_pos, df["fold_change"], color=colors, alpha=0.7)

        # Add significance markers
        for i, (_, row) in enumerate(df.iterrows()):
            if row["significance"] and row["significance"] != "ns":
                ax.text(
                    row["fold_change"] + 0.05,
                    i,
                    row["significance"],
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                )

        ax.axvline(1.0, color="black", linestyle="--", alpha=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df["condition"])
        ax.set_xlabel(f"Fold Change vs {comparison_set.base_condition}")
        ax.set_title(f"{comparison_set.metric}: Fold Change Comparison")

    elif plot_type == "forest":
        # Forest plot with confidence intervals
        y_pos = np.arange(len(df))

        for i, (_, row) in enumerate(df.iterrows()):
            # Calculate 95% CI (approximation)
            ci_lower = row["condition_mean"] - 1.96 * row["condition_std"] / np.sqrt(
                row["condition_n"]
            )
            ci_upper = row["condition_mean"] + 1.96 * row["condition_std"] / np.sqrt(
                row["condition_n"]
            )

            color = "#22c55e" if row["significance"] != "ns" else "#888888"

            ax.errorbar(
                row["condition_mean"],
                i,
                xerr=[
                    [row["condition_mean"] - ci_lower],
                    [ci_upper - row["condition_mean"]],
                ],
                fmt="o",
                color=color,
                capsize=5,
                markersize=8,
            )

        # Add reference line for base mean
        base_mean = df["base_mean"].iloc[0]
        ax.axvline(
            base_mean,
            color="red",
            linestyle="--",
            alpha=0.5,
            label=f"{comparison_set.base_condition} mean",
        )

        ax.set_yticks(y_pos)
        ax.set_yticklabels(df["condition"])
        ax.set_xlabel(comparison_set.metric)
        ax.legend()

    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    return fig


def plot_replicate_correlation(
    df: pd.DataFrame,
    x_metric: str,
    y_metric: str,
    hue_col: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 8),
    ax: Optional[Any] = None,
) -> Any:
    """
    Create scatter plot comparing two metrics with correlation.

    Args:
        df: DataFrame with sample data
        x_metric: Column name for x-axis
        y_metric: Column name for y-axis
        hue_col: Optional column for color grouping
        figsize: Figure size
        ax: Optional matplotlib Axes

    Returns:
        matplotlib Figure object
    """
    plt, sns = _check_plotting_deps()

    try:
        from scipy import stats as scipy_stats
    except ImportError:
        raise ImportError(
            "scipy is required for correlation calculation. "
            "Install with: pip install trace-crispr[visualization]"
        )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Calculate correlation
    valid = df[[x_metric, y_metric]].dropna()
    if len(valid) > 2:
        r, p = scipy_stats.pearsonr(valid[x_metric], valid[y_metric])
    else:
        r, p = np.nan, np.nan

    # Create scatter plot
    if hue_col and hue_col in df.columns:
        sns.scatterplot(data=df, x=x_metric, y=y_metric, hue=hue_col, ax=ax, alpha=0.7)
    else:
        ax.scatter(df[x_metric], df[y_metric], alpha=0.7)

    # Add diagonal line
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.plot(lims, lims, "k--", alpha=0.3, zorder=0)

    # Add correlation annotation
    ax.text(
        0.05,
        0.95,
        f"r = {r:.3f}\np = {p:.2e}",
        transform=ax.transAxes,
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    ax.set_xlabel(x_metric)
    ax.set_ylabel(y_metric)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig


def plot_multi_metric_comparison(
    df: pd.DataFrame,
    condition_col: str,
    metrics: List[str],
    base_condition: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
) -> Any:
    """
    Create a multi-panel figure comparing multiple metrics across conditions.

    Args:
        df: DataFrame with sample data
        condition_col: Column name containing condition labels
        metrics: List of metrics to plot
        base_condition: Reference condition to highlight
        figsize: Figure size (default: auto-sized based on metric count)

    Returns:
        matplotlib Figure object
    """
    plt, sns = _check_plotting_deps()
    from .aggregation import get_condition_stats
    from .statistics import compare_conditions

    n_metrics = len(metrics)
    if figsize is None:
        figsize = (5 * n_metrics, 5)

    fig, axes = plt.subplots(1, n_metrics, figsize=figsize)

    # Handle single metric case
    if n_metrics == 1:
        axes = [axes]

    for i, metric in enumerate(metrics):
        stats_dict = get_condition_stats(df, condition_col, metric)

        # Get comparisons if base_condition specified
        comparisons = None
        if base_condition and base_condition in stats_dict:
            comparisons = compare_conditions(
                stats_dict, base_condition=base_condition, metric=metric
            )

        plot_condition_comparison(
            stats_dict,
            comparisons=comparisons,
            base_condition=base_condition,
            title=metric,
            ylabel=metric,
            ax=axes[i],
        )

    plt.tight_layout()
    return fig
