"""
Type definitions for TRACE analysis module.

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
import pandas as pd
import numpy as np


@dataclass
class ConditionStats:
    """
    Statistics for a single condition group.

    Attributes:
        condition: Name/identifier of the condition
        metric: Name of the metric being summarized (e.g., 'dedup_hdr_pct')
        mean: Mean value across replicates
        std: Standard deviation
        sem: Standard error of the mean
        n: Number of replicates
        values: List of individual values for scatter plotting
        sample_ids: List of sample IDs corresponding to values
        metadata: Additional metadata (e.g., dates, batches)
    """

    condition: str
    metric: str
    mean: float
    std: float
    sem: float
    n: int
    values: List[float] = field(default_factory=list)
    sample_ids: List[str] = field(default_factory=list)
    metadata: Dict[str, List[Any]] = field(default_factory=dict)

    @classmethod
    def from_series(
        cls,
        condition: str,
        metric: str,
        data: pd.Series,
        sample_ids: Optional[List[str]] = None,
        metadata: Optional[Dict[str, List[Any]]] = None,
    ) -> "ConditionStats":
        """Create ConditionStats from a pandas Series."""
        values = data.dropna().tolist()
        n = len(values)
        mean = float(np.mean(values)) if n > 0 else 0.0
        std = float(np.std(values, ddof=1)) if n > 1 else 0.0
        sem = std / np.sqrt(n) if n > 1 else 0.0

        return cls(
            condition=condition,
            metric=metric,
            mean=mean,
            std=std,
            sem=sem,
            n=n,
            values=values,
            sample_ids=sample_ids or [],
            metadata=metadata or {},
        )


@dataclass
class ComparisonResult:
    """
    Result of statistical comparison between conditions.

    Attributes:
        condition: The test condition being compared
        base_condition: The reference/control condition
        metric: The metric being compared
        condition_stats: Stats for the test condition
        base_stats: Stats for the reference condition
        difference: Difference in means (condition - base)
        fold_change: Ratio of means (condition / base)
        t_statistic: T-test statistic
        p_value: Uncorrected p-value
        p_adjusted: FDR-corrected p-value (if applicable)
        significance: Significance level ('***', '**', '*', 'ns')
        test_type: Type of test performed
    """

    condition: str
    base_condition: str
    metric: str
    condition_stats: ConditionStats
    base_stats: ConditionStats
    difference: float
    fold_change: Optional[float]
    t_statistic: float
    p_value: float
    p_adjusted: Optional[float] = None
    significance: str = "ns"
    test_type: str = "ttest_ind"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        return {
            "condition": self.condition,
            "base_condition": self.base_condition,
            "metric": self.metric,
            "condition_mean": self.condition_stats.mean,
            "condition_std": self.condition_stats.std,
            "condition_n": self.condition_stats.n,
            "base_mean": self.base_stats.mean,
            "base_std": self.base_stats.std,
            "base_n": self.base_stats.n,
            "difference": self.difference,
            "fold_change": self.fold_change,
            "t_statistic": self.t_statistic,
            "p_value": self.p_value,
            "p_adjusted": self.p_adjusted,
            "significance": self.significance,
        }


@dataclass
class ComparisonSet:
    """
    Collection of comparison results with summary statistics.

    Attributes:
        comparisons: List of individual ComparisonResults
        metric: The metric being compared
        base_condition: The reference condition used
        fdr_method: FDR correction method applied (if any)
        alpha: Significance threshold used
    """

    comparisons: List[ComparisonResult]
    metric: str
    base_condition: str
    fdr_method: Optional[str] = None
    alpha: float = 0.05

    def to_dataframe(self) -> pd.DataFrame:
        """Convert all comparisons to a DataFrame."""
        rows = [c.to_dict() for c in self.comparisons]
        return pd.DataFrame(rows)

    @property
    def significant_conditions(self) -> List[str]:
        """Return conditions with significant differences."""
        return [c.condition for c in self.comparisons if c.significance != "ns"]

    @property
    def n_significant(self) -> int:
        """Count of significant comparisons."""
        return len(self.significant_conditions)
