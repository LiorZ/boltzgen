"""Base classes for structural backbone filters."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class FilterResult:
    """Result from a single filter applied to a single structure.

    Attributes:
        passed: Whether the structure passed this filter.
        metrics: Numerical metrics computed by the filter (e.g. distances, counts).
        details: Human-readable summary of the result.
    """

    passed: bool
    metrics: dict[str, float] = field(default_factory=dict)
    details: str = ""


class BaseFilter(ABC):
    """Abstract base class for structural backbone filters.

    Subclasses implement `run()` which takes a CIF path and returns a
    FilterResult. All configurable thresholds are set via __init__.
    Filters are stateless — run() has no side effects.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Short identifier for this filter (used as CSV column prefix)."""

    @abstractmethod
    def run(self, cif_path: Path) -> FilterResult:
        """Apply this filter to a single structure.

        Args:
            cif_path: Path to a CIF file containing the complex.

        Returns:
            FilterResult with pass/fail, numerical metrics, and details.
        """
