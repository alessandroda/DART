"""Pipeline abstraction contracts.

Shared infrastructure should depend on these contracts, not on concrete pipelines.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Protocol


class PipelinePaths(Protocol):
    """Marker protocol for pipeline-specific path access APIs."""


class PipelineConfig(Protocol):
    """Marker protocol for pipeline-specific config models."""


class BasePipeline(ABC):
    """Minimal interface implemented by all concrete pipelines."""

    @abstractmethod
    def run_pipeline(self) -> None:
        """Run the full pipeline workflow."""


class PipelineBuilder(ABC):
    """Factory helper for creating pipeline-specific components."""

    name: str

    @classmethod
    @abstractmethod
    def build(cls, config: Any, time_manager: Any) -> BasePipeline:
        """Create a configured pipeline instance."""
