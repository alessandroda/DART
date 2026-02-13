from __future__ import annotations

from paths import PathManager
from pipelines.farm.config import FarmPipelineConfig


class FarmPaths:
    """Pipeline-specific path facade for FARM."""

    def __init__(self, cfg: FarmPipelineConfig):
        self._manager = PathManager.model_from_pipeline(cfg)

    @property
    def manager(self) -> PathManager:
        return self._manager
