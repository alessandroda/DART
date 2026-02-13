from __future__ import annotations

from paths import PathManager
from pipelines.chimere2023.config import Chimere2023PipelineConfig


class Chimere2023Paths:
    """Pipeline-specific path facade for CHIMERE 2023."""

    def __init__(self, cfg: Chimere2023PipelineConfig):
        self._manager = PathManager.model_from_pipeline(cfg)

    @property
    def manager(self) -> PathManager:
        return self._manager
