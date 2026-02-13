from __future__ import annotations

from paths import PathManager
from pipelines.chimere2017.config import Chimere2017PipelineConfig


class Chimere2017Paths:
    """Pipeline-specific path facade for CHIMERE 2017."""

    def __init__(self, cfg: Chimere2017PipelineConfig):
        self._manager = PathManager.model_from_pipeline(cfg)

    @property
    def manager(self) -> PathManager:
        return self._manager
