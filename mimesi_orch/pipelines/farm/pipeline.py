from __future__ import annotations

from pipelines.farm_pipeline import FarmDartPipeline
from pipelines.registry import register_pipeline
from pipelines.farm.config import FarmPipelineConfig
from pipelines.farm.paths import FarmPaths


class FarmPipeline(FarmDartPipeline):
    """Concrete FARM pipeline implemented on top of existing execution logic."""



def _build(config, time_manager):
    p_cfg = FarmPipelineConfig.from_app_config(config)
    p_paths = FarmPaths(p_cfg)
    return FarmPipeline(time_manager, p_paths.manager, config)


register_pipeline("farm", _build)
