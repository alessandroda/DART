from __future__ import annotations

from pipelines.chimere_pipeline import Chimere2017DartPipeline
from pipelines.registry import register_pipeline
from pipelines.chimere2017.config import Chimere2017PipelineConfig
from pipelines.chimere2017.paths import Chimere2017Paths


class Chimere2017Pipeline(Chimere2017DartPipeline):
    """Concrete CHIMERE 2017 pipeline using existing execution logic."""



def _build(config, time_manager):
    p_cfg = Chimere2017PipelineConfig.from_app_config(config)
    p_paths = Chimere2017Paths(p_cfg)
    return Chimere2017Pipeline(time_manager, p_paths.manager, config)


register_pipeline("chimere2017", _build)
