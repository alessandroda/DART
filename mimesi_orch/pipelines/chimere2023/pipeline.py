from __future__ import annotations

from pipelines.chimere_irene_pipeline import ChimereV2023DartPipeline
from pipelines.registry import register_pipeline
from pipelines.chimere2023.config import Chimere2023PipelineConfig
from pipelines.chimere2023.paths import Chimere2023Paths


class Chimere2023Pipeline(ChimereV2023DartPipeline):
    """Concrete CHIMERE 2023 pipeline using existing execution logic."""



def _build(config, time_manager):
    p_cfg = Chimere2023PipelineConfig.from_app_config(config)
    p_paths = Chimere2023Paths(p_cfg)
    return Chimere2023Pipeline(time_manager, p_paths.manager, config)


register_pipeline("chimere2023", _build)
