"""Pipeline selection and instantiation."""

from __future__ import annotations

from config_models import AppConfig
from mimesi_types import ModelType
from pipelines.registry import create_pipeline

# Import modules for registration side effects.
from pipelines.chimere2017 import pipeline as _chimere2017_pipeline  # noqa: F401
from pipelines.chimere2023 import pipeline as _chimere2023_pipeline  # noqa: F401


def infer_pipeline_name(config: AppConfig) -> str:
    if config.pipeline is not None and config.pipeline.name:
        requested = config.pipeline.name.strip().lower()
        if requested not in ["chimere2017", "chimere2023"]:
            raise ValueError(
                f"Unsupported pipeline '{requested}'. Only 'chimere2017' and 'chimere2023' are enabled."
            )
        return requested

    if config.assimilation.model_type != ModelType.CHIMERE:
        raise ValueError(
            f"Unsupported model_type={config.assimilation.model_type}. "
            "Only CHIMERE models are enabled."
        )

    # Infer chimere2023 if the specific 2023 config fields are present
    if config.model_data is not None and config.model_data.control_run_exp_name:
        return "chimere2023"

    # Default fallback
    return "chimere2017"


def build_pipeline(config: AppConfig, time_manager: object):
    name = infer_pipeline_name(config)
    return create_pipeline(name, config, time_manager)
