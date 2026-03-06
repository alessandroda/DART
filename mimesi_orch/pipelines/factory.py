"""Pipeline selection and instantiation."""

from __future__ import annotations

from config_models import AppConfig
from mimesi_types import ModelType
from pipelines.registry import create_pipeline

# Import modules for registration side effects.
from pipelines.chimere2017 import pipeline as _chimere2017_pipeline  # noqa: F401


def infer_pipeline_name(config: AppConfig) -> str:
    if config.pipeline is not None and config.pipeline.name:
        requested = config.pipeline.name.strip().lower()
        if requested != "chimere2017":
            raise ValueError(
                f"Unsupported pipeline '{requested}'. Only 'chimere2017' is enabled."
            )
        return requested

    if config.assimilation.model_type != ModelType.CHIMERE:
        raise ValueError(
            f"Unsupported model_type={config.assimilation.model_type}. "
            "Only CHIMERE 2017 is enabled."
        )

    if config.model_data is not None and config.model_data.control_run_exp_name:
        raise ValueError(
            "CHIMERE 2023-style configuration detected "
            "(model_data.control_run_exp_name is set), "
            "but only CHIMERE 2017 is enabled."
        )
    return "chimere2017"


def build_pipeline(config: AppConfig, time_manager: object):
    name = infer_pipeline_name(config)
    return create_pipeline(name, config, time_manager)
