"""
Pydantic configuration models for MIMESI orchestration.

This module defines the validated, typed configuration schema used by all
pipelines (FARM, CHIMERE, future models).

Design principles:
- Strict separation between configuration and logic
- Early validation of paths, dates, and numeric constraints
- One single root config object (AppConfig)
"""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator
from pathlib import Path
from datetime import datetime
from typing import Any, Literal, Optional
from mimesi_types import ModelType, Scheduler

# ---------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------


class PathsConfig(BaseModel):
    """
    Filesystem paths used by the orchestration.
    All paths are relative to base_path unless explicitly absolute.
    """

    base_path: Path
    env_python: Path

    listing_file: Path

    run_submit_model_template: Optional[Path] = None
    run_submit_farm_template: Optional[Path] = None
    path_submit_bsh: Path
    path_filter: Path
    path_data: Path
    path_control_run: Optional[Path] = None
    chimere_par_template: Optional[Path] = None
    run_submit_replace_perturbations: Optional[Path] = None
    path_perturbed_emi: Optional[Path] = None
    log_directory: Optional[Path] = None

    @field_validator("*", mode="before")
    @classmethod
    def expand_paths(cls, v):
        if isinstance(v, str):
            return Path(v)
        return v


# ---------------------------------------------------------------------
# TIME
# ---------------------------------------------------------------------


class TimeConfig(BaseModel):
    """
    Time control of the simulation–assimilation loop.
    """

    start_time: datetime
    end_time: datetime
    dt_seconds: int = Field(gt=0)

    backup_perturb_days: Optional[int] = Field(default=None, ge=0)
    backup_ic_hours: Optional[int] = Field(default=None, ge=0)
    backup_ic_option: Optional[str] = None

    @field_validator("end_time")
    @classmethod
    def check_time_order(cls, end_time, info):
        start_time = info.data.get("start_time")
        if start_time and end_time <= start_time:
            raise ValueError("end_time must be after start_time")
        return end_time


# ---------------------------------------------------------------------
# ASSIMILATION
# ---------------------------------------------------------------------


class AssimilationConfig(BaseModel):
    """
    Data assimilation settings.
    """

    model_type: ModelType
    ass_var: str
    state_variable_qty: str
    obs_type: str

    case_dir: str
    case_emi_dir: Optional[str] = None
    emi_perturbation_dir : Optional[str] = None
    emi_var: Optional[str] = None

    no_mems: int = Field(gt=0)
    run_assimilation_flag: bool = True

    @field_validator("model_type", mode="before")
    @classmethod
    def normalize_model_type(cls, v):
        if isinstance(v, str):
            return v.lower()
        return v

# ---------------------------------------------------------------------
# CLUSTER
# ---------------------------------------------------------------------


class ClusterConfig(BaseModel):
    cluster_queue: str
    scheduler: Scheduler
    project_name: str
    nb_proc: Optional[int] = Field(default=None, gt=0)
    walltime: Optional[int] = Field(default=None, gt=0)
    mail: Optional[str] = None

    @field_validator("scheduler", mode="before")
    @classmethod
    def normalize_scheduler(cls, v):
        if isinstance(v, str):
            return v.lower()
        return v
# ---------------------------------------------------------------------
# DART
# ---------------------------------------------------------------------


class DartConfig(BaseModel):
    """
    DART-specific paths and directories.
    """

    work_dir: Path
    analysis_dir: Path
    preassim_dir: Path
    posteriors_dir: Path
    obs_converters_dir: Path

    @field_validator("*", mode="before")
    @classmethod
    def expand_paths(cls, v):
        if isinstance(v, str):
            return Path(v)
        return v


# ---------------------------------------------------------------------
# LOGGING
# ---------------------------------------------------------------------


class LoggingConfig(BaseModel):
    """
    Logging configuration.
    """

    level: str = "INFO"
    format: str = "%(asctime)s %(message)s"


# ---------------------------------------------------------------------
# SATELLITE DATA
# ---------------------------------------------------------------------


class SatelliteDataConfig(BaseModel):
    """
    Satellite observation handling.
    """

    search_window_seconds: int = Field(gt=0)

# ---------------------------------------------------------------------
# MODEL DATA
# ---------------------------------------------------------------------


class ModelDataConfig(BaseModel):
    """
    Satellite observation handling.
    """

    control_run_exp_name: Optional[str] = None
    domain: Optional[str] = None
    chimpart: Optional[str] = None
    submit_sequentially: Optional[bool] = None
    restart_from_controlrun: Optional[bool] = None
    perturbation_names: Optional[list] = None
    ensemble_list: Optional[list] = None


class PipelineSelectionConfig(BaseModel):
    """
    Explicit pipeline selection override.
    """

    name: str
    config: Optional[dict[str, Any]] = None


class MonitoringConfig(BaseModel):
    """
    Runtime monitoring and scheduler interaction settings.

    These parameters control how the orchestrator:
    - monitors submitted jobs
    - waits for completion
    - allocates computational resources
    """

    cores: int = Field(gt=0)
    job_check_interval: int = Field(gt=0, description="Polling interval in seconds")


class RuntimeConfig(BaseModel):
    """
    Runtime profile and process-level behaviors.
    """

    profile: str = "mimesi"
    log_dir_mode: Literal["path_data", "cwd"] = "path_data"
    log_dir_name: Optional[str] = None

    @field_validator("profile", mode="before")
    @classmethod
    def normalize_profile(cls, v):
        if isinstance(v, str):
            return v.strip().lower()
        return v

    @field_validator("log_dir_name", mode="before")
    @classmethod
    def normalize_log_dir_name(cls, v):
        if isinstance(v, str):
            v = v.strip()
            return v if v else None
        return v


# ---------------------------------------------------------------------
# ROOT CONFIG
# ---------------------------------------------------------------------


class AppConfig(BaseModel):
    """
    Root application configuration.

    This is the single entry point for all orchestration logic.
    """

    paths: PathsConfig
    time: TimeConfig
    assimilation: AssimilationConfig
    cluster: ClusterConfig

    dart: Optional[DartConfig] = None
    logging: Optional[LoggingConfig] = None
    satellite_data: Optional[SatelliteDataConfig] = None
    model_data: Optional[ModelDataConfig] = None
    monitoring: Optional[MonitoringConfig] = None
    runtime: Optional[RuntimeConfig] = None
    pipeline: Optional[PipelineSelectionConfig] = None

    # runtime-only (not from YAML)
    _config_path: Optional[Path] = None

    class Config:
        extra = "forbid"

AppConfig.model_rebuild()
