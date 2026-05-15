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
from typing import Optional
from mimesi_types import ModelType, Scheduler

# ---------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------


class PathsConfig(BaseModel):
    """
    Filesystem paths used by the orchestration.
    All paths are relative to base_path unless explicitly absolute.
    """

    base_path_ctm: Path
    base_path_DART: Path
    env_python: Path

    listing_file: Path

    run_submit_model_template:  Optional[Path] = None
    path_submit_bsh: Path
    path_filter: Path 
    path_data: Path
    case_name: Optional[Path] = None
    path_control_run: Optional[Path] = None

    run_submit_replace_perturbations: Optional[Path] = None
    path_perturbed_emi: Optional[Path] = None
    path_perturbed_meteo: Optional[Path] = None
    log_directory: Optional[Path] = None

    config_path: Optional[Path] = None

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
    ass_var: Optional[str] = None
    state_variable_qty: Optional[str] = None
    obs_type: str

    case_dir: Optional[str] = None
    case_emi_dir: Optional[str] = None
    emi_var: Optional[str] = None

    no_mems: int = Field(gt=0)
    run_assimilation_flag: bool = True
    update_restart: bool = None

    var_list_3d: Optional[list] = None
    var_list_2d: Optional[list] = None

    obs_filter_negative: Optional[bool] = None

    @field_validator("model_type", mode="before")
    @classmethod
    def normalize_model_type(cls, v):
        if isinstance(v, str):
            return v.lower()
        return v
    
    model_config = {
        "protected_namespaces": ()
    }

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
    obs_name: str
    vertical_ref_height: Optional[int] = None
    superobs: Optional[str] = None
    qa_value: Optional[float] = None

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
    is_control_ensemble: Optional[bool] = None
    perturbation_names: Optional[list] = None
    ensemble_list: Optional[list] = None
    dom_west: Optional[float] = None
    dom_east: Optional[float] = None
    dom_south: Optional[float] = None
    dom_north: Optional[float] = None
    nz: Optional[int] = None
    dlon: Optional[float] = None
    dlat: Optional[float] = None


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

    # runtime-only (not from YAML)
    _config_path: Optional[Path] = None

    model_config = {
        "protected_namespaces": (),
        "extra": "forbid"
    }


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

AppConfig.model_rebuild()
