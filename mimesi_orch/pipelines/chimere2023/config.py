from __future__ import annotations

from pathlib import Path
from pydantic import BaseModel
from datetime import datetime

from config_models import AppConfig


class Chimere2023PipelineConfig(BaseModel):
    base_path_ctm: Path
    base_path_DART: Path

    env_python: Path
    listing_file: Path
    
    case_name: Path
    path_filter: Path
    path_data: Path
    path_control_run: Path
    path_perturbed_emi: Path | None = None
    path_perturbed_meteo: Path | None = None
    
    ass_var: str
    obs_type: str
    no_mems: int
    run_assimilation_flag: bool
    update_restart: bool
    var_list_3d: list | None = None
    var_list_2d: list | None = None
    obs_filter_negative: bool

    cluster_queue: str
    mail: str
    nb_proc: int
    project_name: str
    scheduler: str
    walltime: int

    chimpart: str
    control_run_exp_name: str
    dlat: float
    dlon: float
    dom_east: float
    dom_west: float
    dom_north: float
    dom_south:float
    nz: int
    domain: str
    ensemble_list: list
    perturbation_names: list
    restart_from_controlrun: bool
    submit_sequentially: bool

    obs_name: str
    qa_value: float
    search_window_seconds: int
    superobs: str
    vertical_ref_height: int

    @classmethod
    def from_app_config(cls, app: AppConfig) -> "Chimere2023PipelineConfig":
        return cls(
            base_path_ctm=app.paths.base_path_ctm,
            base_path_DART=app.paths.base_path_DART,

            env_python=app.paths.env_python,
            listing_file=app.paths.listing_file,
            
            path_filter=app.paths.path_filter,
            path_data=app.paths.path_data,
            path_control_run=app.paths.path_control_run,
            path_perturbed_emi=app.paths.path_perturbed_emi,
            path_perturbed_meteo=app.paths.path_perturbed_meteo,
            case_name=app.paths.case_name,
            
            ass_var=app.assimilation.ass_var,
            obs_type=app.assimilation.obs_type,
            no_mems=app.assimilation.no_mems,
            run_assimilation_flag=app.assimilation.run_assimilation_flag,
            update_restart=app.assimilation.update_restart,
            var_list_3d=app.assimilation.var_list_3d,
            var_list_2d=app.assimilation.var_list_2d,
            obs_filter_negative=app.assimilation.obs_filter_negative,

            cluster_queue=app.cluster.cluster_queue,
            mail=app.cluster.mail,
            nb_proc=app.cluster.nb_proc,
            project_name=app.cluster.project_name,
            scheduler=app.cluster.scheduler.value,
            walltime=app.cluster.walltime,

            chimpart=app.model_data.chimpart,
            control_run_exp_name=app.model_data.control_run_exp_name,
            dlat=app.model_data.dlat,
            dlon=app.model_data.dlon,
            dom_east=app.model_data.dom_east,
            dom_west=app.model_data.dom_west,
            dom_north=app.model_data.dom_north,
            dom_south=app.model_data.dom_south,
            nz=app.model_data.nz,
            domain=app.model_data.domain,
            ensemble_list=app.model_data.ensemble_list,
            perturbation_names=app.model_data.perturbation_names,
            restart_from_controlrun=app.model_data.restart_from_controlrun,
            submit_sequentially=app.model_data.submit_sequentially,

            obs_name=app.satellite_data.obs_name,
            qa_value=app.satellite_data.qa_value,
            search_window_seconds=app.satellite_data.search_window_seconds,
            superobs=app.satellite_data.superobs,
            vertical_ref_height=app.satellite_data.vertical_ref_height
            )
