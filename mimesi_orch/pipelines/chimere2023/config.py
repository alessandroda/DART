from __future__ import annotations

from pathlib import Path
from pydantic import BaseModel

from config_models import AppConfig


class Chimere2023PipelineConfig(BaseModel):
    base_path: Path
    env_python: Path
    listing_file: Path
    path_submit_bsh: Path
    path_filter: Path
    path_data: Path
    path_control_run: Path
    path_perturbed_emi: Path | None = None
    path_perturbed_meteo: Path | None = None

    ass_var: str
    state_variable_qty: str
    obs_type: str
    case_dir: str
    case_emi_dir: str | None = None
    emi_var: str | None = None
    no_mems: int
    run_assimilation_flag: bool = True

    cluster_queue: str
    project_name: str
    nb_proc: int | None = None
    walltime: int | None = None
    mail: str | None = None

    backup_perturb_days: int | None = None
    backup_ic_hours: int | None = None
    backup_ic_option: str | None = None

    control_run_exp_name: str | None = None
    domain: str | None = None
    chimpart: str | None = None
    submit_sequentially: bool | None = None
    restart_from_controlrun: bool | None = None
    perturbation_names: list | None = None
    ensemble_list: list | None = None

    @classmethod
    def from_app_config(cls, app: AppConfig) -> "Chimere2023PipelineConfig":
        if app.paths.path_control_run is None:
            raise ValueError("chimere2023 pipeline requires paths.path_control_run")
        return cls(
            base_path=app.paths.base_path,
            env_python=app.paths.env_python,
            listing_file=app.paths.listing_file,
            path_submit_bsh=app.paths.path_submit_bsh,
            path_filter=app.paths.path_filter,
            path_data=app.paths.path_data,
            path_control_run=app.paths.path_control_run,
            path_perturbed_emi=app.paths.path_perturbed_emi,
            path_perturbed_meteo=app.paths.path_perturbed_meteo,
            ass_var=app.assimilation.ass_var,
            state_variable_qty=app.assimilation.state_variable_qty,
            obs_type=app.assimilation.obs_type,
            case_dir=app.assimilation.case_dir,
            case_emi_dir=app.assimilation.case_emi_dir,
            emi_var=app.assimilation.emi_var,
            no_mems=app.assimilation.no_mems,
            run_assimilation_flag=app.assimilation.run_assimilation_flag,
            cluster_queue=app.cluster.cluster_queue,
            project_name=app.cluster.project_name,
            nb_proc=app.cluster.nb_proc,
            walltime=app.cluster.walltime,
            mail=app.cluster.mail,
            backup_perturb_days=app.time.backup_perturb_days,
            backup_ic_hours=app.time.backup_ic_hours,
            backup_ic_option=app.time.backup_ic_option,
            control_run_exp_name=app.model_data.control_run_exp_name if app.model_data else None,
            domain=app.model_data.domain if app.model_data else None,
            chimpart=app.model_data.chimpart if app.model_data else None,
            submit_sequentially=app.model_data.submit_sequentially if app.model_data else None,
            restart_from_controlrun=app.model_data.restart_from_controlrun if app.model_data else None,
            perturbation_names=app.model_data.perturbation_names if app.model_data else None,
            ensemble_list=app.model_data.ensemble_list if app.model_data else None,
        )
