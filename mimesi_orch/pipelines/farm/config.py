from __future__ import annotations

from pathlib import Path
from pydantic import BaseModel

from config_models import AppConfig


class FarmPipelineConfig(BaseModel):
    base_path: Path
    env_python: Path
    listing_file: Path
    run_submit_model_template: Path | None
    path_submit_bsh: Path
    path_filter: Path
    path_data: Path
    run_submit_replace_perturbations: Path | None

    ass_var: str
    state_variable_qty: str
    obs_type: str
    case_dir: str
    case_emi_dir: str | None = None
    emi_var: str | None = None
    no_mems: int
    run_assimilation_flag: bool = True

    cluster_queue: str

    backup_perturb_days: int | None = None
    backup_ic_hours: int | None = None
    backup_ic_option: str | None = None

    @classmethod
    def from_app_config(cls, app: AppConfig) -> "FarmPipelineConfig":
        return cls(
            base_path=app.paths.base_path,
            env_python=app.paths.env_python,
            listing_file=app.paths.listing_file,
            run_submit_model_template=app.paths.run_submit_model_template,
            path_submit_bsh=app.paths.path_submit_bsh,
            path_filter=app.paths.path_filter,
            path_data=app.paths.path_data,
            run_submit_replace_perturbations=app.paths.run_submit_replace_perturbations,
            ass_var=app.assimilation.ass_var,
            state_variable_qty=app.assimilation.state_variable_qty,
            obs_type=app.assimilation.obs_type,
            case_dir=app.assimilation.case_dir,
            case_emi_dir=app.assimilation.case_emi_dir,
            emi_var=app.assimilation.emi_var,
            no_mems=app.assimilation.no_mems,
            run_assimilation_flag=app.assimilation.run_assimilation_flag,
            cluster_queue=app.cluster.cluster_queue,
            backup_perturb_days=app.time.backup_perturb_days,
            backup_ic_hours=app.time.backup_ic_hours,
            backup_ic_option=app.time.backup_ic_option,
        )
