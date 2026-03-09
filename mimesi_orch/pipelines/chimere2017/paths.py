from __future__ import annotations

from pathlib import Path
import pandas as pd

from pipelines.chimere2017.config import Chimere2017PipelineConfig
from mimesi_types import ModelType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from orchestrator_utils import TimeManager

class Chimere2017Paths:
    """Model-scoped path API for CHIMERE 2017."""

    def __init__(self, cfg: Chimere2017PipelineConfig):
        self.cfg = cfg
        self.base_path = cfg.base_path.resolve()
        self.env_python = self._resolve(cfg.env_python)
        self.listing_file = self._resolve(cfg.listing_file)
        self.path_submit_bsh = self._resolve(cfg.path_submit_bsh)
        self.path_filter = self._resolve(cfg.path_filter)
        self.path_data = self._resolve(cfg.path_data)
        self.path_control_run = None
        self.run_submit_model_template = (
            self._resolve(cfg.run_submit_model_template)
            if cfg.run_submit_model_template is not None
            else None
        )
        self.chimere_par_template = (
            self._resolve(cfg.chimere_par_template)
            if cfg.chimere_par_template is not None
            else None
        )

    def _resolve(self, p: Path) -> Path:
        if p.is_absolute():
            return p.resolve()
        return (self.base_path / p).resolve()

    # CHIMERE 2017
    def chimere_name_run_sub_ens_bash(self, string_to_replace_template: str) -> str:
        return f"run_mimesi_ens_{string_to_replace_template}.sh"

    def chimere_output_runs_dir(self, mem: int) -> Path:
        return self.path_data / f"runs_chimere/ITA7/RUNS_{mem}"

    # DART + S5P converter
    def dart_posteriors_dir(self, timestamp: str) -> Path:
        return self.path_data / f"posteriors/{timestamp}"

    def dart_s5p_base(self) -> Path:
        return self.base_path / "DART/observations/obs_converters/S5P_TROPOMI_L3"

    def dart_s5p_work(self) -> Path:
        return self.dart_s5p_base() / "work"

    def dart_s5p_input_template(self) -> Path:
        return self.dart_s5p_work() / "input_template.nml"

    def dart_s5p_input(self) -> Path:
        return self.dart_s5p_work() / "input.nml"

    def dart_s5p_data_dir(self) -> Path:
        return self.base_path / "SAT_OBS"

    def dart_file_s5p_orbit(self, orbit_filename: str) -> Path:
        return self.base_path() / orbit_filename

    def dart_s5p_output_dir(self) -> Path:
        return self.dart_s5p_data_dir() / "C03dart"

    def dart_obs_seq(self, obs_seq_name) -> Path:
        return self.dart_s5p_output_dir() / obs_seq_name

    def get_chimere_output_path(
        self,
        model: ModelType,
        mem: int,
        time_manager : TimeManager,
        prefix: str,
        offset : 1
        ) -> Path:
        if model != ModelType.CHIMERE:
            raise ValueError(f"Chimere2017Paths does not support model type: {model}")
        ts = time_manager.formatted_time(0, "%Y%m%d%H")
        
        ts_p1 = time_manager.formatted_time(offset, "%Y%m%d%H")

        return self.chimere_output_runs_dir(mem) / f"{prefix}.{ts}_{ts_p1}_ITA7.nc"
