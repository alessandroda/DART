from pathlib import Path
import logging
import pandas as pd
from config_models import PathsConfig
from mimesi_types import ModelType

logger = logging.getLogger(__name__)


class PathManager:
    """
    Centralized filesystem API for MIMESI.

    Responsibilities:
    - Resolve all paths from configuration
    - Validate required paths at startup
    - Provide semantic accessors (no string concatenation elsewhere)
    """

    def __init__(self, config: PathsConfig, *, log_paths: bool = True):
        self.cfg = config
        self.log_paths = log_paths

        # ---- resolve base path ----
        self.base_path: Path = config.base_path.resolve()

        # ---- resolve static paths ----
        self.env_python = self._resolve(config.env_python)
        self.listing_file = self._resolve(config.listing_file)
        self.run_submit_model_template = self._resolve(config.run_submit_model_template)
        self.path_submit_bsh = self._resolve(config.path_submit_bsh)
        self.path_filter = self._resolve(config.path_filter)
        self.path_data = self._resolve(config.path_data)
        self.run_submit_replace_perturbations = self._resolve(
            config.run_submit_replace_perturbations
        )
        self.chimere_par_template = (
                self._resolve(config.chimere_par_template)
                if config.chimere_par_template is not None
                else None
        )
        self._check_static_paths()

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------

    def _resolve(self, p: Path) -> Path:
        """Resolve a path relative to base_path."""
        return (self.base_path / p).resolve()

    def _check_static_paths(self):
        logger.info("Checking static paths")

        paths = {
            "base_path": self.base_path,
            "env_python": self.env_python,
            "listing_file": self.listing_file,
            "run_submit_model_template": self.run_submit_model_template,
            "path_submit_bsh": self.path_submit_bsh,
            "path_filter": self.path_filter,
            "path_data": self.path_data,
            "run_submit_replace_perturbations": self.run_submit_replace_perturbations,
            "chimere_par_template" : self.chimere_par_template
        }

        for name, path in paths.items():
            if path is None:
                continue
            if not path.exists():
                raise FileNotFoundError(f"[PathManager] {name} not found: {path}")
            if self.log_paths:
                logger.info(f"[PathManager] {name}: {path}")

    # ------------------------------------------------------------------
    # FARM paths (semantic API)
    # ------------------------------------------------------------------

    def farm_name_run_sub_ens_bash(self, string_to_replace_template) -> Path:
        return self.path_submit_bsh / f"run_submit_ens{string_to_replace_template}"

    def farm_output_dir(self, mem: int) -> Path:
        return self.path_data / f"OUTPUT_{mem}/OUT"

    def farm_ic(self, mem: int, timestamp: str) -> Path:
        return self.farm_output_dir(mem) / f"ic_g1_{timestamp}.nc"

    def farm_prior_dir(self, mem: int) -> Path:
        return self.farm_output_dir(mem) / "prior"

    def farm_output_runs_dir(self, mem: int, file_name: str) -> Path:
        return self.path_data / f"OUTPUT_{mem}/OUT/{file_name}"

    # ------------------------------------------------------------------
    # CHIMERE paths
    #
    def chimere_name_run_sub_ens_bash(self, string_to_replace_template) -> Path:
        return f"run_mimesi_ens_{string_to_replace_template}.sh"

    def chimere_output_runs_dir(self, mem: int) -> Path:
        return self.path_data / f"runs_chimere/ITA7/RUNS_{mem}"

    # ------------------------------------------------------------------
    # HERMES / emissions
    # ------------------------------------------------------------------

    def hermes_emission(self, mem: int, date_ymd00: str) -> Path:
        return self.path_data / f"INPUT/HERMES/emi_{mem}/HERMESv3_{date_ymd00}.nc"

    def hermes_emission_dir(self) -> Path:
        return self.path_data / "INPUT/HERMES/"

    # ------------------------------------------------------------------
    # DART paths
    # ------------------------------------------------------------------

    def dart_to_dart_dir(self) -> Path:
        return self.path_data / "to_DART"

    def dart_ic(self, mem: int, seconds: int, days: int) -> Path:
        return self.dart_to_dart_dir() / f"ic_g1_{seconds}_{days}_{mem}.nc"

    def dart_posteriors_dir(self, timestamp: str) -> Path:
        return self.path_data / f"posteriors/{timestamp}"

    def dart_analysis_dir(self, timestamp: str) -> Path:
        return self.path_data / f"analysis/{timestamp}"

    def dart_preassim_dir(self, timestamp: str) -> Path:
        return self.path_data / f"preassim/{timestamp}"

    # ------------------------------------------------------------------
    # DART – obs converters (S5P)
    # ------------------------------------------------------------------

    def dart_s5p_base(self) -> Path:
        return self.base_path / "SAT_OBS"

    def dart_s5p_work(self) -> Path:
        return self.base_path() / "DART/observations/obs_converters/S5P_TROPOMI_L3/"

    def dart_s5p_input_template(self) -> Path:
        return self.dart_s5p_work() / "input_template.nml"

    def dart_s5p_input(self) -> Path:
        return self.dart_s5p_work() / "input.nml"

    def dart_s5p_data_dir(self) -> Path:
        return self.dart_s5p_base() / "DART_obs"

    def dart_file_s5p_orbit(self, orbit_filename: str) -> Path:
        return self.dart_s5p_base() / orbit_filename

    def dart_s5p_output_dir(self) -> Path:
        return self.dart_s5p_data_dir() 

    def dart_obs_seq_name(
        self,
        *,
        sat_obs: pd.Timestamp | None = None,
        seconds: int | None = None,
        days: int | None = None,
    ) -> str:
        if sat_obs is not None:
            return f"obs_seq_{sat_obs.strftime('%Y%m%dT%H%M%S')}.out"
        if seconds is not None and days is not None:
            return f"obs_seq_{seconds}_{days}.out"
        raise ValueError("Provide sat_obs or both seconds and days for obs_seq naming")

    def dart_obs_seq(
        self,
        *,
        sat_obs: pd.Timestamp | None = None,
        seconds: int | None = None,
        days: int | None = None,
    ) -> Path:
        return self.dart_s5p_output_dir() / self.dart_obs_seq_name(
            sat_obs=sat_obs,
            seconds=seconds,
            days=days,
        )
    def get_ic_g1_path(
        self,
        model: ModelType,
        mem: int,
        timestamp: pd.Timestamp,
    ) -> Path:

        ts = timestamp.strftime("%Y%m%d%H")
        file_name = f"ic_g1_{ts}.nc"

        if model == ModelType.FARM:
            return self.farm_output_runs_dir(mem, file_name)

        elif model == ModelType.CHIMERE:
            return self.chimere_output_runs_dir(mem) / file_name

        else:
            raise ValueError(f"Unsupported model type: {model}")
