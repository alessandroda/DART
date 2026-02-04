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
        self.path_control_run = self._resolve(config.path_control_run)
        self.path_perturbed_emi = self._resolve(config.path_perturbed_emi)
        self.path_perturbed_meteo = self._resolve(config.path_perturbed_meteo)

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
            "path_control_run": self.path_control_run,
            "path_perturbed_emi": self.path_perturbed_emi,
            "path_perturbed_meteo": self.path_perturbed_meteo
        }

        for name, path in paths.items():
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
    # CHIMERE v2017 paths
    # ------------------------------------------------------------------
    def chimere_name_run_sub_ens_bash(self, string_to_replace_template) -> Path:
        return (
            self.path_submit_bsh / f"run_mimesi_ens_member_{string_to_replace_template}"
        )

    def chimere_output_runs_dir(self, mem: int) -> Path:
        return self.path_data / f"runs_chimere/ITA7/RUNS_{mem}"
    
    # ------------------------------------------------------------------
    # CHIMERE v2023 paths
    # ------------------------------------------------------------------

    def chimere2023_run_dir(self, mem: int) -> Path:
        return self.path_data / f"ENS{mem}"

    def chimere2023_WPS(self) -> Path:
        return self.path_control_run / "WPS" 
    
    def chimere2023_run_dir_WPS(self, mem: int) -> Path:
        return self.chimere2023_run_dir(mem) / "WPS"
    
    def chimere2023_ibc_src_dir(self) -> Path:
        return self.path_control_run / "IBC"
    
    def chimere2023_ibc_dir(self, mem: int) -> Path:
        return self.chimere2023_run_dir(mem) / "IBC"
    
    def chimere2023_BOUN_LIST_SRC(self, date_ymd: str, control_run_exp_name: str) -> Path:
        return self.chimere2023_ibc_src_dir() / f"BOUN_CONCS.{date_ymd}00_24_{control_run_exp_name}_EUROCOMEX3.list"

    def chimere2023_BOUN_LIST(self, date_ymd: str, mem: int) -> Path:
        return self.chimere2023_ibc_dir(mem) / f"BOUN_CONCS.{date_ymd}00_24_ENS{mem}_EUROCOMEX3.list"
    
    def chimere2023_END_FILE_SRC(self, control_run_exp_name: str, date_ymd: str) -> Path:
        return self.path_control_run / f"end.{date_ymd}00_24_{control_run_exp_name}.nc"
    
    def chimere2023_END_FILE(self, mem: int, date_ymd: str) -> Path:
        return self.chimere2023_run_dir(mem) / f"end.{date_ymd}00_24_ENS{mem}.nc"
    
    def chimere2023_PAR_BASE_TEMPLATE(self) -> Path:
        return self.base_path / "chimere.ensemble.par"
    
    def chimere2023_PAR_FILE(self, mem: int) -> Path:
        return self.base_path / f"chimere.ENS{mem}.par"
    
    def chimere2023_PAR_FILE_RUN_DIR(self, mem: int) -> Path:
        return self.chimere2023_run_dir(mem)/ f"chimere.ENS{mem}.par"
    
    def chimere2023_EMIS_FILE_SRC(self, is_pert: bool, domain: str, month: int, weekday: str, emi_id: int) -> Path:
        if is_pert:
            return self.path_perturbed_emi / f"EMIS.{domain}.{month}.{weekday}.s.ens0{emi_id}.nc"
        else: 
            return self.path_control_run / f"EMIS.{domain}.{month}.{weekday}.s.nc"
    
    def chimere2023_EMI_FILE(self, mem: int, domain: str, month: int, weekday) -> Path:
        return self.chimere2023_run_dir(mem) / f"EMIS.{domain}.{month}.{weekday}.s.nc"
    
    def chimere2023_METEO_FILE_SRC(self, is_pert: bool, domain: str, date_ymd: str, meteo_id: int) -> Path:
        if is_pert:
            return self.path_perturbed_meteo / f"exdomout_{date_ymd}00_24_{domain}.ens0{meteo_id}.nc"
        else:
            return self.path_control_run / f"exdomout_{date_ymd}00_24_{domain}.nc"
    
    def chimere2023_METEO_FILE(self, mem: int, domain: str, date_ymd: str) -> Path:
        return self.chimere2023_run_dir(mem) / f"exdomout_{date_ymd}00_24_{domain}.nc"
    
    def chimere2023_BASH_SUBMIT_SCRIPT(self, mem: int) -> Path:
        return self.base_path / f"submit_p_{mem}.sh"
    


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
        return self.base_path / "DART/observations/obs_converters/S5P_TROPOMI_L3"

    def dart_s5p_work(self) -> Path:
        return self.dart_s5p_base() / "work"

    def dart_s5p_input_template(self) -> Path:
        return self.dart_s5p_work() / "input_template.nml"

    def dart_s5p_input(self) -> Path:
        return self.dart_s5p_work() / "input.nml"

    def dart_s5p_data_dir(self) -> Path:
        return self.dart_s5p_base() / "data/SO2-COBRA"

    def dart_file_s5p_orbit(self, orbit_filename: str) -> Path:
        return self.dart_s5p_data_dir() / orbit_filename

    def dart_s5p_output_dir(self) -> Path:
        return self.dart_s5p_data_dir() / "C03dart"

    def dart_obs_seq(self, seconds: int, days: int) -> Path:
        return self.dart_s5p_output_dir() / f"obs_seq_{seconds}_{days}.out"

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
