from __future__ import annotations

from pathlib import Path
import pandas as pd

from pipelines.chimere2023.config import Chimere2023PipelineConfig
from mimesi_types import ModelType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pipeline_time import TimeManager

class Chimere2023Paths:
    """Model-scoped path API for CHIMERE 2023."""

    def __init__(self, cfg: Chimere2023PipelineConfig):
        self.cfg = cfg
        self.base_path_ctm = cfg.base_path_ctm.resolve()
        self.base_path_DART = cfg.base_path_DART.resolve()

        # Core Orchestrator Paths
        self.env_python = self._resolve(cfg.env_python)
        #self.case_name = self._resolve(cfg.case_name)
        self.listing_file = self._resolve(cfg.listing_file)
        self.path_data = self._resolve(cfg.path_data)
        self.path_filter = self._resolve(cfg.path_filter)
        
        # Source Model Paths
        self.path_control_run = self._resolve(cfg.path_control_run)
        self.path_perturbed_emi = self._resolve(cfg.path_perturbed_emi)
        self.path_perturbed_meteo = self._resolve(cfg.path_perturbed_meteo)

        self.path_data_dart: Path = self.path_data / 'OUT_DART' / cfg.case_name
        self.path_data_ctm: Path = self.path_data / 'OUT_Chimere' / cfg.case_name

    def _resolve(self, p: Path | None) -> Path:
        if p is None:
            return None
        if p.is_absolute():
            return p.resolve()
        return (self.base_path_DART / p).resolve()

    # ------------------------------------------------------------------
    # CHIMERE v2023 paths
    # ------------------------------------------------------------------

    def chimere2023_run_dir(self, mem: int) -> Path:
        return self.path_data_ctm / f"ENS{mem}"

    def chimere2023_WPS(self) -> Path:
        return self.path_control_run / "WPS" 
    
    def chimere2023_run_dir_WPS(self, mem: int) -> Path:
        return self.chimere2023_run_dir(mem) / "WPS"
    
    def chimere2023_ibc_src_dir(self) -> Path:
        return self.path_control_run / "IBC"
    
    def chimere2023_ibc_dir(self, mem: int) -> Path:
        return self.chimere2023_run_dir(mem) / "IBC"
    
    def chimere2023_BOUN_LIST_SRC(self, date_ymd: str, control_run_exp_name: str, domain: str) -> Path:
        return self.chimere2023_ibc_src_dir() / f"BOUN_CONCS.{date_ymd}00_24_{control_run_exp_name}_{domain}.list"

    def chimere2023_BOUN_LIST(self, date_ymdH: str, NHOURS: int, mem: int, domain: str) -> Path:
        return self.chimere2023_ibc_dir(mem) / f"BOUN_CONCS.{date_ymdH}_{NHOURS}_ENS{mem}_{domain}.list"
    
    def chimere2023_BOUN_NC(self, date_ymdH: str, NHOURS: int, mem: int, control_run_exp_name: str, domain: str) -> Path:
        return self.chimere2023_ibc_dir(mem) / f"BOUN_CONCS.{date_ymdH}_{NHOURS}_{control_run_exp_name}_{domain}.nc"
    
    def chimere2023_END_FILE_SRC(self, control_run_exp_name: str, date_ymd: str) -> Path:
        return self.path_control_run / f"end.{date_ymd}00_24_{control_run_exp_name}.nc"
    
    def chimere2023_END_FILE(self, mem: int, date_ymdH: str, NHOURS: int) -> Path:
        return self.chimere2023_run_dir(mem) / f"end.{date_ymdH}_{NHOURS}_ENS{mem}.nc"
    
    def chimere2023_PAR_BASE_TEMPLATE(self) -> Path:
        return self.base_path_ctm / "chimere.template_ensemble.par"
    
    def chimere2023_PAR_FILE(self, mem: int) -> Path:
        return self.base_path_ctm / f"chimere.ENS{mem}.par"
    
    def chimere2023_PAR_FILE_RUN_DIR(self, mem: int, START_DATEHOUR: str) -> Path:
        return self.chimere2023_run_dir(mem)/ f"chimere.ENS{mem}_{START_DATEHOUR}.par"
    
    def chimere2023_BASH_SUBMIT_SCRIPT_TEMPLATE(self) -> Path:
        return self.base_path_ctm / f"submit_p_template.sh"
    
    def chimere2023_BASH_SUBMIT_SCRIPT(self, mem: int) -> Path:
        return self.base_path_ctm / f"submit_p_{mem}.sh"
    
    def chimere2023_BASH_SUBMIT_SCRIPT_RUN_DIR(self, mem: int, START_DATEHOUR: str) -> Path:
        return self.chimere2023_run_dir(mem) / f"submit_p_{mem}_{START_DATEHOUR}.sh"

    def chimere2023_EMIS_FILE_SRC(self, is_pert: bool, domain: str, month: int, weekday: str, emi_id: Optional[int]) -> Path:
        if is_pert:
            return self.path_perturbed_emi / f"EMIS.{domain}.{month}.{weekday}.s.ens{emi_id:02d}.nc"
        else: 
            return self.path_control_run / f"EMIS.{domain}.{month}.{weekday}.s.nc"
    
    def chimere2023_EMI_FILE(self, mem: int, domain: str, month: int, weekday) -> Path:
        return self.chimere2023_run_dir(mem) / f"EMIS.{domain}.{month}.{weekday}.s.nc"
    
    def chimere2023_METEO_FILE_SRC(self, is_pert: bool, domain: str, date_ymd: str, meteo_id: Optional[int]) -> Path:
        if is_pert:
            return self.path_perturbed_meteo / f"exdomout_{date_ymd}00_24_{domain}.ens{meteo_id:02d}.nc"
        else:
            return self.path_control_run / f"exdomout_{date_ymd}00_24_{domain}.nc"
    
    def chimere2023_METEO_FILE(self, mem: int, domain: str, date_ymdH: str, NHOURS: int) -> Path:
        return self.chimere2023_run_dir(mem) / f"exdomout_{date_ymdH}_{NHOURS}_{domain}.nc"
    
    def chimere2023_out_file(self, mem: int, date_ymdH: str, NHOURS: int) -> Path:
        return self.chimere2023_run_dir(mem) / f"chim_ENS{mem}_{date_ymdH}_{NHOURS}_out.nc"

    def dart_filter_input_list_file(self, mem: int, date_ymdH: str) -> Path:
        return self.chimere2023_run_dir(mem) / f"chim_ENS{mem}_{date_ymdH}_1_out_toDART.nc"


    # ------------------------------------------------------------------
    # DART paths
    # ------------------------------------------------------------------

    def dart_posteriors_dir(self, date_ymdH: str) -> Path:
        return self.path_data_dart / f"posteriors/{date_ymdH}"

    def dart_analysis_dir(self, date_ymdH: str) -> Path:
        return self.path_data_dart / f"analysis/{date_ymdH}"

    def dart_preassim_dir(self, date_ymdH: str) -> Path:
        return self.path_data_dart / f"preassim/{date_ymdH}"
    
    def ratio_memory_file(self, mem: int, pollutant: str) -> Path:
        return self.path_data_dart / "ratio_memory" / f"ratio_memory_file_{pollutant}_ENS{mem}.nc"

    # ------------------------------------------------------------------
    # DART – filter
    # ------------------------------------------------------------------
    
    def dart_filter_input_template(self) -> Path:
        return self.base_path_DART / self.path_filter / "input_template.nml"
    
    def dart_filter_input(self) -> Path:
        return self.base_path_DART / self.path_filter / "input.nml"
    
    def dart_filter_input_list(self) -> Path:
        return self.base_path_DART / self.path_filter / "filter_input_list.txt"
    
    def dart_filter_output_list(self) -> Path:
        return self.base_path_DART / self.path_filter / "filter_output_list.txt"
    
    def dart_run_filter_template(self) -> Path:
        return self.base_path_DART / self.path_filter / "run_filter.template.bsh"
    
    def dart_run_filter(self) -> Path:
        return self.base_path_DART / self.path_filter / "run_filter.bsh"
    
    def dart_filter_output_list_file(self, mem: int, date_ymdH: str) -> Path:
        return self.dart_posteriors_dir(date_ymdH) / f"chim_ENS{mem}_{date_ymdH}_1_out_fromDART.nc"

    # ------------------------------------------------------------------
    # DART – obs converters (S5P)
    # ------------------------------------------------------------------

    def dart_s5p_base(self) -> Path:
        return self.base_path_DART / "observations/obs_converters/S5P_TROPOMI_L3"

    def dart_s5p_work(self) -> Path:
        return self.dart_s5p_base() / "work"

    def dart_s5p_input_template(self) -> Path:
        return self.dart_s5p_work() / "input_template.nml"

    def dart_s5p_input(self) -> Path:
        return self.dart_s5p_work() / "input.nml"

    def dart_s5p_data_dir(self, obs_name: str) -> Path:
        #return self.dart_s5p_base() / "data/SO2-COBRA"
        return self.dart_s5p_base() / "data" / obs_name
    
    def dart_file_s5p_orbit(self, orbit_filename: str, obs_name: str) -> Path:
        return self.dart_s5p_data_dir(obs_name) / orbit_filename
    
    def dart_file_s5p_orbit_filtered(self, orbit_filename: str, obs_name: str) -> Path:
        return self.dart_s5p_data_dir(obs_name) / "negative_filtered" / orbit_filename
    
    def dart_obs_seq(self, orbit_filename: str, obs_name: str, obs_seq_name: str) -> Path:
        return self.dart_file_s5p_orbit(orbit_filename, obs_name).parent / obs_seq_name
    
    def dart_obs_seq_filtered(self, orbit_filename: str, obs_name: str, obs_seq_name: str) -> Path:
        return self.dart_file_s5p_orbit_filtered(orbit_filename, obs_name).parent / obs_seq_name


    
    
 
   



