from datetime import timedelta
from pathlib import Path
import shutil
import time
import numpy as np
import re
from datetime import datetime
from pipeline_errors import FatalPipelineError, ModelRunError
from mimesi_types import Scheduler
from config_models import AppConfig
from pipelines.base_pipeline import BaseAssimilationPipeline
import logging
import os
import pandas as pd
import subprocess
from pipelines.chimere2023.paths import Chimere2023Paths
import xarray as xr
from time_utils import set_date_gregorian
from orchestrator_utils import (
    CommandSpec,
    modify_yaml_date,
    replace_nml_template,
    searchFile,
    run_command_in_directory,
    safe_symlink,
    check_and_clean_broken_links,
    from_liststr_to_listdict,
    submit_irene,
    check_restart_files_exist_irene,
    cut_block,
    add_missing_variable,
    write_dart_filter_list,
    update_pollutant_in_end,
    save_diff,
    remove_negative_values,
    monitor_job_status
)
from pipeline_time import AssimWindow, TimeManager


logger = logging.getLogger(__name__)

class Chimere2023DartPipeline(BaseAssimilationPipeline):
    def __init__(
        self,
        time_manager: TimeManager,
        paths: Chimere2023Paths,
        config: AppConfig,
    ):
        super().__init__(time_manager)
        # Load configuration from YAML file
        self.config = config
        # Setup path manager and time manager with loaded configurations
        self.paths = paths
        self.listing = pd.read_csv(self.paths.listing_file, sep=";")
        self.listing["start_time"] = pd.to_datetime(self.listing["start_time"])
        self.days_obs = 0
        self.seconds_obs = 0
        self.days_model = 0
        self.seconds_model = 0
        self.output_sim_folder = None
        a = self.config.assimilation
        self.current_end_file = {}
        self.satdata_found = False

        self.model_type = a.model_type
        logger.info(f"Running assimilation with model_type={self.model_type}")

        self.ass_var = a.ass_var
        self.emi_var = a.emi_var
        self.no_mems = a.no_mems
        self.obs_type = a.obs_type
        self.run_assimilation_flag = a.run_assimilation_flag
        self.update_restart = a.update_restart
        self.case_emi_dir = a.case_emi_dir
        self.var_list_3d = a.var_list_3d
        self.var_list_2d = a.var_list_2d
        self.obs_filt = a.obs_filter_negative

        self.queue = self.config.cluster.cluster_queue
        self.project_name = self.config.cluster.project_name
        self.nb_proc = self.config.cluster.nb_proc
        self.walltime = self.config.cluster.walltime
        self.mail = self.config.cluster.mail

        self.ensemble_list = self.config.model_data.ensemble_list
        self.perturbation_names = self.config.model_data.perturbation_names
        self.control_run_exp_name = self.config.model_data.control_run_exp_name
        self.domain = self.config.model_data.domain
        self.chimpart = self.config.model_data.chimpart
        self.submit_sequentially = self.config.model_data.submit_sequentially
        self.restart_from_controlrun = self.config.model_data.restart_from_controlrun
        self.is_control_ensemble = self.config.model_data.is_control_ensemble
        self.dom_west = self.config.model_data.dom_west
        self.dom_east = self.config.model_data.dom_east
        self.dom_south = self.config.model_data.dom_south
        self.dom_north = self.config.model_data.dom_north
        self.nz = self.config.model_data.nz
        self.dlon = self.config.model_data.dlon
        self.dlat = self.config.model_data.dlat

        self.search_window_seconds = self.config.satellite_data.search_window_seconds
        self.obs_name = self.config.satellite_data.obs_name
        self.vertical_ref_height = self.config.satellite_data.vertical_ref_height
        self.superobs = self.config.satellite_data.superobs
        self.qa_value = self.config.satellite_data.qa_value


        self.scheduler = self.config.cluster.scheduler
        logger.info(f"Using scheduler={self.scheduler}, queue={self.queue}")
        self._pending_orbit_filename = None
        self._pending_orbit_start = None
        self._pending_orbit_time = None
        self._generated_daily_emission_files: list[Path] = []
        self._generated_daily_emission_stamp: str | None = None
        self._previous_window_start: pd.Timestamp | None = None
        self._previous_window_end: pd.Timestamp | None = None
        self._cleanup_completed_cycles: int = 0

    def build_assim_window(self) -> AssimWindow:
        """
        Build the current cycle window.

        Windows start and end on model hour boundaries. Observations are
        assigned to the nearest cycle end using the +/- 30 minute rule.
        """
        start_time = self.time_manager.current_time
        
        if start_time.hour == 0:
            return AssimWindow(
                start_time=start_time,
                end_time=start_time + self.time_manager.dt,
                run_hours=1,
                has_assimilation=False,   # optional
            )
        half_dt = self.time_manager.dt / 2
        day_end = start_time.replace(
            hour=0,
            minute=0,
            second=0,
            microsecond=0,
        ) + timedelta(days=1)
        run_limit = min(day_end, self.time_manager.end_time + self.time_manager.dt)

        slot_time = start_time + self.time_manager.dt
        while slot_time <= run_limit:
            orbit_matches = searchFile(
                slot_time - half_dt,
                slot_time + half_dt,
                self.listing,
            )
            if not orbit_matches.empty:
                orbit_matches = orbit_matches.copy()
                orbit_matches["start_time"] = pd.to_datetime(
                    orbit_matches["start_time"]
                )
                orbit_matches = orbit_matches.sort_values("start_time")
                orbit_row = orbit_matches.iloc[0]
                run_hours = int((slot_time - start_time) / self.time_manager.dt)
                return AssimWindow(
                    start_time=start_time,
                    end_time=slot_time,
                    run_hours=run_hours,
                    has_assimilation=True,
                    obs_time=pd.to_datetime(orbit_row["start_time"]),
                    orbit_filename=orbit_row["filename"],
                )
            slot_time += self.time_manager.dt

        if run_limit <= start_time:
            raise FatalPipelineError(
                f"Invalid run limit for assimilation window: start={start_time} end={run_limit}"
            )

        return AssimWindow(
            start_time=start_time,
            end_time=run_limit,
            run_hours=int((run_limit - start_time) / self.time_manager.dt),
            has_assimilation=False,
        )

    def before_step(self):
        dict_mem_list = from_liststr_to_listdict(self.ensemble_list, self.perturbation_names)
        if len(dict_mem_list) != self.no_mems:
            logger.info(f"NB The ensemble list values ({len(dict_mem_list)}) are not equal in number to the requested assimilation ensemble ({self.no_mems})")

        for dict_mem in dict_mem_list: #loop on mems
            self.paths.chimere2023_run_dir(dict_mem["MemberID"]).mkdir(parents=True, exist_ok=True)
            logger.info(f"Creating directories and links for ENS{dict_mem['MemberID']} to run chimere's parallel part")
            #link BOUN, EMIS and METEO (from either control run or perturbed dir))
            
            date_month = self.time_manager.current_time.strftime("%m")
            date_weekday = self.time_manager.current_time.strftime("%A")
            logger.info("Linking EMIS ...")
            safe_symlink(self.paths.chimere2023_EMIS_FILE_SRC("EmisID" in dict_mem.keys(), self.domain, date_month, date_weekday, dict_mem["EmisID"]),
                            self.paths.chimere2023_EMI_FILE(dict_mem["MemberID"], self.domain, date_month, date_weekday))
            if self.time_manager.slot_time.strftime("%H") == '00':
                safe_symlink(self.paths.chimere2023_EMIS_FILE_SRC("EmisID" in dict_mem.keys(), self.domain, date_month, (self.time_manager.current_time + timedelta(days=1)).strftime("%A"), dict_mem["EmisID"]),
                            self.paths.chimere2023_EMI_FILE(dict_mem["MemberID"], self.domain, date_month, (self.time_manager.current_time + timedelta(days=1)).strftime("%A")))
            if self.restart_from_controlrun and self.time_manager.current_time == self.time_manager.start_time:
                #link first end file
                logger.info("Linking END ...")
                end_file_date_control_run = self.time_manager.start_time - timedelta(days=1)
                safe_symlink(self.paths.chimere2023_END_FILE_SRC(self.control_run_exp_name, end_file_date_control_run.strftime("%Y%m%d")), 
                             self.paths.chimere2023_END_FILE(dict_mem["MemberID"], end_file_date_control_run.strftime("%Y%m%d00"), 24))
                self.current_end_file[dict_mem["MemberID"]]=self.paths.chimere2023_END_FILE(dict_mem["MemberID"], end_file_date_control_run.strftime("%Y%m%d00"), 24)
            else:
                self.current_end_file[dict_mem["MemberID"]]=None
            if check_and_clean_broken_links(self.paths.chimere2023_run_dir(dict_mem["MemberID"])):
                raise FatalPipelineError(f"Broken links detected. Please clean up. The target needs to exists.")
            else:
                logger.info(f">> All links are good for ENS{dict_mem['MemberID']}  ...")
            
            #date_ymd00 = self.time_manager.current_time.strftime("%Y%m%d00")
            #logger.info("Linking METEO ...")
            #safe_symlink(self.paths.chimere2023_METEO_FILE_SRC("MeteoID" in dict_mem.keys(), self.domain, date_ymd, dict_mem["MeteoID"]), 
            #                self.paths.chimere2023_METEO_FILE(dict_mem["MemberID"], self.domain, date_ymd00, 24))
            #logger.info("Linking BOUN ...")
            #safe_symlink(self.paths.chimere2023_BOUN_LIST_SRC(date_ymd, self.control_run_exp_name, self.domain), 
            #                self.paths.chimere2023_BOUN_LIST(date_ymd00, 24, dict_mem["MemberID"], self.domain))
           
            logger.info("Computing BOUNs and exdomouts for the specific hours to run ...")
            date_ymd = self.time_manager.current_time.strftime("%Y%m%d")
            
            boun_src = self.paths.chimere2023_BOUN_LIST_SRC(date_ymd, self.control_run_exp_name, self.domain).read_text().splitlines()[1]
            meteo_src = self.paths.chimere2023_METEO_FILE_SRC("MeteoID" in dict_mem.keys(), self.domain, date_ymd, dict_mem["MeteoID"] if "MeteoID" in dict_mem.keys() else None)
            
            with xr.open_dataset(boun_src) as ds_boun, xr.open_dataset(str(meteo_src)) as ds_meteo:
                    date_ymdH = self.time_manager.current_time.strftime("%Y%m%d%H")
                    date_H = self.time_manager.current_time.strftime("%H")

                    save_path_boun = self.paths.chimere2023_BOUN_NC(date_ymdH, self.time_manager.run_hours, dict_mem["MemberID"], self.control_run_exp_name, self.domain)
                    save_path_boun_list = self.paths.chimere2023_BOUN_LIST(date_ymdH, self.time_manager.run_hours, dict_mem["MemberID"], self.domain)
                    cut_block(ds_boun, int(date_H), self.time_manager.run_hours, save_path_boun, save_path_boun_list)

                    save_path_meteo = self.paths.chimere2023_METEO_FILE(dict_mem["MemberID"], self.domain, date_ymdH, self.time_manager.run_hours)
                    cut_block(ds_meteo, int(date_H), self.time_manager.run_hours, save_path_meteo)
                   
                    #filter negative values in obs e magari anche di piu, c'é da vedere PUM

    def run_model(self):
        """
        Run CHIMERE for the current time step: from current_time to current_time + 1 hour.
        """
        logger.info(f"---------->>> Running CHIMERE model from {self.time_manager.current_time} to {self.time_manager.slot_time}")
        date_ymdH = self.time_manager.current_time.strftime("%Y%m%d%H")

        job_ids = []
        for mem in range(1, self.no_mems + 1):
            try:
                logger.info("Replacing @TOKENS in CHIMERE .par template file ...")
                logger.info(f"The output directory (run_dir) is: {self.paths.chimere2023_run_dir(mem)}")
                if self.current_end_file[mem] is None:
                    if self.is_control_ensemble and (self.time_manager.current_time == self.time_manager.start_time):
                        endfile_f = self.paths.chimere2023_END_FILE(mem, self.time_manager.end_file_datetime.strftime("%Y%m%d%H")[:-2]+"00", 24) 
                    elif self.time_manager.prev_run_hours is not None:
                        endfile_f = self.paths.chimere2023_END_FILE(mem, self.time_manager.end_file_datetime.strftime("%Y%m%d%H"), self.time_manager.prev_run_hours)
                    else:
                        endfile_f = self.paths.chimere2023_END_FILE(mem, "2020020613", 1)

                else:
                    endfile_f = self.current_end_file[mem]
                logger.info(f"The END file used for ENS{mem} is: {endfile_f}")
                replace_nml_template(
                    input_nml_path=self.paths.chimere2023_PAR_BASE_TEMPLATE(),
                    entries_tbr_dict={
                        "@LAB": f"ENS{mem}",
                        "@SIMULDIR": self.paths.chimere2023_run_dir(mem),
                        "@IUSEINI": "2",
                        "@ENDFILE": endfile_f,
                        "@EMISSDIR": self.paths.chimere2023_run_dir(mem)
                    },
                    output_nml_path=self.paths.chimere2023_PAR_FILE(mem)
                )
                shutil.copy(self.paths.chimere2023_PAR_FILE(mem), self.paths.chimere2023_PAR_FILE_RUN_DIR(mem, date_ymdH))
            except Exception as e:
                raise FatalPipelineError(f"Failed to prepare CHIMERE .par file: {e}")
            try:
                logger.info("Replacing @TOKENS in CHIMERE template submit script ...")
                replace_nml_template(
                    input_nml_path=self.paths.chimere2023_BASH_SUBMIT_SCRIPT_TEMPLATE(),
                    entries_tbr_dict={
                        "@JOBNAME": f"ENS{mem}",
                        "@NPROC": f"{self.nb_proc}",
                        "@QUEUED_NODES": f"{self.queue}",
                        "@PROJECT": f"{self.project_name}",
                        "@WALLTIME": f"{self.walltime}",
                        "@RUN_DIR": f"{self.paths.chimere2023_run_dir(mem)}",
                        "@MAIL": f"{self.mail}",
                        "@PARFILE": f"{self.paths.chimere2023_PAR_FILE(mem).name}",
                        "@START_DATEHOUR": f"{date_ymdH}",
                        "@NHOURS_forward": f"{self.time_manager.run_hours}",
                        "@CHIMPART": f"{self.chimpart}"
                    },
                    output_nml_path=self.paths.chimere2023_BASH_SUBMIT_SCRIPT(mem)
                )
                shutil.copy(self.paths.chimere2023_BASH_SUBMIT_SCRIPT(mem), self.paths.chimere2023_BASH_SUBMIT_SCRIPT_RUN_DIR(mem, date_ymdH))
                subprocess.run(["chmod", "+x", str(self.paths.chimere2023_BASH_SUBMIT_SCRIPT(mem))], check=True)
            except Exception as e:
                raise FatalPipelineError(f"Failed to prepare CHIMERE submit script: {e}")
             
            logger.info(f"Queuing job for member {mem}...")
            if len(job_ids)==0 or not self.submit_sequentially:
                submit_command = f"ccc_msub ./{self.paths.chimere2023_BASH_SUBMIT_SCRIPT(mem).name}"
            else:
                submit_command = f"ccc_msub -a {job_ids[-1]} ./{self.paths.chimere2023_BASH_SUBMIT_SCRIPT(mem).name}"

            job_ids.append(submit_irene(CommandSpec(
                command=submit_command,
                directory=self.paths.base_path_ctm
            )))
            
        logger.info(f"Checking job status ...")
        monitor_job_status(job_ids, self.scheduler, self.model_type)
        
        logger.info(f"Checking restart files were created ...")
        date_ymdH = self.time_manager.current_time.strftime("%Y%m%d%H")
        mems_to_rerun = check_restart_files_exist_irene(
            ic_paths = [self.paths.chimere2023_END_FILE(memb, date_ymdH, self.time_manager.run_hours) for memb in range(1, self.no_mems + 1)], #self.paths.chimere2023_END_FILE(mem, date_ymdH, 1),
            model=self.model_type,
            no_mems=self.no_mems,
            ) #if old and not updated it returns a false positive
        if mems_to_rerun:
            logger.info(f"Check chimere log file at: {self.paths.chimere2023_run_dir(mem)}/ENS{mem}_{date_ymdH}.out")
            raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
        logger.info(f" Run_model() completed successfully.")

    def after_model(self):
        logger.info("Saving chimere's output files title timestamp (that is the starting time of the run)")
        self.time_manager.end_file_datetime = self.time_manager.current_time #(equivalent to self.time_manager.simulated_time - timedelta(hours=1))

    def run_assimilation_if_needed(self):
        """
        Run DART assimilation only if:
        - assimilation is enabled
        - satellite observations are available
        """
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled by config")
            return
        
        logger.info(f"---------->>> Running process_satellite_data()")
        orbit_filename = self.process_satellite_data()
        if not orbit_filename:
            logger.info("[DART] No satellite data found, skipping assimilation")
            return

        logger.info(f"---------->>> Running run_obs_converter()")
        obs_seq_name = self.run_obs_converter(orbit_filename)

        obs_path = self.paths.dart_obs_seq(orbit_filename, self.obs_name, obs_seq_name)
        if not obs_path.exists():
            logger.info("[DART] convertion was run but all obs were excluded and obs_seq was not created, skipping assimilation")
            return
        else:
            logger.info(f"[DART] obs_seq created: {obs_path}")
        self.satdata_found = True
        self.run_dart(obs_seq_name, obs_path)

    def process_satellite_data(self):

        orbit_filename = searchFile(
            self.time_manager.current_time - 0.5 * self.time_manager.dt,
            self.time_manager.current_time + 0.5 * self.time_manager.dt,
            self.listing,
        )
        if orbit_filename.empty:
            return False

        orbit_filename["start_time"] = pd.to_datetime(orbit_filename["start_time"])
        
        #orbit_filename = orbit_filename[orbit_filename["start_time"].dt.hour >= 10]
        #if orbit_filename.empty:
        #    logger.info(f'No valid orbit file found after 10 AM.')
        #    return False

        logger.info(f"Orbit file found: {orbit_filename['filename'].values[0]}")
        self.time_manager.sat_obs = pd.to_datetime(
            orbit_filename["start_time"].values[0]
        )
        return orbit_filename["filename"].values[0]

    def run_obs_converter(self, orbit_filename):
        self.seconds_obs, self.days_obs = set_date_gregorian(
            self.time_manager.sat_obs.year,
            self.time_manager.sat_obs.month,
            self.time_manager.sat_obs.day,
            self.time_manager.sat_obs.hour,
            self.time_manager.sat_obs.minute,
            self.time_manager.sat_obs.second,
        )
        obs_seq_name = f"obs_seq_{self.seconds_obs}_{self.days_obs}.out"
        obs_seq_path = self.paths.dart_obs_seq(orbit_filename, self.obs_name, obs_seq_name)
        
        # Skip submission if file already exists
        if Path(obs_seq_path).exists():
            logger.info(f"Obs sequence file already exists: {obs_seq_path}")
            return obs_seq_name
        
        #filter negative values in satellite observations if requested (and save the filtered file in a different location to keep the original data)
        if self.obs_filt:
            logger.info(f"---------->>> Removing negative values in satellite observations ...")
            remove_negative_values(self.paths.dart_file_s5p_orbit(orbit_filename, self.obs_name), 
                                   self.paths.dart_file_s5p_orbit_filtered(orbit_filename, self.obs_name))
        logger.info(f"CHECK file_path_s5p: {self.paths.dart_file_s5p_orbit_filtered(orbit_filename, self.obs_name)}") if self.obs_filt else logger.info(f"file_path_s5p: {self.paths.dart_file_s5p_orbit(orbit_filename, self.obs_name)}")
        
        
        replace_nml_template(
            input_nml_path=self.paths.dart_s5p_input_template(),
            entries_tbr_dict={
                "$file_path_s5p": self.paths.dart_file_s5p_orbit_filtered(orbit_filename, self.obs_name) if self.obs_filt else self.paths.dart_file_s5p_orbit(orbit_filename, self.obs_name),
                "$file_out": obs_seq_path,
                "$obs_type": self.obs_type,
                "$dom_west": self.dom_west,
                "$dom_east": self.dom_east,
                "$dom_south": self.dom_south,
                "$dom_north": self.dom_north,
                "$nz": self.nz,
                "$dlon": self.dlon,
                "$dlat": self.dlat,
                "$vertical_ref_height": self.vertical_ref_height,
                "$superobs": self.superobs,
                "$qa_value": self.qa_value
            },
            output_nml_path=self.paths.dart_s5p_input(),
        )

        try:
            spec = CommandSpec(
                command="convert_s5p_tropomi_l3",
                directory=self.paths.dart_s5p_work(),
                )
            rc, _ = run_command_in_directory(spec)
            if rc != 0:
                raise SchedulerError(
                    f"Submission command failed: {spec.command} " f"(return code {rc})"
                )
        except Exception as e:
            logger.error(f"Error running obs converter: {e}")
            return False
        return obs_seq_name

    def run_dart(self, obs_seq_name, obs_path):
        logger.info("---------->>> Running DART")
        """
        All chimere files in the title have the starting time of the run and now current_time=simulated_time=final simulation time
        """
        date_ymdH = self.time_manager.end_file_datetime.strftime("%Y%m%d%H")
        date_ymdHMS = self.time_manager.end_file_datetime.strftime("%Y%m%d_%H%M%S") #tbm it has 00:00 always because it computed from current_time; 
        # It is used in the obs_seq.final name: need to decide what to put; either chim out that has HH:30, or the satellite obs first aquisition time or the restart time)
        logger.info(f"The timestamp in DART results' titles follows chimere's logic: it is used the simulation starting time (that was previously saved and is {date_ymdH}), not the simulated_time")

        self.paths.dart_posteriors_dir(date_ymdH).mkdir(parents=True, exist_ok=True)

        add_missing_variable(
            no_mems=self.no_mems, 
            var_to_add='psfc', 
            domain=self.domain, 
            out_file_func=self.paths.chimere2023_out_file,
            orig_file_func=self.paths.chimere2023_METEO_FILE,
            date_ymdH=date_ymdH, 
            NHOURS=self.time_manager.run_hours
        )

        replace_nml_template(
            input_nml_path=self.paths.dart_filter_input_template(),
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": self.paths.dart_posteriors_dir(date_ymdH),
                "$folder_obs_path": obs_path.parent if not self.obs_filt else obs_path.parent,
                "$date_assim": date_ymdHMS,
                "$init_time_days": str(self.days_model), #computed in base_pipeline
                "$init_time_seconds": str(self.seconds_model), #computed in base_pipeline
                "$first_obs_days": str(self.days_obs), #computed when creating obs_seq.out 
                "$first_obs_seconds": str(self.seconds_obs), #computed when creating obs_seq.out
                "$no_mems": str(self.no_mems),
                "$obs_type": str(self.obs_type),
                "$num_3d": len(self.var_list_3d),
                "$list_3d": ', '.join("'" + v + "'" for v in self.var_list_3d),
                "$list_2d": ', '.join("'" + v + "'" for v in self.var_list_2d),
                "$num_2d": len(self.var_list_2d),
                
            },
            output_nml_path=self.paths.dart_filter_input(),
        )

        write_dart_filter_list(
            list_file_func=self.paths.dart_filter_input_list(), 
            out_file_func=self.paths.chimere2023_out_file, 
            no_mems=self.no_mems, 
            date_ymdH=date_ymdH, 
            NHOURS=self.time_manager.run_hours)

        write_dart_filter_list(
            list_file_func=self.paths.dart_filter_output_list(), 
            out_file_func=self.paths.dart_filter_output_list_file, 
            no_mems=self.no_mems, 
            date_ymdH=date_ymdH, 
            NHOURS=self.time_manager.run_hours)
        
        replace_nml_template(
            input_nml_path=self.paths.dart_run_filter_template(),
            entries_tbr_dict={
                "@NPROC": "1", #f"{self.nb_proc}",
                "@QUEUED_NODES": f"{self.queue}",
                "@PROJECT": f"{self.project_name}",
                "@WALLTIME": f"{self.walltime}",
                "@MAIL": f"{self.mail}",
                "@CURRENT_DATE": date_ymdH,
                "@DEST_LOG_PATH": self.paths.dart_posteriors_dir(date_ymdH),
                "@WORKDIR": self.paths.dart_run_filter().parent
            },
            output_nml_path=self.paths.dart_run_filter(),
        )
        spec = CommandSpec(
                command= f"./{self.paths.dart_run_filter().name}", #f"ccc_msub ./{self.paths.dart_run_filter().name}",
                directory=self.paths.dart_run_filter().parent
            )
        job_id = submit_irene(spec)
        if job_id:
            monitor_job_status([job_id], self.scheduler)
        self.move_analysis_files()

        logger.info(f"Computing differences between analysis/preassim means (ana - preassim)...")
        # 1. Difference: analysis_mean.nc - preassim_mean.nc
        analysis_mean = self.paths.dart_analysis_dir(date_ymdH) / "analysis_mean.nc"
        preassim_mean = self.paths.dart_preassim_dir(date_ymdH) / "preassim_mean.nc"
        diff_mean_out = self.paths.dart_analysis_dir(date_ymdH) / "analysis_increment_mean.nc"

        if analysis_mean.exists() and preassim_mean.exists():
            save_diff(analysis_mean, preassim_mean, diff_mean_out, "Mean Analysis Increment")
        else:
            raise FatalPipelineError(f"DART failed to produce analysis and/or preassim means for {date_ymdH}, cannot compute differences.")

        logger.info(f"run_dart() is DONE.")
    
    def move_analysis_files(self):
        date_ymdH = self.time_manager.end_file_datetime.strftime("%Y%m%d%H")
        self.paths.dart_analysis_dir(date_ymdH).mkdir(parents=True, exist_ok=True)
        self.paths.dart_preassim_dir(date_ymdH).mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Moving DART output files to analysis and preassim directories for date {date_ymdH} if present ...")
        for filename in os.listdir(f"{self.paths.path_filter}"):
            if filename.startswith("analysis_"):
                try:
                    shutil.move(
                        os.path.join(self.paths.path_filter, filename),
                        os.path.join(self.paths.dart_analysis_dir(date_ymdH), filename))
                    logger.info(f"Moved '{filename}' to '{self.paths.dart_analysis_dir(date_ymdH)}'")
                except shutil.Error:
                    logger.error(f"Failed to move '{filename}' to '{self.paths.dart_analysis_dir(date_ymdH)}' because it already exists.")
            elif filename.startswith("preassim_"):
                try:
                    shutil.move(
                        os.path.join(self.paths.path_filter, filename),
                        os.path.join(self.paths.dart_preassim_dir(date_ymdH), filename))
                    logger.info(f"Moved '{filename}' to '{self.paths.dart_preassim_dir(date_ymdH)}'")
                except shutil.Error:
                    logger.error(f"Failed to move '{filename}' to '{self.paths.dart_preassim_dir(date_ymdH)}' because it already exists.")


    def after_assimilation(self):
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled by config")
            return
        if not self.satdata_found:
            logger.info("after_assimilation() skipped")
            return
        if not self.update_restart:
            logger.info("No cycling flag was selected, skipping restart file update (Data Fusion)")
            return
        logger.info("---------->>> Running update_pollutant_in_end()")
        date_ymdH = self.time_manager.end_file_datetime.strftime("%Y%m%d%H")
        for mem in range(1, self.no_mems + 1):
            update_pollutant_in_end(dart_file=self.paths.dart_filter_output_list_file(mem, date_ymdH, self.time_manager.run_hours), 
                                    end_file=self.paths.chimere2023_END_FILE(mem, date_ymdH, self.time_manager.run_hours),
                                    out_file=self.paths.chimere2023_out_file(mem, date_ymdH, self.time_manager.run_hours),
                                    pollutant='NO2')
            logger.info(f"Computing differences between posterior vs. original CHIMERE outputs ...")
            # 2. Difference: Posterior (DART) - Original (CHIMERE) for each ensemble member
            posterior_file = self.paths.dart_filter_output_list_file(mem, date_ymdH, self.time_manager.run_hours)
            chimere_file = self.paths.chimere2023_out_file(mem, date_ymdH, self.time_manager.run_hours)
            # Saving in the posteriors directory for the specific date
            diff_mem_out = self.paths.dart_posteriors_dir(date_ymdH) / f"diff_posterior_ENS{mem}_{date_ymdH}.nc"

            if posterior_file.exists() and chimere_file.exists():
                save_diff(posterior_file, chimere_file, diff_mem_out, f"Posterior Diff ENS{mem}")
            else:
                logger.debug(f"Skipping member {mem} diff: files missing.")
            
    def finalize_step(self):
        """
        Cleanup + YAML update.
        """

        modify_yaml_date(
            self.config._config_path,
            self.time_manager.simulated_time.strftime("%Y-%m-%d %H:00:00"),
        )
        self.time_manager.prev_run_hours = self.time_manager.run_hours
        
        self.satdata_found = False
        logger.info("Cycle is DONE; starting a new loop!")


from pipelines.registry import register_pipeline
from pipelines.chimere2023.config import Chimere2023PipelineConfig


class Chimere2023Pipeline(Chimere2023DartPipeline):
    """Concrete CHIMERE 2023 pipeline with in-folder implementation."""


def _build(config, time_manager):
    p_cfg = Chimere2023PipelineConfig.from_app_config(config)
    p_paths = Chimere2023Paths(p_cfg)
    return Chimere2023Pipeline(time_manager, p_paths, config)


register_pipeline("chimere2023", _build)
