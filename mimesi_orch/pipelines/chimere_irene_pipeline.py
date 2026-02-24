from datetime import timedelta
from pathlib import Path
import shutil
import subprocess
import time

from mimesi_types import Scheduler
from config_models import AppConfig
from paths import PathManager
from pipelines.base_pipeline import BaseAssimilationPipeline
from pipeline_errors import FatalPipelineError, ModelRunError, SchedulerError
import logging
import os
import pandas as pd
from orchestrator_utils import (
    CommandSpec,
    check_job_status_cresco,
    modify_yaml_date,
    filter_dates,
    TimeManager,
    replace_nml_template,
    replace_priorinflation,
    searchFile,
    set_date_gregorian,
    submit_and_wait,
    run_command_in_directory,
    run_command_in_directory_bsub,
    submit_and_wait_slurm,
    safe_symlink,
    check_and_clean_broken_links,
    from_liststr_to_listdict,
    submit_irene,
    monitor_job_status,
    check_restart_files_exist,
    compute_hourly,
    monitor_job_dart,
    add_missing_variable,
    write_dart_filter_list,
    update_pollutant_in_end
)


logger = logging.getLogger(__name__)


class ChimereV2023DartPipeline(BaseAssimilationPipeline):
    def __init__(
        self,
        time_manager: TimeManager,
        path_manager: PathManager,
        config: AppConfig,
    ):
        super().__init__(time_manager)
        # Load configuration from YAML file
        self.config = config
        # Setup path manager and time manager with loaded configurations
        self.path_manager = path_manager
        self.listing = pd.read_csv(self.path_manager.listing_file, sep=";")
        self.listing["start_time"] = pd.to_datetime(self.listing["start_time"])
        self.days_obs = 0
        self.seconds_obs = 0
        self.days_model = 0
        self.seconds_model = 0
        self.output_sim_folder = None
        self.current_end_file = {}
        a = self.config.assimilation

        self.model_type = a.model_type
        logger.info(f"Running assimilation with model_type={self.model_type}")

        self.ass_var = a.ass_var
        self.emi_var = a.emi_var
        self.no_mems = a.no_mems
        self.case_dir = a.case_dir
        self.obs_type = a.obs_type
        self.state_variable_conc = a.ass_var
        self.state_variable_qty = a.state_variable_qty
        self.run_assimilation_flag = a.run_assimilation_flag
        self.case_emi_dir = a.case_emi_dir

        self.queue = self.config.cluster.cluster_queue
        self.project_name = self.config.cluster.project_name
        self.nb_proc = self.config.cluster.nb_proc
        self.walltime = self.config.cluster.walltime
        self.mail = self.config.cluster.mail

        self.backup_perturb_days = self.config.time.backup_perturb_days
        self.backup_ic_hours = self.config.time.backup_ic_hours
        self.backup_ic_option = self.config.time.backup_ic_option

        self.ensemble_list = self.config.model_data.ensemble_list
        self.perturbation_names = self.config.model_data.perturbation_names
        self.control_run_exp_name = self.config.model_data.control_run_exp_name
        self.domain = self.config.model_data.domain
        self.chimpart = self.config.model_data.chimpart
        self.submit_sequentially = self.config.model_data.submit_sequentially
        self.restart_from_controlrun = self.config.model_data.restart_from_controlrun
        self.dom_west = self.config.model_data.dom_west
        self.dom_east = self.config.model_data.dom_east
        self.dom_south = self.config.model_data.dom_south
        self.dom_north = self.config.model_data.dom_north
        self.nz = self.config.model_data.nz
        self.dlon = self.config.model_data.dlon
        self.dlat = self.config.model_data.dlat

        self.search_window_seconds = self.config.satellite_data.search_window_seconds
        self.obs_name = self.config.satellite_data.obs_name
        self.collection = self.config.satellite_data.collection
        self.vertical_ref_height = self.config.satellite_data.vertical_ref_height
        self.superobs = self.config.satellite_data.superobs

        self.scheduler = self.config.cluster.scheduler
        logger.info(f"Using scheduler={self.scheduler}, queue={self.queue}")

    def before_step(self):
        logger.info(f"---------->>> current time: {self.time_manager.current_time}; start: {self.time_manager.start_time}; end: {self.time_manager.end_time}")
        dict_mem_list = from_liststr_to_listdict(self.ensemble_list, self.perturbation_names)
        if len(dict_mem_list) != self.no_mems:
            logger.info(f"NB The ensemble list values ({len(dict_mem_list)}) are not equal in number to the requested assimilation ensemble ({self.no_mems})")

        for dict_mem in dict_mem_list: #loop on mems
            self.path_manager.chimere2023_run_dir(dict_mem["MemberID"]).mkdir(parents=True, exist_ok=True)
            logger.info(f"Creating directories and links for ENS{dict_mem['MemberID']} to run chimere's parallel part")
            #link WPS
            #TODO simply replace template (move in run_model)
            #logger.info("Linking WPS ...")
            #safe_symlink(self.path_manager.chimere2023_WPS(), 
            #             self.path_manager.chimere2023_run_dir_WPS(dict_mem["MemberID"]))
            #link BOUN, EMIS and METEO (from either control run or perturbed dir))
            date_month = self.time_manager.current_time.strftime("%m")
            date_weekday = self.time_manager.current_time.strftime("%A")
            logger.info("Linking EMIS ...")
            safe_symlink(self.path_manager.chimere2023_EMIS_FILE_SRC("EmisID" in dict_mem.keys(), self.domain, date_month, date_weekday, dict_mem["EmisID"]),
                            self.path_manager.chimere2023_EMI_FILE(dict_mem["MemberID"], self.domain, date_month, date_weekday))
            if self.time_manager.current_time.strftime("%H") == '23':
                safe_symlink(self.path_manager.chimere2023_EMIS_FILE_SRC("EmisID" in dict_mem.keys(), self.domain, date_month, (self.time_manager.current_time + timedelta(days=1)).strftime("%A"), dict_mem["EmisID"]),
                            self.path_manager.chimere2023_EMI_FILE(dict_mem["MemberID"], self.domain, date_month, (self.time_manager.current_time + timedelta(days=1)).strftime("%A")))
            if self.restart_from_controlrun and self.time_manager.current_time == self.time_manager.start_time:
                #link first end file
                logger.info("Linking END ...")
                safe_symlink(self.path_manager.chimere2023_END_FILE_SRC(self.control_run_exp_name, self.time_manager.end_file_date_control_run.strftime("%Y%m%d")), 
                             self.path_manager.chimere2023_END_FILE(dict_mem["MemberID"], self.time_manager.end_file_date_control_run.strftime("%Y%m%d00"), 24))
                self.current_end_file[dict_mem["MemberID"]]=self.path_manager.chimere2023_END_FILE(dict_mem["MemberID"], self.time_manager.end_file_date_control_run.strftime("%Y%m%d00"), 24)
            else:
                self.current_end_file[dict_mem["MemberID"]]=None
            if check_and_clean_broken_links(self.path_manager.chimere2023_run_dir(dict_mem["MemberID"])):
                raise FatalPipelineError(f"Broken links detected. Please clean up. The target needs to exists.")
            else:
                logger.info(f">> All links are good for ENS{dict_mem['MemberID']}  ...")
            
            date_ymd = self.time_manager.current_time.strftime("%Y%m%d")
            date_ymdH = self.time_manager.current_time.strftime("%Y%m%d%H")
            date_H = self.time_manager.current_time.strftime("%H")
            #date_ymd00 = self.time_manager.current_time.strftime("%Y%m%d00")
            #logger.info("Linking METEO ...")
            #safe_symlink(self.path_manager.chimere2023_METEO_FILE_SRC("MeteoID" in dict_mem.keys(), self.domain, date_ymd, dict_mem["MeteoID"]), 
            #                self.path_manager.chimere2023_METEO_FILE(dict_mem["MemberID"], self.domain, date_ymd00, 24))
            #logger.info("Linking BOUN ...")
            #safe_symlink(self.path_manager.chimere2023_BOUN_LIST_SRC(date_ymd, self.control_run_exp_name, self.domain), 
            #                self.path_manager.chimere2023_BOUN_LIST(date_ymd00, 24, dict_mem["MemberID"], self.domain))
            
            logger.info("Computing BOUNs for the specific hours to run")
            #sarà da cambiare quando gireremo su piu ore
            compute_hourly(self.path_manager.chimere2023_BOUN_LIST_SRC(date_ymd, self.control_run_exp_name, self.domain).read_text().splitlines()[1], 
                                   int(date_H), 
                                   self.path_manager.chimere2023_BOUN_NC(date_ymdH, 1, dict_mem["MemberID"], self.control_run_exp_name, self.domain), 
                                   self.path_manager.chimere2023_BOUN_LIST(date_ymdH, 1, dict_mem["MemberID"], self.domain))
            logger.info("Computing exdomouts for the specific hours to run")
            #sarà da cambiare quando gireremo su piu ore
            compute_hourly(str(self.path_manager.chimere2023_METEO_FILE_SRC("MeteoID" in dict_mem.keys(), self.domain, date_ymd, dict_mem["MeteoID"])), 
                                   int(date_H), 
                                   self.path_manager.chimere2023_METEO_FILE(dict_mem["MemberID"], self.domain, date_ymdH, 1))

    def run_model(self):
        """
        Run CHIMERE for the current time step.
        """
        logger.info(f"---------->>> Running CHIMERE model at {self.time_manager.current_time}")
    
        date_ymdH = self.time_manager.current_time.strftime("%Y%m%d%H")
        self.time_manager.end_file_date = self.time_manager.current_time - timedelta(hours=1)
        end_file_date_ymdH = self.time_manager.end_file_date.strftime("%Y%m%d%H")

        job_ids = []
        for mem in range(self.no_mems):
            try:
                logger.info("Replacing @TOKENS in CHIMERE .par template file ...")
                logger.info(f"The output directory (run_dir) is: {self.path_manager.chimere2023_run_dir(mem)}")
                replace_nml_template(
                    input_nml_path=self.path_manager.chimere2023_PAR_BASE_TEMPLATE(),
                    entries_tbr_dict={
                        "@LAB": f"ENS{mem}",
                        "@SIMULDIR": self.path_manager.chimere2023_run_dir(mem),
                        "@IUSEINI": "2",
                        "@ENDFILE": self.path_manager.chimere2023_END_FILE(mem, end_file_date_ymdH, 1) if self.current_end_file[mem] is None else self.current_end_file[mem],
                        "@EMISSDIR": self.path_manager.chimere2023_run_dir(mem)
                    },
                    output_nml_path=self.path_manager.chimere2023_PAR_FILE(mem)
                )
                shutil.copy(self.path_manager.chimere2023_PAR_FILE(mem), self.path_manager.chimere2023_PAR_FILE_RUN_DIR(mem, date_ymdH))
            except Exception as e:
                raise FatalPipelineError(f"Failed to prepare CHIMERE .par file: {e}")
            try:
                logger.info("Replacing @TOKENS in CHIMERE template submit script ...")
                replace_nml_template(
                    input_nml_path=self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT_TEMPLATE(),
                    entries_tbr_dict={
                        "@JOBNAME": f"ENS{mem}",
                        "@NPROC": f"{self.nb_proc}",
                        "@QUEUED_NODES": f"{self.queue}",
                        "@PROJECT": f"{self.project_name}",
                        "@WALLTIME": f"{self.walltime}",
                        "@RUN_DIR": f"{self.path_manager.chimere2023_run_dir(mem)}",
                        "@MAIL": f"{self.mail}",
                        "@PARFILE": f"{self.path_manager.chimere2023_PAR_FILE(mem).name}",
                        "@START_DATEHOUR": f"{date_ymdH}",
                        "@NHOURS_forward": "1",
                        "@CHIMPART": f"{self.chimpart}"
                    },
                    output_nml_path=self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem)
                )
                shutil.copy(self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem), self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT_RUN_DIR(mem, date_ymdH))
                subprocess.run(["chmod", "+x", str(self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem))], check=True)
            except Exception as e:
                raise FatalPipelineError(f"Failed to prepare CHIMERE submit script: {e}")
            
            logger.info(f"Queuing job for member {mem}...")
            if len(job_ids)==0 or not self.submit_sequentially:
                submit_command = f"ccc_msub ./{self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem).name}"
            else:
                submit_command = f"ccc_msub -a {job_ids[-1]} ./{self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem).name}"

            job_ids.append(submit_irene(CommandSpec(
                command=submit_command,
                directory=self.path_manager.base_path_ctm
            )))
            
        logger.info(f"Checking job status ...")
        monitor_job_status(job_ids, self.scheduler, self.model_type)
        mems_to_rerun = check_restart_files_exist(
            ic_path = self.path_manager.chimere2023_END_FILE(mem, date_ymdH, 1),
            model=self.model_type,
            no_mems=self.no_mems,
            )
        if mems_to_rerun:
            logger.info(f"Check chimere log file at: {self.path_manager.chimere2023_run_dir(mem)}/ENS{mem}_{date_ymdH}.out")
            raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
        logger.info(f" Run_model() compled successfully.")

    def after_model(self):
        pass

    def run_assimilation_if_needed(self):
        """
        Run DART assimilation only if:
        - assimilation is enabled
        - satellite observations are available
        """
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled by config")
            return
        logger.info(f"---------->>> Running process_satellite_data() for {self.time_manager.current_time}")
        orbit_filename = self.process_satellite_data()
        if not orbit_filename:
            logger.info("[DART] No satellite data found, skipping assimilation")
            return

        logger.info(f"---------->>> Running run_obs_converter() for {self.time_manager.current_time}")
        obs_seq_name = self.run_obs_converter(orbit_filename)

        obs_path = (self.path_manager.dart_s5p_output_dir(self.obs_name, self.collection) / Path(obs_seq_name))
        if not obs_path.exists():
            logger.info("[DART] obs_seq not found, skipping assimilation")
            return
        else:
            logger.info(f"[DART] obs_seq created: {obs_path}")
        
        self.run_dart(obs_seq_name)

    def process_satellite_data(self):

        orbit_filename = searchFile(
            self.time_manager.current_time - 0.5 * self.time_manager.dt,
            self.time_manager.current_time + 0.5 * self.time_manager.dt,
            self.listing,
        )
        if orbit_filename.empty:
            return False

        orbit_filename["start_time"] = pd.to_datetime(orbit_filename["start_time"])
        orbit_filename = orbit_filename[orbit_filename["start_time"].dt.hour >= 10]

        if orbit_filename.empty:
            logger.info(f'No valid orbit file found after 10 AM.{orbit_filename["start_time"]}')
            return False

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

        replace_nml_template(
            input_nml_path=self.path_manager.dart_s5p_input_template(),
            entries_tbr_dict={
                "$file_path_s5p": self.path_manager.dart_file_s5p_orbit(orbit_filename, self.obs_name),
                "$file_out": self.path_manager.dart_obs_seq(
                    self.seconds_obs, self.days_obs, self.obs_name, self.collection
                ),
                "$obs_type": self.obs_type,
                "$dom_west": self.dom_west,
                "$dom_east": self.dom_east,
                "$dom_south": self.dom_south,
                "$dom_north": self.dom_north,
                "$nz": self.nz,
                "$dlon": self.dlon,
                "$dlat": self.dlat,
                "$vertical_ref_height": self.vertical_ref_height,
                "$superobs": self.superobs
            },
            output_nml_path=self.path_manager.dart_s5p_input(),
        )
        try:
            spec = CommandSpec(
                command="convert_s5p_tropomi_l3",
                directory=self.path_manager.dart_s5p_work(),
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

    def run_dart(self, obs_seq_name):
        logger.info("---------->>> Running DART")
        """
        self.time_manager.end_file_date = self.time_manager.current_time - timedelta(hours=1)
        end_file_date_ymdH = self.time_manager.end_file_date.strftime("%Y%m%d%H")

        create:
        @property
        def end_file_date(self):
            return self.current_time - timedelta(hours=1)

        because all chimere files in the title have the starting time ==> change the date_ymdH used in the following
        """
        
        date_ymdH = self.time_manager.simulated_time.strftime("%Y%m%d%H")
        date_ymdHMS = self.time_manager.simulated_time.strftime("%Y%m%d_%H%M%S")
        self.path_manager.dart_posteriors_dir(date_ymdH).mkdir(parents=True, exist_ok=True)

        add_missing_variable(
            no_mems=self.no_mems, 
            var_to_add='psfc', 
            domain=self.domain, 
            out_file_func=self.path_manager.chimere2023_out_file,
            orig_file_func=self.path_manager.chimere2023_METEO_FILE,
            date_ymdH=date_ymdH, 
            NHOURS=1
        )

        replace_nml_template(
            input_nml_path=self.path_manager.dart_filter_input_template(),
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": self.path_manager.dart_posteriors_dir(date_ymdH),
                "$folder_obs_path": self.path_manager.dart_s5p_output_dir(self.obs_name, self.collection),
                "$date_assim": date_ymdHMS,
                "$init_time_days": str(self.days_model), #computed in base_pipeline
                "$init_time_seconds": str(self.seconds_model), #computed in base_pipeline
                "$first_obs_days": str(self.days_obs), #computed when creating obs_seq.out 
                "$first_obs_seconds": str(self.seconds_obs), #computed when creating obs_seq.out
                "$no_mems": str(self.no_mems),
                "$obs_type": str(self.obs_type),
            },
            output_nml_path=self.path_manager.dart_filter_input(),
        )

        write_dart_filter_list(
            list_file_func=self.path_manager.dart_filter_input_list(), 
            out_file_func=self.path_manager.chimere2023_out_file, 
            no_mems=self.no_mems, 
            date_ymdH=date_ymdH, 
            NHOURS=1)

        write_dart_filter_list(
            list_file_func=self.path_manager.dart_filter_output_list(), 
            out_file_func=self.path_manager.dart_filter_output_list_file, 
            no_mems=self.no_mems, 
            date_ymdH=date_ymdH, 
            NHOURS=1)
        
        replace_nml_template(
            input_nml_path=self.path_manager.dart_run_filter_template(),
            entries_tbr_dict={
                "@NPROC": f"{self.nb_proc}",
                "@QUEUED_NODES": f"{self.queue}",
                "@PROJECT": f"{self.project_name}",
                "@WALLTIME": f"{self.walltime}",
                "@MAIL": f"{self.mail}",
                "@CURRENT_DATE": date_ymdH,
                "@DEST_LOG_PATH": self.path_manager.dart_posteriors_dir(date_ymdH),
            },
            output_nml_path=self.path_manager.dart_run_filter(),
        )
        spec = CommandSpec(
                command=f"ccc_msub ./{self.path_manager.dart_run_filter().name}",
                directory=self.path_manager.dart_run_filter().parent
            )
        job_id = submit_irene(spec)
        monitor_job_status([job_id], self.scheduler, self.model_type)
        ### check filter succeded

    def after_assimilation(self):
        """
        Optional hook executed after the assimilation step.

        Typical use cases:
        - mapping analysis fields back to the model format
        - updating boundary or initial conditions
        """
        logger.info("---------->>> Running after_assimilation()")
        self.move_analysis_files()
        logger.info("---------->>> Running update_pollutant_in_end()")

        #finish the following and test 
        """date_ymdH = self.time_manager.simulated_time.strftime("%Y%m%d%H")
        for mem in range(self.no_mems):
            update_pollutant_in_end(self.path_manager.dart_filter_output_list_file(mem, date_ymdH, 1), self.path_manager.chimere2023_END_FILE(mem, date_ymdH, 1))
        """
    
    def move_analysis_files(self):
        date_ymdH = self.time_manager.simulated_time.strftime("%Y%m%d%H")
        self.path_manager.dart_analysis_dir(date_ymdH).mkdir(parents=True, exist_ok=True)
        self.path_manager.dart_preassim_dir(date_ymdH).mkdir(parents=True, exist_ok=True)
        
        for filename in os.listdir(f"{self.path_manager.path_filter}"):
            if filename.startswith("analysis_"):
                try:
                    shutil.move(
                        os.path.join(self.path_manager.path_filter, filename),
                        os.path.join(self.path_manager.dart_analysis_dir(date_ymdH), filename))
                    logger.info(f"Moved '{filename}' to '{self.path_manager.dart_analysis_dir(date_ymdH)}'")
                except shutil.Error:
                    logger.error(f"Failed to move '{filename}' to '{self.path_manager.dart_analysis_dir(date_ymdH)}' because it already exists.")
            elif filename.startswith("preassim_"):
                try:
                    shutil.move(
                        os.path.join(self.path_manager.path_filter, filename),
                        os.path.join(self.path_manager.dart_preassim_dir(date_ymdH), filename))
                    logger.error(f"Moved '{filename}' to '{self.path_manager.dart_preassim_dir(date_ymdH)}'")
                except shutil.Error:
                    logger.error(f"Failed to move '{filename}' to '{self.path_manager.dart_preassim_dir(date_ymdH)}' because it already exists.")

    def finalize_step(self):
        """
        Cleanup + YAML update.
        """

        #modify_yaml_date(
        #    self.config["_config_path"],
        #    self.time_manager.simulated_time.strftime("%Y-%m-%d %H:00:00"),
        #)
        #self.cleanup_FARM()

        pass

