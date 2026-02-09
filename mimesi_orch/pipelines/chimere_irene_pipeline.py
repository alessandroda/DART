from datetime import timedelta
from pathlib import Path
import shutil
import subprocess
import time

from mimesi_types import Scheduler
from config_models import AppConfig
from paths import PathManager
from pipelines.base_pipeline import BaseAssimilationPipeline
from pipeline_errors import FatalPipelineError, ModelRunError
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
    get_list_mems_to_rerun,
    compute_hourly
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

        self.scheduler = self.config.cluster.scheduler
        logger.info(f"Using scheduler={self.scheduler}, queue={self.queue}")


    def before_step(self):
        logger.info(f">>> current time: {self.time_manager.current_time}; start: {self.time_manager.start_time}; end: {self.time_manager.end_time}")
        dict_mem_list = from_liststr_to_listdict(self.ensemble_list, self.perturbation_names)
        if len(dict_mem_list) != self.no_mems:
            logger.info(f"NB The ensemble list values ({len(dict_mem_list)}) are not equal in number to the requested assimilation ensemble ({self.no_mems})")

        for dict_mem in dict_mem_list: #loop on mems
            self.path_manager.chimere2023_run_dir(dict_mem["MemberID"]).mkdir(parents=True, exist_ok=True)
            logger.info(f"Creating directories and links for ENS{dict_mem['MemberID']} to run chimere's parallel part")
            #link WPS
            #TODO simply replace template (move in run_model)
            safe_symlink(self.path_manager.chimere2023_WPS(), 
                         self.path_manager.chimere2023_run_dir_WPS(dict_mem["MemberID"]))
            #link BOUN, EMIS and METEO (from either control run or perturbed dir))
            for date in pd.date_range(self.time_manager.start_time, self.time_manager.end_time, freq='D'):
                date_ymd = date.strftime("%Y%m%d")
                date_ymd00 = date.strftime("%Y%m%d00")
                date_month = date.strftime("%m")
                date_weekday = date.strftime("%A")
                safe_symlink(self.path_manager.chimere2023_EMIS_FILE_SRC("EmisID" in dict_mem.keys(), self.domain, date_month, date_weekday, dict_mem["EmisID"]),
                             self.path_manager.chimere2023_EMI_FILE(dict_mem["MemberID"], self.domain, date_month, date_weekday))
                safe_symlink(self.path_manager.chimere2023_METEO_FILE_SRC("MeteoID" in dict_mem.keys(), self.domain, date_ymd, dict_mem["MeteoID"]), 
                             self.path_manager.chimere2023_METEO_FILE(dict_mem["MemberID"], self.domain, date_ymd00, 24))
                safe_symlink(self.path_manager.chimere2023_BOUN_LIST_SRC(date_ymd, self.control_run_exp_name, self.domain), 
                             self.path_manager.chimere2023_BOUN_LIST(date_ymd00, 24, dict_mem["MemberID"], self.domain))
            #link first end file (if it is a continuation run, the end file in the folder is kept)
            safe_symlink(self.path_manager.chimere2023_END_FILE_SRC(self.control_run_exp_name, self.time_manager.end_file_date_control_run.strftime("%Y%m%d")), 
                         self.path_manager.chimere2023_END_FILE(dict_mem["MemberID"], self.time_manager.end_file_date_control_run.strftime("%Y%m%d")))
            #link anche a emis del giorno succ se ci si ferma a 00 e non ma le si vogliono
            
            if check_and_clean_broken_links(self.path_manager.chimere2023_run_dir(dict_mem["MemberID"])):
                raise FatalPipelineError(f"Broken links detected. Please clean up. The target needs to exists.")
            else:
                logger.info(f">> All links are good for ENS{dict_mem['MemberID']}  ...")

    def run_model(self):
        """
        Run CHIMERE for the current time step.
        """
        logger.info(f"[STEP] Running CHIMERE model at {self.time_manager.current_time}")
        
        date_ymd = self.time_manager.current_time.strftime("%Y%m%d")
        date_ymdH = self.time_manager.current_time.strftime("%Y%m%d%H")
        date_H = self.time_manager.current_time.strftime("%H")
        self.time_manager.end_file_date = self.time_manager.current_time - timedelta(days=1)
        end_file_date_ymd = self.time_manager.end_file_date.strftime("%Y%m%d")

        job_ids = []
        for mem in range(self.no_mems):
            logger.info("Computing BOUNs for the specific hours to run")
            #sarà da cambiare quando gireremo su piu ore
            compute_hourly(self.path_manager.chimere2023_BOUN_LIST(date_ymdH, 24, mem, self.domain).read_text().splitlines()[1], 
                                   int(date_H), 
                                   self.path_manager.chimere2023_BOUN_NC(date_ymdH, 1, mem, self.control_run_exp_name, self.domain), 
                                   self.path_manager.chimere2023_BOUN_LIST(date_ymdH, 1, mem, self.domain))
            logger.info("Computing exdomouts for the specific hours to run")
            #sarà da cambiare quando gireremo su piu ore
            compute_hourly(str(self.path_manager.chimere2023_METEO_FILE(mem, self.domain, date_ymdH, 24)), 
                                   int(date_H), 
                                   self.path_manager.chimere2023_METEO_FILE(mem, self.domain, date_ymdH, 1))
            try:
                logger.info("Replacing @TOKENS in CHIMERE .par template file ...")
                logger.info(f"The output directory (run_dir) is: {self.path_manager.chimere2023_run_dir(mem)}")
                replace_nml_template(
                    input_nml_path=self.path_manager.chimere2023_PAR_BASE_TEMPLATE(),
                    entries_tbr_dict={
                        "@LAB": f"ENS{mem}",
                        "@SIMULDIR": self.path_manager.chimere2023_run_dir(mem),
                        "@IUSEINI": "2",
                        "@ENDFILE": self.path_manager.chimere2023_END_FILE(mem, end_file_date_ymd), #per ora c'é 24 dentro 
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
            if len(job_ids)==0 or not self.config.model_data.submit_sequentially:
                submit_command = f"ccc_msub ./{self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem).name}"
            else:
                submit_command = f"ccc_msub -a {job_ids[-1]} ./{self.path_manager.chimere2023_BASH_SUBMIT_SCRIPT(mem).name}"

            job_ids.append(submit_irene(CommandSpec(
                command=submit_command,
                directory=self.path_manager.base_path
            )))
            
        logger.info(f"Checking job status ...")
        mems_to_rerun = get_list_mems_to_rerun(
            job_ids,
            self.scheduler,
            self.model_type,
            )
        #if not mems_to_rerun:
        #    raise ModelRunError("submission stopped")
        if mems_to_rerun:
            raise ModelRunError(f"Ensemble members failed: {mems_to_rerun}")
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

        orbit_filename = self.process_satellite_data()
        if not orbit_filename:
            logger.info("[DART] No satellite data found, skipping assimilation")
            return

        obs_seq_name = self.run_obs_converter(orbit_filename)

        obs_path = (
            self.path_manager.base_path
            / f"DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/{obs_seq_name}"
        )
        if not obs_path.exists():
            logger.info("[DART] obs_seq not found, skipping assimilation")
            return

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
            logger.info(
                f'No valid orbit file found after 10 AM.{orbit_filename["start_time"]}'
            )
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
            self.path_manager.dart_s5p_input_template(),
            entries_tbr_dict={
                "$file_path_s5p": self.path_manager.dart_file_s5p_orbit(orbit_filename),
                "$file_out": self.path_manager.dart_obs_seq(
                    self.seconds_obs, self.days_obs
                ),
                "$obs_type": self.obs_type,
            },
            output_nml_path=self.path_manager.dart_s5p_input(),
        )
        try:
            run_command_in_directory(
                "convert_s5p_tropomi_l3",
                self.path_manager.dart_s5p_work(),
            )
        except Exception as e:
            logger.error(f"Error running obs converter: {e}")
            return False
        return obs_seq_name

    def run_dart(self, obs_seq_name):
        logger.info("Running DART")

        self.output_sim_folder = self.path_manager.dart_posteriors_dir(
            self.time_manager.simulated_time.strftime("%Y%m%d%H")
        )
        Path(self.output_sim_folder).mkdir(parents=True, exist_ok=True)

        replace_nml_template(
            self.path_manager.base_path
            / self.path_manager.path_filter
            / "input_template.nml",
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": self.output_sim_folder,
                "$folder_obs_path": self.path_manager.dart_s5p_output_dir(),
                "$date_assim": self.time_manager.current_time.strftime("%Y%m%d_%H%M%S"),
                "$template_farm": self.path_manager.path_data
                / f"to_DART/ic_g1_{self.seconds_model}_{self.days_model}_0.nc",
                "$init_time_days": str(self.days_model),
                "$init_time_seconds": str(self.seconds_model),
                "$first_obs_days": str(self.days_obs),
                "$first_obs_seconds": str(self.seconds_obs),
                "$no_mems": str(self.no_mems),
                "$obs_type": str(self.obs_type),
                "$state_variable_conc": str(self.state_variable_conc),
                "$state_variable_qty": str(self.state_variable_qty),
            },
            output_nml_path=self.path_manager.base_path
            / self.path_manager.path_filter
            / "input.nml",
        )
        # FILTER_INPUT_LIST.TXT
        replace_nml_template(
            self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_input_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.path_manager.path_data / f"to_DART/",
                "$days": str(self.seconds_model),
                "$seconds": str(self.days_model),
            },
            output_nml_path=self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_input_list.txt",
        )

        # FILTER_OUTPUT_LIST.TXT
        replace_nml_template(
            self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_output_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.output_sim_folder,
                "$date": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
            },
            output_nml_path=self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_output_list.txt",
        )
        # SUBMIT_FILTER.BSH
        replace_nml_template(
            self.path_manager.base_path
            / "RUN/script/templates/submit_filter.template.bsh",
            entries_tbr_dict={
                "CURRENT_DATE": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                "CORES": str(20),
                "QUEUE": self.cresco_queue,
                "DEST_LOG_PATH": self.output_sim_folder,
            },
            output_nml_path=self.path_manager.path_submit_bsh / "submit_filter.bsh",
        )
        # RUN_FILTER.BSH
        replace_nml_template(
            self.path_manager.base_path
            / "RUN/script/templates/run_filter.template.bsh",
            entries_tbr_dict={
                "CORES": str(20),
                "@ABS_FILTER_PATH": self.path_manager.base_path
                / self.path_manager.path_filter,
            },
            output_nml_path=self.path_manager.path_submit_bsh / "run_filter.bsh",
        )

        job_id = run_command_in_directory_bsub( #da portare a run_command_in_directory
            "./submit_filter.bsh", self.path_manager.path_submit_bsh, farm=False
        )
        time.sleep(10)
        self.monitor_job_dart(job_id)

    def monitor_job_dart(self, job_id):
        logger.info(f"Monitoring job {job_id}")
        job_id = job_id.strip()[1:-1]

        while True:
            if check_job_status_cresco(job_id, which_run="FARM"):
                print("Job completed successfully.")
                # Handle successful job completion: move files
                self.move_analysis_files()
                replace_priorinflation(
                    self.path_manager,
                    self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                )
                break
            else:
                print("Job is still running. Waiting...")
                time.sleep(10)

    def move_analysis_files(self):
        analysis_sim_folder = (
            self.path_manager.path_data
            / f"analysis/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
        preassim_sim_folder = (
            self.path_manager.path_data
            / f"preassim/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)
        for filename in os.listdir(f"{self.path_manager.path_filter}"):
            if filename.startswith("analysis_"):
                try:
                    shutil.move(
                        os.path.join(
                            self.path_manager.path_filter,
                            filename,
                        ),
                        os.path.join(analysis_sim_folder, filename),
                    )
                except shutil.Error:
                    print(
                        f"Failed to move '{filename}' to '{analysis_sim_folder}' because it already exists."
                    )
            elif filename.startswith("preassim_"):
                try:
                    shutil.move(
                        os.path.join(
                            self.path_manager.path_filter,
                            filename,
                        ),
                        os.path.join(preassim_sim_folder, filename),
                    )
                except shutil.Error:
                    print(
                        f"Failed to move '{filename}' to '{preassim_sim_folder}' because it already exists."
                    )

    
    def after_assimilation(self):
        """
        Optional hook executed after the assimilation step.

        Typical use cases:
        - mapping analysis fields back to the model format
        - updating boundary or initial conditions
        """
        pass

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
